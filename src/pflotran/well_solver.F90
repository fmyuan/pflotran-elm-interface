module Well_Solver_module

!  This well solver finds the current state of a well, automatically switching
!  to the controlling target mode.
!  The calculated well flows are place in the Pflotran residual
!  Well data is obtained from     the Well_Data class
!  Solution results are stored in the Well_Data class

#include "petsc/finclude/petscsys.h"

  use petscsys
  use PFLOTRAN_Constants_module
  use Well_Type_class
  use Well_Data_class
  use PM_TOWG_Aux_module
  use Option_module

  implicit none

  private

    PetscInt, parameter :: max_iterTT  = 20  ! Max iterations over target types
    PetscInt, parameter :: max_iterPW  = 20  ! Max iterations for well solution
    PetscInt, parameter :: max_iterWBS = 20  ! Max iterations for well solution

    PetscReal, parameter :: ws_atm   = 1.01325d5  ! Used for default bhp limits
    PetscReal, parameter :: ws_pwmin = 0.1*ws_atm

    PetscBool :: ws_is_black_oil = PETSC_FALSE ! Is black oil with bubble point

    PetscReal :: ws_gravity                    ! Vertical gravity constant

    PetscBool :: ws_isothermal = PETSC_TRUE    ! Case is isothermal

    ! Bools indicating phases are present

    PetscBool :: ws_oil
    PetscBool :: ws_gas
    PetscBool :: ws_wat
    PetscBool :: ws_slv

    ! Pointers into well solution

    PetscInt :: ws_loc_soil = -1
    PetscInt :: ws_loc_sgpb = -1
    PetscInt :: ws_loc_trel = -1
    PetscInt :: ws_loc_sslv = -1

    PetscInt :: ws_nphase             ! Number of phases
    PetscInt :: ws_ncomp              ! Number of components
    PetscInt :: ws_ncompe             ! Number of components plus energy
    PetscInt :: ws_nxw                ! Number of variables in well bore soln
    PetscInt :: ws_ndof               ! Number of primary variables/cell

    PetscInt :: w_type                ! Well type (producer, gas injector...)

    PetscBool :: ws_isproducer    = PETSC_TRUE  ! Flags for prod. or inje.
    PetscBool :: ws_isinjector    = PETSC_FALSE

    PetscBool :: ws_isoilinjector = PETSC_FALSE ! Flags for injector type
    PetscBool :: ws_isgasinjector = PETSC_FALSE
    PetscBool :: ws_iswatinjector = PETSC_FALSE
    PetscBool :: ws_isslvinjector = PETSC_FALSE

    PetscInt :: w_ncmpl  = 0             ! Number of completions (this proc)
    PetscInt :: w_ncmplg = 0             ! Number of completions (all procs)

    PetscInt :: ws_xpr = 1               ! Location of pressure in ndof

    PetscBool :: w_crossproc = PETSC_FALSE ! Flag indicating cross-proc well

    PetscReal :: w_sign                    ! Sign (+ve prod, -ve inje)

    PetscReal :: w_trel                    ! Well rel. temp. (deg C)
    PetscReal :: w_pb                      ! Well bubble point pressure
    PetscBool :: w_issat                   ! Well state flag

    PetscReal, allocatable :: w_sp(:)      ! Well saturations by phase
    PetscReal, allocatable :: w_spxw(:,:)  ! ...wrt well soln

    ! Well molar densities (moles/unit res. vol) by phase
    PetscReal, allocatable :: w_mdp(:)
    PetscReal, allocatable :: w_mdppw(:)   ! ...wrt well pressure
    PetscReal, allocatable :: w_mdpxw(:,:) ! ...wrt well soln

    ! Well gravity densities (G.Kg/unit res. vol) by phase
    PetscReal, allocatable :: w_kgdp(:)     ! Well grav dens by phase
    PetscReal, allocatable :: w_kgdppw(:)   ! ...wrt well pressure
    PetscReal, allocatable :: w_kgdpxw(:,:) ! ...wrt well soln

    ! Well component molar densities (moles/unit res. vol.) by component
    PetscReal, allocatable :: w_mdc(:)      ! Well comp. mol. dens.
    PetscReal, allocatable :: w_mdcpw(:)    ! ...wrt well pressure
    PetscReal, allocatable :: w_mdcxw(:,:)  ! ...wrt well soln

    ! Well component mole fractions by component
    PetscReal, allocatable :: w_zmf(:)      ! Well mole fractions
    PetscReal, allocatable :: w_zmfpw(:)    ! ...wrt well pressure
    PetscReal, allocatable :: w_zmfxw(:,:)  ! ...wrt well soln

    ! Wellbore phase enthalpy/unit res vol
    PetscReal, allocatable :: w_hp(:)
    PetscReal, allocatable :: w_hpdpw(:)    ! ...wrt well pressure
    PetscReal, allocatable :: w_hpdxw(:,:)  ! ...wrt well soln

    ! Wellbore total enthalpy/unit res vol

    PetscReal :: w_h
    PetscReal :: w_hdpw                 ! ...wrt well pressure
    PetscReal, allocatable :: w_hdxw(:) ! ...wrt well soln

    ! Wellbore total enthalpy/mole

    PetscReal :: w_hpm
    PetscReal :: w_hpmdpw                 ! ...wrt well pressure
    PetscReal, allocatable :: w_hpmdxw(:) ! ...wrt well soln

    PetscReal :: w_zref                  ! Reference elevation for well

    ! Gravity density (G*well mass density)

    PetscReal :: w_gd                   ! Grav. den.
    PetscReal :: w_gdpw                 ! Grav. den. wrt Pw
    PetscReal, allocatable :: w_gdxw(:) ! Grav. den. wrt Xw

    PetscReal :: ws_targets(N_WELL_TT)          ! Well targets by target type
    PetscReal :: ws_actuals(N_WELL_TT)          ! Well actuals by target type

    PetscReal, allocatable :: ws_svpm(:) ! Well sv/mole   by comp
    PetscReal, allocatable :: ws_mspm(:) ! Well mass/mole by comp

    PetscInt , allocatable :: c_ghosted_id(:) ! Cmpl cell ghosted id
    PetscInt , allocatable :: c_local_id(:)   ! Cmpl cell local   id
    PetscInt , allocatable :: c_to_cg(:)      ! Cmpl local to cmpl global

    ! Global cmpl cell ghosted id
    PetscInt , allocatable :: cg_ghosted_id(:)

    PetscReal, allocatable :: c_ccf(:)        ! CCF by completion
    PetscReal, allocatable :: c_z  (:)        ! Cmpl elev. by compl.

    ! Completion flows for this proc
    PetscReal, allocatable :: c_flows  (    :,:) ! Comp & heat flows
    PetscReal, allocatable :: c_flowspw(    :,:) ! ..wrt Pw
    PetscReal, allocatable :: c_flowsxw(  :,:,:) ! ..wrt Xw
    PetscReal, allocatable :: c_flowsxc(:,:,:,:) ! ..wrt Xc

    ! Well flows for this proc
    PetscReal, allocatable :: w_flows  (    :)   ! Molar   well flows
    PetscReal, allocatable :: w_flowspw(    :)   ! ..wrt Pw
    PetscReal, allocatable :: w_flowsxw(  :,:)   ! ..wrt Xw
    PetscReal, allocatable :: w_flowsxc(:,:,:)   ! ..wrt Xc

    ! Well flows for all procs
    PetscReal, allocatable :: w_flowsG  (    :)  ! Molar   well flows
    PetscReal, allocatable :: w_flowsGpw(    :)  ! ..wrt Pw
    PetscReal, allocatable :: w_flowsGxw(  :,:)  ! ..wrt Xw
    PetscReal, allocatable :: w_flowsGxc(:,:,:)  ! ..wrt Xc

    ! Derivative of Xw wrt Pw
    PetscReal, allocatable :: w_tactxc(  :,:)
    PetscReal, allocatable :: w_rwxc  (  :,:)

    ! Derivative of Xw wrt Pw
    PetscReal, allocatable :: w_dxwdpw(    :)
    PetscReal, allocatable :: w_dpwdxc(  :,:)
    PetscReal, allocatable :: w_dxwdxc(:,:,:)

    ! Cmpl values this proc
    PetscReal, allocatable :: c_p   (  :)    ! Cmpl pressures
    PetscReal, allocatable :: c_t   (  :)    ! Cmpl temperature
    PetscReal, allocatable :: c_pb  (  :)    ! Cmpl bubble points
    PetscReal, allocatable :: c_mdp (:,:)    ! Cmpl molar densities
    PetscReal, allocatable :: c_sp  (:,:)    ! Cmpl saturations
    PetscReal, allocatable :: c_mob (:,:)    ! Cmpl mobilities K/vsc
    PetscReal, allocatable :: c_xo  (  :)    ! Cmpl oil mole fraction
    PetscReal, allocatable :: c_xg  (  :)    ! Cmpl gas mole fraction
    PetscReal, allocatable :: c_kgdp(:,:)    ! Cmpl mass density
    PetscReal, allocatable :: c_hp  (:,:)    ! Cmpl enthalpy density

    PetscReal, allocatable :: c_mdpX (:,:,:) ! Cmpl molar dens. wrt Xc
    PetscReal, allocatable :: c_spX  (:,:,:) ! Cmpl saturations wrt Xc
    PetscReal, allocatable :: c_mobX (:,:,:) ! Cmpl mobs~K/vsc  wrt Xc
    PetscReal, allocatable :: c_xoX  (  :,:) ! Cmpl oil mol frc wrt Xc
    PetscReal, allocatable :: c_xgX  (  :,:) ! Cmpl gas mol frc wrt Xc
    PetscReal, allocatable :: c_kgdpX(:,:,:) ! Cmpl mass dens.  wrt Xc
    PetscReal, allocatable :: c_hpX  (:,:,:) ! Cmpl enth. dens. wrt Xc

    character(len = MAXSTRINGLENGTH) :: ws_name       ! Well name

    MPI_Comm :: w_comm                                ! Well MPI communicator
    PetscBool :: w_ismp=PETSC_FALSE                   ! Is multi-proc

    PetscReal, allocatable :: xwbs    (:    )  ! Wellbore sol.
    PetscReal, allocatable :: rwbs    (:    )  ! Wellbore res.
    PetscReal, allocatable :: rwbshold(:    )  ! Stored wellbore res.
    PetscReal, allocatable :: xwbshold(:    )  ! Stored wellbore sol.
    PetscReal, allocatable :: rwbspw  (:    )  ! Wellbore res. wrt Pw
    PetscReal, allocatable :: rwbsxc  (:,:,:)  ! Wellbore res. wrt Xc
    PetscReal, allocatable :: jwbs    (:,:  )  ! Welbore Jac.
    PetscReal, allocatable :: jwbsi   (:,:  )  ! Welbore Jac. (inv)
    PetscReal, allocatable :: dxwbs   (:    )  ! Wellbore sol. change

    PetscBool :: is_solvent

    PetscInt :: ixwso
    PetscInt :: ixwsg
    PetscInt :: ixwpb
    PetscInt :: ixwss
    PetscInt :: ixwt

    PetscReal :: ws_injection_p = 1.01325d5
    PetscReal :: ws_injection_t = 15.0
    PetscReal :: ws_injection_h =  0.0

    PetscInt :: ws_unconv_t = 0
    PetscInt :: ws_unconv_w = 0
    PetscInt :: ws_unconv_b = 0
    PetscInt :: ws_MPI_errs = 0

  public ::  SolveWell, InitialiseWell, doWellMPISetup

contains

! ************************************************************************** !

subroutine initialiseWell(well_data, grid, material_auxvars, option)
  !
  ! Do well initialisation, including completion connection factors
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  use Well_Data_class
  use Grid_Grdecl_module, only : FindWellIndex
  use Grid_module
  use Material_Aux_class
  use PM_TOilIms_Aux_module

  implicit none

  class(well_data_type) , pointer :: well_data
  type (grid_type)      , pointer :: grid
  type (material_auxvar_type), pointer :: material_auxvars(:)
  type (option_type), pointer :: option

  PetscInt  :: icmpl, icmplg, ncmpl
  PetscInt  :: drilling_direction
  PetscReal :: dx, dy, dz, ccf, dx1, dx2, dh, k1, k2, &
               radius, skinfactor, thetafactor, r0

  PetscInt  :: iw
  PetscBool :: onproc, found
  PetscInt  :: local_id, ghosted_id
  character(len = MAXSTRINGLENGTH) :: ws_name       ! Well name

  ! Check surface densities have been correctly set

  call checkSurfaceDensities(option)

  ! Initialise

  iw = -1
  call well_data%getWellName(ws_name)

  found = FindWellIndex(ws_name, iw)

  if (towg_miscibility_model == TOWG_SOLVENT_TL) then
    is_solvent = PETSC_TRUE
  else
    is_solvent = PETSC_FALSE
  endif

  ! Store gravity constant

   ws_gravity = abs(option%gravity(3))

  ! Fill in completion data (from values set for wells etc.)

  call well_data%FillCmplData(iw, grid)

  nCmpl = well_data%getNCmpl()

  ! Loop over completions

  do icmpl = 1, ncmpl

    ! Get the completion location

    call well_data%GetCmplLocation(icmpl, local_id, ghosted_id, onproc, icmplg)

    ! if not on-proc then mark the completion for deletion

    if (.not.onproc) then
      call well_data%markCmplForDeletion(icmpl)
    endif
 enddo

  ! Now delete the marked completions

  call well_data%deleteMarkedCompletions()

  ! Get new number of completions and generate ccfs

  ncmpl = well_data%getNCmpl()
  do icmpl = 1, ncmpl

    ! Get completion location

    call well_data%GetCmplLocation(icmpl, local_id, ghosted_id, onproc, icmplg)

    ! Get cell dimensions

    dx = well_data%GetCmplDx(icmpl)
    dy = well_data%GetCmplDy(icmpl)
    dz = well_data%GetCmplDz(icmpl)

    ! Get cell directional properties

    drilling_direction = well_data%getCmplDrillingDirection(icmpl)
    radius             = well_data%getCmplRadius           (icmpl)
    skinfactor         = well_data%getCmplSkinFactor       (icmpl)
    thetafactor        = well_data%getCmplThetaFactor      (icmpl)

    dx1 = 0.0
    dx2 = 0.0
    dh  = 0.0
    k1  = 0.0
    k2  = 0.0

    select case(drilling_direction)
      case(X_DIRECTION)
        dx1 = dy
        dx2 = dz
        dh  = dx
        k1  = material_auxvars(ghosted_id)%permeability(perm_yy_index)
        k2  = material_auxvars(ghosted_id)%permeability(perm_zz_index)
      case(Y_DIRECTION)
        dx1 = dx
        dx2 = dz
        dh  = dy
        k1  = material_auxvars(ghosted_id)%permeability(perm_xx_index)
        k2  = material_auxvars(ghosted_id)%permeability(perm_zz_index)
      case(Z_DIRECTION)
        dx1 = dx
        dx2 = dy
        dh  = dz
        k1  = material_auxvars(ghosted_id)%permeability(perm_xx_index)
        k2  = material_auxvars(ghosted_id)%permeability(perm_yy_index)
     end select

     if ((k1 > 0.0) .and. (k2 > 0.0)) then
      r0 = (dx1**2.d0 * (k2/k1)**0.5d0 + dx2**2.d0 * (k1/k2)**0.5d0)**0.5d0 * &
           0.28d0 / ((k2/k1)**0.25d0 + (k1/k2)**0.25d0)

      ccf = 2.0d0 * PI * dh * dsqrt(k1*k2) * &
                      thetaFactor / &
                      ( dlog(r0/radius) + skinfactor )
     else
       ccf = 0.0
     endif

      call well_data%setCCF(icmpl, ccf)

  enddo

  ! Fill in well reference z value

  call well_data%setZRef(option)

  if (option%iflowmode == TOWG_MODE) then
    ws_isothermal = towg_isothermal
  elseif (option%iflowmode == TOIL_IMS_MODE) then
    ws_isothermal = toil_ims_isothermal
  else
    ws_isothermal = option%use_isothermal
  endif

  call CheckSolverAvailable(option)

end subroutine InitialiseWell

! *************************************************************************** !

subroutine doWellMPISetup(option, num_well, well_data_list)
  !
  ! Setup the connectors for the wells
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18
  !

  implicit none

  type(option_type), pointer :: option

  PetscInt, intent(in) :: num_well
  type(well_data_list_type), pointer :: well_data_list

  PetscInt, allocatable :: ncpr   (:)
  PetscInt, allocatable :: ncprall(:)
  PetscInt, allocatable :: rfortw (:)

  PetscInt :: comm_size, iwell, ncmpl, nipr, irankp
  PetscMPIInt :: myrank, ierr, nrankw, nCmplG
  PetscBool :: ismp
  MPI_Group :: group,world_group
  MPI_Comm  :: comm

  comm_size = option%mycommsize
  myrank    = option%myrank
  world_group = 0

  ! Allocate the completion count buffers

  allocate(ncpr   (comm_size))
  allocate(ncprall(comm_size))
  allocate(rfortw (comm_size))

  ! Items per rank

  nipr = 1

  ! Loop over wells

  do iwell = 1, num_well

  ! Find which procs this well exists on

  ! First set all rank completion counts to zero

    ncpr    = 0
    ncprall = 0
    rfortw  = 0

  ! Now increment the completion count for this rank

    ncmpl = getWellNCmpl(iwell, well_data_list)
    ncpr(myrank+1) = ncmpl

  ! Do a gather to find out where the other rank completions are

    ierr = 0
    call MPI_Allreduce(ncpr, ncprall, comm_size, &
                       MPI_INTEGER, MPI_SUM, option%mycomm, ierr)
    call checkErr(ierr)

  ! Is this a multi-proc well?

    ncmplg = 0
    nrankw = 0

    ismp = PETSC_FALSE
    do irankp = 1, comm_size
      if (ncprall(irankp) > 0) then
        ncmplg = ncmplg+ncprall(irankp)
        nrankw = nrankw+1
        rfortw(nrankw) = irankp-1
      endif
    enddo

  ! If more than two ranks see the well, is multi-proc

    if (nrankw >= 2) ismp = PETSC_TRUE

  ! If multi-proc, create the group and communicator

    group=0
    comm =0
    if (ismp) then
      call MPI_Comm_group (option%mycomm, world_group,ierr);
      call checkErr(ierr)
      call MPI_Group_incl (world_group  , nrankw, rfortw, group, ierr)
      call checkErr(ierr)
      call MPI_Comm_create(option%mycomm, group, comm, ierr)
      call checkErr(ierr)
    endif

    ! Set global number of completions, multi-proc flag, group and communicator

    call WellSetGlobalInfo(iwell, nrankw, ncmplg, ismp, &
                           group, comm, well_data_list)

  enddo

  call WellSetGlobalInfoSet()

  ! Free completion count buffers

  deallocate(ncpr   )
  deallocate(ncprall)
  deallocate(rfortw )

end subroutine doWellMPISetup

! *************************************************************************** !

subroutine SolveWell(aux, option, well_data, r_p)
  !
  ! Supervise well solution to obtain residual and Jacobian contributions
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Well_Data_class
  use EOS_Oil_module
  use EOS_Gas_module
  use EOS_Water_module
  use EOS_Slv_module
  use Auxiliary_module

  implicit none

  type(auxiliary_type) :: aux
  type(option_type), pointer :: option
  class(well_data_type), pointer :: well_data
  PetscReal, pointer :: r_p(:)

  PetscInt  :: local_id   = 1
  PetscInt  :: ghosted_id = 1
  PetscInt  :: icmpl, icmplg, ipoil, ipgas, ipwat, ipslv, loc
  PetscBool :: finished

  PetscBool :: welllocationfound
  PetscBool :: welltargetsfound
  PetscReal :: mw, sd, pw, jwwi, tactpw
  PetscInt  :: itertt, itt, jtt, rank, iphase, status
  PetscBool :: possible, well_solution_set, onproc, wellisdead

  ! Basic setup

  ws_unconv_t = 0
  ws_unconv_w = 0
  ws_unconv_b = 0
  ws_MPI_errs = 0

  welllocationfound = PETSC_FALSE
  welltargetsfound  = PETSC_FALSE
  possible          = PETSC_FALSE
  wellisdead        = PETSC_FALSE

  ws_nphase = option%nphase
  ws_ncomp  = option%nphase
  ws_ncompe = option%nphase+1

  if (ws_isothermal) then
    ws_nxw = ws_ncompe-2
  else
    ws_nxw = ws_ncompe-1
  endif

  w_ismp = PETSC_FALSE

  ws_ndof = option%nflowdof

  ipoil = option%oil_phase
  ipgas = option%gas_phase
  ipwat = option%liquid_phase
  ipslv = option%solvent_phase

  ws_oil = PETSC_FALSE
  ws_gas = PETSC_FALSE
  ws_wat = PETSC_FALSE
  ws_slv = PETSC_FALSE

  do iphase = 1, ws_nphase
    if (iphase == ipoil) ws_oil = PETSC_TRUE
    if (iphase == ipgas) ws_gas = PETSC_TRUE
    if (iphase == ipwat) ws_wat = PETSC_TRUE
    if (iphase == ipslv) ws_slv = PETSC_TRUE
  enddo

  ! Pointers into the wellbore solution array

  loc = 0
  if (ws_oil) then
    loc = loc+1
    ws_loc_soil = loc
  endif
  if (ws_gas) then
    loc = loc+1
    ws_loc_sgpb = loc
  endif
  if (.not.ws_isothermal) then
    loc = loc+1
    ws_loc_trel = loc
  endif
  if (ws_slv) then
    loc = loc+1
    ws_loc_sslv = loc
  endif

  ! Check for black oil cases (with dissolved gas)

  ws_is_black_oil = PETSC_FALSE
  if (     ( towg_miscibility_model == TOWG_BLACK_OIL  )  &
      .or. ( towg_miscibility_model == TOWG_SOLVENT_TL )) &
    ws_is_black_oil = PETSC_TRUE

  w_type = well_data%getType()

  call well_data%getWellInjectionPAndT(ws_injection_p, ws_injection_t)

  ws_isoilinjector = PETSC_FALSE
  ws_isgasinjector = PETSC_FALSE
  ws_iswatinjector = PETSC_FALSE
  ws_isslvinjector = PETSC_FALSE

  if (w_type == PROD_WELL_TYPE) then
    ws_isproducer  = PETSC_TRUE
    ws_isinjector  = PETSC_FALSE
    ws_injection_h = 0.0
  else
    ws_isproducer = PETSC_FALSE
    ws_isinjector = PETSC_TRUE
    if (w_type == OIL_INJ_WELL_TYPE) ws_isoilinjector = PETSC_TRUE
    if (w_type == GAS_INJ_WELL_TYPE) ws_isgasinjector = PETSC_TRUE
    if (w_type == WAT_INJ_WELL_TYPE) ws_iswatinjector = PETSC_TRUE
    if (w_type == SLV_INJ_WELL_TYPE) ws_isslvinjector = PETSC_TRUE
    call findInjectionEnthalpy()
  endif

  if (w_type == PROD_WELL_TYPE) then
    w_sign = 1.0
  else
    w_sign = -1.0
  endif

  w_zref = well_data%getZRef()

  ! Solve well if required

  w_ncmpl  = well_data%getNCmpl ()
  w_ncmplg = well_data%getNCmplG()
  w_comm   = well_data%getWellComm(w_ismp)

  if (w_ncmplg > w_ncmpl) then
    w_crossproc = PETSC_TRUE
  else
    w_crossproc = PETSC_FALSE
  endif

  status = well_data%GetWellStatus()

  if (w_ncmpl > 0 .and. status == W_STATUS_OPEN) then

    ! Debug

    call well_data%getWellName(ws_name)
    rank = option%myrank

    ! Allocate the work arrays

    call allocateWorkArrays()

    do icmpl = 1, w_ncmpl
      call well_data%GetCmplLocation(icmpl, local_id, &
                                     ghosted_id, onproc, icmplg)
      c_local_id  (icmpl) = local_id
      c_ghosted_id(icmpl) = ghosted_id
      c_to_cg     (icmpl) = icmplg
      c_ccf       (icmpl) = well_data%getCCF  (icmpl)
      c_z         (icmpl) = well_data%getCmplZ(icmpl)
    enddo

    do icmplg = 1, w_ncmplg
      call well_data%GetCmplLocationG(icmplg, ghosted_id)
      cg_ghosted_id(icmplg) = ghosted_id
    enddo

    ! Load the cell properties into arrays indexed by completion

    if (option%iflowmode == TOWG_MODE) then
      call wellSolverLoaderTOWG(aux, option)
    elseif (option%iflowmode == TOIL_IMS_MODE) then
      call wellSolverLoaderTOIL(aux, option)
    endif

    ! Setup well solution

    well_solution_set = well_data%getWellSolutionSet()
    if (well_solution_set) then
      call well_data%GetWellSolution(pw, w_pb, w_sp, w_issat, w_trel)
    else

      ! Take initial solution from first completion

      pw  = c_p(1)
      ! Set up pressure to get correct initial flow direction
      if (ws_isproducer) pw = 0.9*c_p(1) ! Positive starting drawdown
      if (ws_isinjector) pw = 1.1*c_p(1) ! Negative starting drawdown

      ! Bubble point if required

      if (ws_is_black_oil) then
        w_pb  = c_pb(1)
      endif

      ! Treat producer and injector separately

      if (ws_isproducer) then

        ! Case of producer

        do iphase = 1, ws_nphase
          w_sp(iphase) = c_sp(1, iphase)
        enddo

        if (ws_oil .and. ws_gas) then
          if (w_sp(ipoil)>1.0E-6 .and. w_sp(ipgas)<1.0E-6) then
             w_issat = PETSC_FALSE
          else
             w_issat = PETSC_TRUE
          endif
        endif

      else

       ! Case of injector

        w_sp = 0.0
        if (ws_isoilinjector .and. ws_oil) then
          w_sp(ipoil) = 1.0
          w_issat     = PETSC_FALSE
        endif
        if (ws_isgasinjector .and. ws_gas) then
          w_sp(ipgas) = 1.0
          w_issat     = PETSC_TRUE
        endif
        if (ws_iswatinjector .and. ws_wat) then
          w_sp(ipwat) = 1.0
          w_issat     = PETSC_TRUE
        endif
        if (ws_isslvinjector .and. ws_slv) then
          w_sp(ipslv) = 1.0
          w_issat     = PETSC_TRUE
        endif
      endif

      w_trel = c_t (1)

    endif

    ! Set up surface volume/mole values

    do iphase = 1, ws_nphase

      mw = 1.0
      sd = 1.0

      if (iphase == option%oil_phase    ) mw = EOSOilGetFMW()
      if (iphase == option%gas_phase    ) mw = EOSGasGetFMW()
      if (iphase == option%liquid_phase ) mw = FMWH2O
      if (iphase == option%solvent_phase) mw = EOSSlvGetFMW()

      if (iphase == option%oil_phase    ) sd = EOSOilGetSurfaceDensity()
      if (iphase == option%gas_phase    ) sd = EOSGasGetSurfaceDensity()
      if (iphase == option%liquid_phase ) sd = EOSWaterGetSurfaceDensity()
      if (iphase == option%solvent_phase) sd = EOSSlvGetSurfaceDensity()

      ws_svpm(iphase) = mw/sd
      ws_mspm(iphase) = mw
    enddo

    ! Extract well status data

    welltargetsfound = well_data%getTargets(ws_targets)
    itt = well_data%getTT()

    ! Loop over attempts to find the well target type

    do iterTT = 1, max_iterTT
      finished = solveForWellTarget(well_data,pw, option, itt, jwwi,wellisdead)
      if (finished) exit
    enddo

    ! Check for lack of convergence

    if (.not.finished) then
      ws_unconv_t = ws_unconv_t + 1
    endif

    if (wellisdead) then
      call well_data%zeroActuals()
      call well_data%zeroCmplFlows(w_ncmpl , ws_ncompe, &
                                   w_ncmplg, ws_ndof  )
    else

    ! Find full derivative of flows

      call findFullFlowDerivatives(jwwi)

    ! Update residual

      call updateMainResidual(option, r_p)

    ! Get actuals for each target for this proc
    ! (will be globalised in well_data)

      do jtt = 1, N_WELL_TT
        ws_actuals(jtt) = &
        getActualFlowForTargetType(pw, option, jtt, possible, tactpw)
      enddo

    ! Update stored well status

      call well_data%setTT(itt)
      call well_data%setNComp(ws_ncomp)
      call well_data%setCmplFlows(c_flows, c_flowsxc, &
                                  w_ncmpl, ws_ncompe, &
                                  w_ncmplg, ws_ndof, ws_isothermal)
      call well_data%setActuals(ws_actuals)

    ! Update stored well solution

      call IncrementWellWarningCount(ws_unconv_t, ws_unconv_w, &
                                     ws_unconv_b, ws_MPI_errs)
      call well_data%SetWellSolution(pw, w_pb, w_sp, w_issat, w_trel)
      call well_data%SetWellSolutionSet()

    endif

    ! Free the work arrays

    call freeWorkArrays()

  else
    call well_data%zeroActuals()
    call well_data%zeroCmplFlows(w_ncmpl , ws_ncompe, &
                                 w_ncmplg, ws_ndof  )
  endif

end subroutine SolveWell

! *************************************************************************** !

function solveForWellTarget(well_data, pw, option, itt, jwwi, wellisdead)
  !
  ! Solve for a given well target
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscBool :: solveForWellTarget
  class(well_data_type), pointer :: well_data
  PetscReal, intent(inout) :: pw
  type(option_type), pointer :: option
  PetscInt, intent(inout) :: itt
  PetscReal, intent(out)  :: jwwi
  PetscBool, intent(out) :: wellisdead

  PetscInt  :: iterpw, jtt, nconv, ierr
  PetscBool :: finished, possible, usediff
  ! fpsc: flow with positive sign convention
  PetscReal :: conv_crit, rw, jww, fpsc, ftarg, pwtarg, &
               tactpw, eps, anal, diff, rat, rwe, &
               bl(1), bg(1), ftot

  ! Initialize

  ierr      = 0
  conv_crit = 1.0E-6
  eps       = 1.0E-2
  rw        = 0.0
  jww       = 1.0
  usediff   = PETSC_FALSE
  wellisdead = PETSC_FALSE
  solveForWellTarget = PETSC_FALSE

  pwtarg = ws_targets(W_BHP_LIMIT)
  if (pwtarg<0.0) then
    if (ws_isproducer) pwtarg =        ws_atm
    if (ws_isinjector) pwtarg = 1000.0*ws_atm
  endif

  ! Set up a prediction value of the well gravity density

  call findWellboreGravityDensityPredictor()

  ! Check that the well is active by finding total rate (allowing for sign)

  call getRwAndJw(option, W_BHP_LIMIT, pwtarg, rw, jww)
  ftot = w_sign*( getActualFlowForTargetType(pwtarg, option, W_TARG_OSV, &
                  possible, tactpw) &
                 +getActualFlowForTargetType(pwtarg, option, W_TARG_GSV, &
                  possible, tactpw) &
                 +getActualFlowForTargetType(pwtarg, option, W_TARG_WSV, &
                  possible, tactpw) &
                 +getActualFlowForTargetType(pwtarg, option, W_TARG_SSV, &
                  possible, tactpw) )
  if (ftot < 0.0) wellisdead = PETSC_TRUE

  ! Check that we all agree

  wellisdead = allTrue(wellisdead)

  if (wellisdead) then
    solveForWellTarget = PETSC_TRUE
  else

  ! Do the iteration in pw for the well solution when on this target

    finished = PETSC_FALSE
    nconv    = 0
    do iterpw = 1, max_iterPW

      if (usediff) then
        call getRwAndJw(option, itt, pw+eps, rwe, jww)
        call getRwAndJw(option, itt, pw    , rw , jww)
        anal = jww
        diff = (rwe-rw)/eps
        rat = (anal+1.0E-10)/(diff+1.0E-10)
        if (abs(rat-1.0)>0.1) then
          print *, 'ws_name, itt, pw, anal, diff, rat, rw ', &
                    trim(ws_name), itt, pw, anal, diff, rat, rw
        endif
      else
        call getRwAndJw(option, itt, pw, rw, jww)
      endif

      if (w_crossproc) then
        ierr = 0
        bl(1) = rw
        call MPI_Allreduce(bl, bg, ONE_INTEGER, &
                           MPI_DOUBLE_PRECISION, MPI_MAX, w_comm, ierr)
        call checkErr(ierr)
        rw = bg(1)
      endif

      if (abs(rw) < conv_crit) nconv = nconv+1

      !  Check if finished
      finished = PETSC_FALSE
      if (nconv > 1)  finished = PETSC_TRUE

      !  Check that we all agree
      finished = allTrue(finished)

      !  Exit iteration iff all finished
      if (finished) exit

      !  Not finished - continue iterating

      if (abs(jww)>0.0) then
        jwwi = 1.0/jww
      else
        jwwi = 0.0
      endif
      pw = pw-rw*jwwi
      if (pw < ws_pwmin) pw = ws_pwmin

    enddo

  ! Check for convergence error

    if (.not.finished) then
      ws_unconv_w = ws_unconv_w + 1
    endif

  ! Check if all other targets satisfied

    solveForWellTarget = PETSC_TRUE
    do jtt = 1, N_WELL_TT
      ! No need to check current type
      if (jtt == itt) cycle
      ! Case of bhp control
      if (jtt == W_BHP_LIMIT) then
        if (      (ws_isproducer .and. ( pw < 0.999999*pwtarg )) &
             .or. (ws_isinjector .and. ( pw > 1.000001*pwtarg ))) then
          ! Switch to bhp control
          itt = jtt
          solveForWellTarget = PETSC_FALSE
        endif
      else
      ! Case of rate control - get positive convention flow
        fpsc = &
          w_sign*getActualFlowForTargetType(pw, option, jtt, possible, tactpw)

        if (possible) then
          ftarg = ws_targets(jtt)
          if (ftarg >-0.5) then
            ! If non-defaulted flow target exceeded - switch to this target
            if (fpsc > 1.000001*ftarg) then
              itt = jtt
              solveForWellTarget = PETSC_FALSE
            endif
          endif
        endif
      endif
    enddo

    ! Check that we all agree

    solveForWellTarget = allTrue(solveForWellTarget)

  endif

end function solveForWellTarget

! *************************************************************************** !

subroutine getRwAndJw(option, itt, pw, rw, jww)
  !
  ! Obtain residual and Jacobian required to solve for Pw
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  type(option_type), pointer :: option
  PetscInt, intent(in)   :: itt
  PetscReal, intent(in)  :: pw
  PetscReal, intent(out) :: rw, jww

  ! Find the well bore solution and well flows

  call findWellboreSolution(pw, option)

  ! Build up the residual for this target and solution

  jww    = 0.0
  w_rwxc = 0.0

  rw = extractResidualandJacobian(pw, option, itt, jww)

end subroutine getRwAndJw

!**************************************************************************** !

function extractResidualandJacobian(pw, option, itt, jww)
  !
  ! Get residual and Jacobian for well solver
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscReal :: extractResidualandJacobian
  PetscReal, intent(in) :: pw
  type(option_type), pointer :: option
  PetscInt, intent(in) :: itt
  PetscReal, intent(out) :: jww

  PetscBool :: possible

  PetscReal :: tact, tactpw, treq, rw

  rw       = 0.0
  w_tactxc = 0.0

  possible = PETSC_FALSE

  if (itt == W_BHP_LIMIT) then
    treq = ws_targets(itt)
    rw   = pw- treq
    jww  = 1.0
  else
    tact = getActualFlowForTargetType(pw, option, itt, possible, tactpw)
    treq = ws_targets(itt)
    rw   = w_sign*tact  -treq
    jww  = w_sign*tactpw
    extractResidualandJacobian = rw
    w_rwxc = w_sign*w_tactxc
  endif

  extractResidualandJacobian = rw

end function extractResidualandJacobian

! *************************************************************************** !

function getActualFlowForTargetType(pw, option, itt, possible, tactpw)
  !
  ! Get actual flow for a given target type
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscReal :: getActualFlowForTargetType
  PetscReal, intent(in) :: pw
  type(option_type), pointer :: option
  PetscInt, intent(in) :: itt
  PetscBool, intent(out) :: possible
  PetscReal, intent(out) :: tactpw

  PetscReal :: tact, c, co, cw
  PetscInt  :: ipoil, ipgas, ipwat, ipslv

  ! Initialise

  tact   = 0.0
  tactpw = 0.0
  possible = PETSC_FALSE

  ipoil = option%oil_phase
  ipgas = option%gas_phase
  ipwat = option%liquid_phase
  ipslv = option%solvent_phase

  ! Case of bhp limit

  if (itt == W_BHP_LIMIT) then
    tact   = pw
    tactpw = 1.0
  endif

  ! Case of various flow targets

  if ((itt == W_TARG_OSV) .and. ws_oil) then
    c = ws_svpm(ipoil)
    call incrTactX(ipoil, c, tact, tactpw)
    possible = PETSC_TRUE
  endif

  if ((itt == W_TARG_GSV) .and. ws_gas) then
    c = ws_svpm(ipgas)
    call incrTactX(ipgas, c, tact, tactpw)
    possible = PETSC_TRUE
  endif

  if ((itt == W_TARG_WSV) .and. ws_wat) then
    c = ws_svpm(ipwat)
    call incrTactX(ipwat, c, tact, tactpw)
    possible = PETSC_TRUE
  endif

  if ((itt == W_TARG_SSV) .and. ws_slv) then
    c = ws_svpm(ipslv)
    call incrTactX(ipslv, c, tact, tactpw)
    possible = PETSC_TRUE
  endif

  if ((itt == W_TARG_LSV) .and. ws_oil .and. ws_wat) then
    co = ws_svpm(ipoil)
    cw = ws_svpm(ipwat)
    call incrTactX(ipoil, co, tact, tactpw)
    call incrTactX(ipoil, cw, tact, tactpw)
    possible = PETSC_TRUE
  endif

  if ((itt == W_TARG_OM) .and. ws_oil) then
    c = ws_mspm(ipoil)
    call incrTactX(ipoil, c, tact, tactpw)
    possible = PETSC_TRUE
  endif

  if ((itt == W_TARG_GM) .and. ws_gas) then
    c = ws_mspm(ipgas)
    call incrTactX(ipgas, c, tact, tactpw)
    possible = PETSC_TRUE
  endif

  if ((itt == W_TARG_WM) .and. ws_wat) then
    c = ws_mspm(ipwat)
    call incrTactX(ipwat, c, tact, tactpw)
    possible = PETSC_TRUE
  endif

  if ((itt == W_TARG_SM).and. ws_slv) then
    c = ws_mspm(ipslv)
    call incrTactX(ipslv, c, tact, tactpw)
    possible = PETSC_TRUE
  endif

  getActualFlowForTargetType = tact

end function getActualFlowForTargetType

! *************************************************************************** !

subroutine updateMainResidual(option, r_p)
  !
  ! Insert the well flow contributiosn into the main Pflotran residual
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  type(option_type), pointer :: option
  PetscReal, pointer :: r_p(:)

  PetscReal :: res(option%nflowdof)

  PetscInt :: icmpl, icompe, local_id, istart, iend

  ! Loop over completions

  do icmpl = 1, w_ncmpl

    res = 0.0

    local_id = c_local_id(icmpl)
    do icompe = 1, ws_ncompe
      res(icompe) = c_flows(icmpl, icompe)
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
  ! Date  : 09/19/18

  implicit none

  allocate(c_local_id  (w_ncmpl          ))
  allocate(c_ghosted_id(w_ncmpl          ))
  allocate(c_to_cg     (w_ncmpl          ))

  allocate(cg_ghosted_id(w_ncmplg        ))

  allocate(c_ccf       (w_ncmpl          ))
  allocate(c_z         (w_ncmpl          ))

  allocate(ws_svpm     (        ws_ncomp ))
  allocate(ws_mspm     (        ws_ncomp ))

  allocate(c_flows     (w_ncmpl, ws_ncompe)) ! Note extra energy location
  allocate(c_flowspw   (w_ncmpl, ws_ncompe))
  allocate(c_flowsxw   (w_ncmpl, ws_ncompe, ws_nxw ))
  allocate(c_flowsxc   (w_ncmpl, ws_ncompe, w_ncmplg, ws_ndof))

  allocate(w_flows     (ws_ncompe        ))
  allocate(w_flowspw   (ws_ncompe        ))
  allocate(w_flowsxw   (ws_ncompe, ws_nxw ))
  allocate(w_flowsxc   (ws_ncompe, w_ncmplg, ws_ndof))

  allocate(w_rwxc      (w_ncmplg, ws_ndof))
  allocate(w_tactxc    (w_ncmplg, ws_ndof))

  allocate(w_flowsG    (ws_ncompe        ))
  allocate(w_flowsGpw  (ws_ncompe        ))
  allocate(w_flowsGxw  (ws_ncompe, ws_nxw ))
  allocate(w_flowsGxc  (ws_ncompe, w_ncmplg, ws_ndof))

  allocate(w_dxwdpw    (ws_nxw))
  allocate(w_dpwdxc    (        w_ncmplg, ws_ndof))
  allocate(w_dxwdxc    (ws_nxw, w_ncmplg, ws_ndof))

  allocate(c_p         (w_ncmpl          ))
  allocate(c_t         (w_ncmpl          ))
  if (ws_is_black_oil) then
    allocate(c_pb      (w_ncmpl          ))
  endif
  allocate(c_mdp       (w_ncmpl, ws_nphase))
  allocate(c_sp        (w_ncmpl, ws_nphase))
  allocate(c_mob       (w_ncmpl, ws_nphase))
  allocate(c_xo        (w_ncmpl           ))
  allocate(c_xg        (w_ncmpl           ))
  allocate(c_kgdp      (w_ncmpl, ws_nphase))
  allocate(c_hp        (w_ncmpl, ws_nphase))

  allocate(c_mdpX      (w_ncmpl, ws_nphase, ws_ndof))
  allocate(c_spX       (w_ncmpl, ws_nphase, ws_ndof))
  allocate(c_mobX      (w_ncmpl, ws_nphase, ws_ndof))
  allocate(c_xoX       (w_ncmpl           , ws_ndof))
  allocate(c_xgX       (w_ncmpl           , ws_ndof))
  allocate(c_kgdpX     (w_ncmpl, ws_nphase, ws_ndof))
  allocate(c_hpX       (w_ncmpl, ws_nphase, ws_ndof))

  ! Initialise (will stay at zero if numerical derivativess)

  c_mdpX  = 0.0
  c_spX   = 0.0
  c_mobX  = 0.0
  c_xoX   = 0.0
  c_xgX   = 0.0
  c_kgdpX = 0.0
  c_hpX   = 0.0

  allocate(w_sp        (ws_nphase        ))
  allocate(w_spxw      (ws_nphase, ws_nxw))

  allocate(w_mdp       (ws_nphase        ))
  allocate(w_mdppw     (ws_nphase        ))
  allocate(w_mdpxw     (ws_nphase, ws_nxw))

  allocate(w_mdc       (ws_ncomp         ))
  allocate(w_mdcpw     (ws_ncomp         ))
  allocate(w_mdcxw     (ws_ncomp , ws_nxw))

  allocate(w_zmf       (ws_ncomp         ))
  allocate(w_zmfpw     (ws_ncomp         ))
  allocate(w_zmfxw     (ws_ncomp , ws_nxw))

  allocate(w_kgdp      (ws_nphase        ))
  allocate(w_kgdppw    (ws_nphase        ))
  allocate(w_kgdpxw    (ws_nphase, ws_nxw))

  allocate(w_hp        (ws_nphase        ))
  allocate(w_hpdpw     (ws_nphase        ))
  allocate(w_hpdxw     (ws_nphase, ws_nxw))

  allocate(w_hdxw      (           ws_nxw))
  allocate(w_hpmdxw    (           ws_nxw))

  allocate(w_gdxw      (           ws_nxw))

  ! Wellbore solution values

  allocate(xwbs        (        ws_nxw))
  allocate(rwbs        (        ws_nxw))
  allocate(rwbshold    (        ws_nxw))
  allocate(xwbshold    (        ws_nxw))
  allocate(rwbspw      (        ws_nxw))
  allocate(rwbsxc      (ws_nxw, w_ncmplg, ws_ndof))
  allocate(jwbs        (ws_nxw, ws_nxw))
  allocate(jwbsi       (ws_nxw, ws_nxw))
  allocate(dxwbs       (        ws_nxw))

end subroutine allocateWorkArrays

!*****************************************************************************!

subroutine freeWorkArrays()
  !
  ! De-allocate well solve work arrays
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  deallocate(c_local_id  )
  deallocate(c_ghosted_id)
  deallocate(c_to_cg     )

  deallocate(cg_ghosted_id)

  deallocate(c_ccf       )
  deallocate(c_z         )

  deallocate(ws_svpm     )
  deallocate(ws_mspm     )

  deallocate(c_flows     )
  deallocate(c_flowspw   )
  deallocate(c_flowsxw   )
  deallocate(c_flowsxc   )

  deallocate(w_flows     )
  deallocate(w_flowspw   )
  deallocate(w_flowsxw   )
  deallocate(w_flowsxc   )

  deallocate(w_rwxc      )
  deallocate(w_tactxc    )

  deallocate(w_flowsG    )
  deallocate(w_flowsGpw  )
  deallocate(w_flowsGxw  )
  deallocate(w_flowsGxc  )

  deallocate(w_dxwdpw    )
  deallocate(w_dxwdxc    )
  deallocate(w_dpwdxc    )

  deallocate(c_p         )
  deallocate(c_t         )
  if (ws_is_black_oil) then
    deallocate(c_pb        )
  endif
  deallocate(c_mdp       )
  deallocate(c_sp        )
  deallocate(c_mob       )
  deallocate(c_xo        )
  deallocate(c_xg        )
  deallocate(c_kgdp      )
  deallocate(c_hp        )

  deallocate(c_mdpX      )
  deallocate(c_spX       )
  deallocate(c_mobX      )
  deallocate(c_xoX       )
  deallocate(c_xgX       )
  deallocate(c_kgdpX     )
  deallocate(c_hpX       )

  deallocate(w_sp        )
  deallocate(w_spxw      )

  deallocate(w_mdp       )
  deallocate(w_mdppw     )
  deallocate(w_mdpxw     )

  deallocate(w_mdc       )
  deallocate(w_mdcpw     )
  deallocate(w_mdcxw     )

  deallocate(w_zmf       )
  deallocate(w_zmfpw     )
  deallocate(w_zmfxw     )

  deallocate(w_kgdp      )
  deallocate(w_kgdppw    )
  deallocate(w_kgdpxw    )

  deallocate(w_hp        )
  deallocate(w_hpdpw     )
  deallocate(w_hpdxw     )

  deallocate(w_hdxw      )
  deallocate(w_hpmdxw    )

  deallocate(w_gdxw      )

  deallocate(xwbs        )
  deallocate(rwbs        )
  deallocate(rwbshold    )
  deallocate(xwbshold    )
  deallocate(rwbspw      )
  deallocate(rwbsxc      )
  deallocate(jwbs        )
  deallocate(jwbsi       )
  deallocate(dxwbs       )

end subroutine freeWorkArrays

! *************************************************************************** !

subroutine findAllCompletionFlows(pw, option)
  !
  ! Find flows for all completions
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscReal, intent(in) :: pw
  type(option_type), pointer :: option
  PetscInt :: icmpl, icmplg

  c_flows   = 0.0
  c_flowspw = 0.0
  c_flowsxw = 0.0
  c_flowsxc = 0.0
  do icmpl = 1, w_ncmpl
   icmplg = c_to_cg(icmpl)
   call findCompletionFlows(pw, option, icmpl, icmplg)
  enddo

end subroutine findAllCompletionFlows

! *************************************************************************** !

subroutine findCompletionFlows(pw, option, icmpl, icmplg)
  !
  ! Find completion flows for one completion
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscReal, intent(in) :: pw
  type(option_type), pointer :: option
  PetscInt , intent(in) :: icmpl
  PetscInt , intent(in) :: icmplg

  PetscInt  :: jdof

  PetscBool :: componentInPhase, is_oil_in_oil, is_gas_in_oil
  PetscInt  :: iphase, jphase, icomp, eid, ixw
  PetscReal :: xmf, ccf, z, dhh, pdd, pddpw, pddpc, mob, mobt, mdp, &
               flow, flowpw, flowxw, flowpc, flowxc, &
               hp, hpP, mdc, mdcpw
  PetscReal :: mobX, mdpX, xmfX, hpX
  PetscReal, allocatable :: pddxw(:)

  allocate(pddxw(ws_nxw))
  pddxw = 0.0

  ! Initialise

  eid = option%energy_id

  ! Extract completion connection factor and drawdown

  ccf = c_ccf(icmpl)
  z   = c_z  (icmpl)
  dhh = z-w_zref

  pdd   = c_P  (icmpl)-pw +w_gd  *dhh
  pddpw =             -1.0+w_gdpw*dhh
  pddpc = 1.0

  do ixw = 1, ws_nxw
    pddxw(ixw) = w_gdxw(ixw)*dhh
  enddo

  mobt = 0.0
  do iphase = 1, ws_nphase
    mobt = mobt+c_mob(icmpl, iphase)
  enddo

  ! Loop over phases and components building up component flows

  if (ws_is_black_oil) then
    do iphase = 1, ws_nphase
      do icomp = 1, ws_ncomp

        componentInPhase = PETSC_FALSE
        is_oil_in_oil    = PETSC_FALSE
        is_gas_in_oil    = PETSC_FALSE

        ! Oil phase in black oil mode is special case

        if (iphase == option%oil_phase) then
          if (icomp == option%oil_phase) is_oil_in_oil = PETSC_TRUE
          if (icomp == option%gas_phase) is_gas_in_oil = PETSC_TRUE
        endif

        ! OK if phase and component match or dissolved gas in producer

        if ((iphase == icomp) .or. &
             (is_gas_in_oil.and.ws_isproducer)) componentInPhase = PETSC_TRUE

        if (componentInPhase) then

          xmf = 1.0d0
          if (Pdd > 0.0) then
            if (is_oil_in_oil) xmf = c_xo(icmpl)
            if (is_gas_in_oil) xmf = c_xg(icmpl)
          endif

          mob  = 0.0
          mdp  = 0.0
          hp   = 0.0
          hpP  = 0.0

          ! Increment molar flow by component

          if (pdd > 0.0) then

            ! Producing completion: mobility and molar density
            ! are cell functions

            mob   = c_mob(icmpl, iphase)
            mdp   = c_mdp(icmpl, iphase)

            flow   = ccf*xmf*mob*mdp*pdd
            flowpw = ccf*xmf*mob*mdp*pddpw

            c_flows  (icmpl, icomp) = c_flows  (icmpl, icomp)+flow
            c_flowspw(icmpl, icomp) = c_flowspw(icmpl, icomp)+flowpw

            ! For producing case:
            ! Mobility and molar density are Xc functions
            ! Drawdowns                  are Xw functions

            ! Wellbore solution derivatives
            do ixw = 1, ws_nxw
              flowxw = ccf*xmf*mob*mdp*pddxw(ixw)
              c_flowsxw(icmpl, icomp, ixw) = &
              c_flowsxw(icmpl, icomp, ixw)+flowxw
            enddo

            ! Cell solution derivatives (diagonal in completion index)

            ! Wrt cell pressure

            flowpc = ccf*xmf*mob*mdp*pddpc
            c_flowsxc(icmpl, icomp, icmplg, ws_xpr) &
          = c_flowsxc(icmpl, icomp, icmplg, ws_xpr) + flowpc

            ! Wrt cell properties
            do jdof = 1, ws_ndof

              ! Wrt molar density and mobility
              mobX   = c_mobX(icmpl, iphase, jdof)
              mdpX   = c_mdpX(icmpl, iphase, jdof)
              flowxc = ccf*xmf*(mobX*mdp+mob*mdpX)*pdd

              c_flowsxc(icmpl, icomp, icmplg, jdof) &
            = c_flowsxc(icmpl, icomp, icmplg, jdof) + flowxc

              ! Wrt mole fraction
              xmfX = 0.0
              if (is_oil_in_oil) xmfX = c_xoX(icmpl, jdof)
              if (is_gas_in_oil) xmfX = c_xgX(icmpl, jdof)
              flowxc = ccf*xmfX*mob*mdp*pdd

              c_flowsxc(icmpl, icomp, icmplg, jdof) &
            = c_flowsxc(icmpl, icomp, icmplg, jdof) + flowxc

            enddo

          else

            ! For injecting case:
            ! Mobility and molar density are Xc functions
            ! Drawdowns                  are Xw functions

            mdc   = w_mdc  (icomp)
            mdcpw = w_mdcpw(icomp)

            flow   = ccf*mobt*  mdc  *pdd
            flowpw = ccf*mobt*( mdcpw*pdd   &
                               +mdc  *pddpw )

            c_flows  (icmpl, icomp) = c_flows  (icmpl, icomp)+flow
            c_flowspw(icmpl, icomp) = c_flowspw(icmpl, icomp)+flowpw

            ! Wellbore solution derivatives
            do ixw = 1, ws_nxw
              flowxw = ccf*mobt*( w_mdcxw(icomp, ixw)*pdd        &
                                 +mdc                *pddxw(ixw) )
              c_flowsxw(icmpl, icomp, ixw) = &
              c_flowsxw(icmpl, icomp, ixw)+flowxw
            enddo

            ! Cell solution derivatives (diagonal in completion index)

            ! Wrt cell pressure via drawdown

            flowpc = ccf*mobt*mdc*pddpc
            c_flowsxc(icmpl, icomp, icmplg, ws_xpr) &
          = c_flowsxc(icmpl, icomp, icmplg, ws_xpr) + flowpc

            ! Wrt cell voidage mobility via mobt

            do jphase = 1, ws_nphase
              do jdof = 1, ws_ndof
                  flowxc = ccf*c_mobX(icmpl, jphase, jdof)*mdc*pdd
                  c_flowsxc(icmpl, icomp, icmplg, jdof) &
                = c_flowsxc(icmpl, icomp, icmplg, jdof)+flowxc
              enddo
            enddo

          endif

        endif

      enddo
    enddo

  else

    ! Non-black-oil case

    do iphase = 1, ws_nphase

      icomp = iphase

      ! Increment molar flow by component

      if (pdd > 0.0) then

        ! Producing completion: mobility and molar density are cell functions

        mob = c_mob(icmpl, iphase)
        mdp = c_mdp(icmpl, iphase)

        flow   = ccf*mob*mdp*pdd
        flowpw = ccf*mob*mdp*pddpw

        c_flows  (icmpl, icomp) = c_flows  (icmpl, icomp)+flow
        c_flowspw(icmpl, icomp) = c_flowspw(icmpl, icomp)+flowpw

        ! For producing case, mobility and molar density are cell functions
        ! Only the drawdown is a pw function

        do ixw = 1, ws_nxw
          flowxw = ccf*mob*mdp*pddxw(ixw)
          c_flowsxw(icmpl, icomp, ixw) = c_flowsxw(icmpl, icomp, ixw)+flowxw
        enddo

        ! Cell solution derivatives (diagonal in completion index)

        ! Wrt cell pressure
        flowpc = ccf*mob*mdp*pddpc
        c_flowsxc(icmpl, icomp, icmplg, ws_xpr) &
      = c_flowsxc(icmpl, icomp, icmplg, ws_xpr) + flowpc

        ! Wrt cell properties
        do jdof = 1, ws_ndof

          ! Wrt molar density and mobility (this completion)
          mobX   = c_mobX(icmpl, iphase, jdof)
          mdpX   = c_mdpX(icmpl, iphase, jdof)
          flowxc = ccf*(mobX*mdp+mob*mdpX)*pdd

          c_flowsxc(icmpl, icomp, icmplg, jdof) &
        = c_flowsxc(icmpl, icomp, icmplg, jdof) + flowxc

        enddo

      else

        ! Injecting completion: molar densities are well functions

        mdc   = w_mdc  (icomp)
        mdcpw = w_mdcpw(icomp)

        flow   = ccf*mobt*  mdc  *pdd
        flowpw = ccf*mobt*( mdcpw*pdd   &
                           +mdc  *pddpw )

        c_flows  (icmpl, icomp) = c_flows  (icmpl, icomp) + flow
        c_flowspw(icmpl, icomp) = c_flowspw(icmpl, icomp) + flowpw

        do ixw = 1, ws_nxw
          flowxw = ccf*mobt*( w_mdcxw(icomp, ixw)*pdd        &
                             +mdc                *pddxw(ixw) )
          c_flowsxw(icmpl, icomp, ixw) = c_flowsxw(icmpl, icomp, ixw) + flowxw
        enddo

        ! Cell solution derivatives (diagonal in completion index)

         flowpc = ccf*mobt*mdc*pddpc
         c_flowsxc(icmpl, icomp, icmplg, ws_xpr) &
       = c_flowsxc(icmpl, icomp, icmplg, ws_xpr) + flowpc

         ! Sum over elements of mobt

         do jphase = 1, ws_nphase
          do jdof = 1, ws_ndof
            flowxc = ccf*c_mobX(icmpl, jphase, jdof)*mdc*pdd
            c_flowsxc(icmpl, icomp, icmplg, jdof) &
          = c_flowsxc(icmpl, icomp, icmplg, jdof) + flowxc
          enddo
        enddo

      endif

    enddo

  endif

  ! Thermal part

  if (pdd > 0.0) then

    ! Thermal producing completion

    do iphase = 1, ws_nphase

      mob = c_mob(icmpl, iphase)
      mdp = c_mdp(icmpl, iphase)
      hp  = c_hp (icmpl, iphase)

      flow   = ccf*mob*mdp*pdd  *hp
      flowpw = ccf*mob*mdp*pddpw*hp

      c_flows  (icmpl, eid) = c_flows  (icmpl, eid) + flow
      c_flowspw(icmpl, eid) = c_flowspw(icmpl, eid) + flowpw

      do ixw = 1, ws_nxw
        flowxw = ccf*mob*mdp*pddxw(ixw)*hp
        c_flowsxw(icmpl, eid, ixw) = c_flowsxw(icmpl, eid, ixw) + flowxw
      enddo

      ! Cell solution derivatives (diagonal in completion index)

      ! Wrt cell pressure
      flowpc = ccf*mob*mdp*pddpc*hp
      c_flowsxc(icmpl, eid, icmplg, ws_xpr) &
    = c_flowsxc(icmpl, eid, icmplg, ws_xpr) + flowpc

        ! Wrt cell properties
        do jdof = 1, ws_ndof

          ! Wrt molar density, mobility and specific enthalpy
          mobX   = c_mobX(icmpl, iphase, jdof)
          mdpX   = c_mdpX(icmpl, iphase, jdof)
          hpX    = c_hpX (icmpl, iphase, jdof)
          flowxc = ccf*(mobX*mdp*hp+mob*mdpX*hp+mob*mdp*hpX)*pdd

          c_flowsxc(icmpl, eid, icmplg, jdof) &
        = c_flowsxc(icmpl, eid, icmplg, jdof) + flowxc

        enddo

    enddo
  else

    ! Thermal injecting completion

    flow   = ccf*mobt* pdd  *w_h
    flowpw = ccf*mobt*(pddpw*w_h+pdd*w_hdpw)

    c_flows  (icmpl, eid) = c_flows  (icmpl, eid) + flow
    c_flowspw(icmpl, eid) = c_flowspw(icmpl, eid) + flowpw

    do ixw = 1, ws_nxw
      flowxw = ccf*mobt*(pddxw(ixw)*w_h+pdd*w_hdxw(ixw))
      c_flowsxw(icmpl, eid, ixw) = c_flowsxw(icmpl, eid, ixw) + flowxw
    enddo

    ! Cell solution derivatives (diagonal in completion index)

    flowpc = ccf*mobt*pddpc*w_h
    c_flowsxc(icmpl, eid, icmplg, ws_xpr) &
  = c_flowsxc(icmpl, eid, icmplg, ws_xpr) + flowpc

    ! Sum over elements of mobt

    do jphase = 1, ws_nphase
      do jdof = 1, ws_ndof
        flowxc = ccf*c_mobX(icmpl, jphase, jdof)*pdd*w_h
        c_flowsxc(icmpl, eid, icmplg, jdof) &
      = c_flowsxc(icmpl, eid, icmplg, jdof) + flowxc
      enddo
    enddo

  endif

  deallocate(pddxw)

end subroutine findCompletionFlows

! *************************************************************************** !

subroutine buildWellFlows1()
  !
  ! Build the well flows from the completion flows
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscInt :: icompe, icmpl, ierr, ixw, jcmplg, jdof

  ierr = 0

  w_flows     = 0.0
  w_flowspw   = 0.0
  w_flowsxw   = 0.0
  w_flowsxc   = 0.0

  w_flowsG    = 0.0
  w_flowsGpw  = 0.0
  w_flowsGxw  = 0.0
  w_flowsGxc  = 0.0

  do icompe = 1, ws_ncompe
    do icmpl = 1, w_ncmpl
      w_flows  (icompe) = w_flows  (icompe) + c_flows  (icmpl, icompe)
      w_flowspw(icompe) = w_flowspw(icompe) + c_flowspw(icmpl, icompe)
      do ixw = 1, ws_nxw
        w_flowsxw(icompe, ixw) &
      = w_flowsxw(icompe, ixw) + c_flowsxw(icmpl, icompe, ixw)
      enddo
      do jcmplg = 1, w_ncmplg
        do jdof = 1, ws_ndof
           w_flowsxc(icompe, jcmplg, jdof) &
        =  w_flowsxc(icompe, jcmplg, jdof) &
         + c_flowsxc(icmpl, icompe, jcmplg, jdof)
        enddo
      enddo
    enddo
  enddo

  ! Build well flows across procs if required

  if (w_crossproc) then
    call MPI_Allreduce( w_flows   , w_flowsG   , ws_ncompe,                  &
                        MPI_DOUBLE_PRECISION, MPI_SUM, w_comm, ierr          )
    call checkErr(ierr)
    call MPI_Allreduce( w_flowspw , w_flowsGpw , ws_ncompe,                  &
                        MPI_DOUBLE_PRECISION, MPI_SUM, w_comm, ierr          )
    call checkErr(ierr)
    call MPI_Allreduce( w_flowsxw , w_flowsGxw , ws_ncompe*ws_nxw,           &
                        MPI_DOUBLE_PRECISION, MPI_SUM, w_comm, ierr          )
    call checkErr(ierr)
    call MPI_Allreduce( w_flowsxc , w_flowsGxc , ws_ncompe*w_ncmplg*ws_ndof, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, w_comm, ierr          )
    call checkErr(ierr)
  else
    w_flowsG    = w_flows
    w_flowsGpw  = w_flowspw
    w_flowsGxw  = w_flowsxw
    w_flowsGxc  = w_flowsxc
  endif

end subroutine buildWellFlows1

! *************************************************************************** !

subroutine buildWellFlows2()
  !
  ! Build the well flows from the completion flows
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscInt :: icompe, icmpl, ierr, jcmplg, jdof

  ierr = 0

  w_flows     = 0.0
  w_flowspw   = 0.0
  w_flowsxc   = 0.0

  w_flowsG    = 0.0
  w_flowsGpw  = 0.0
  w_flowsGxc  = 0.0

  do icompe = 1, ws_ncompe
    do icmpl = 1, w_ncmpl
      w_flows  (icompe) = w_flows  (icompe) + c_flows  (icmpl, icompe)
      w_flowspw(icompe) = w_flowspw(icompe) + c_flowspw(icmpl, icompe)
      do jcmplg = 1, w_ncmplg
        do jdof = 1, ws_ndof
           w_flowsxc(icompe, jcmplg, jdof) &
        =  w_flowsxc(icompe, jcmplg, jdof) &
         + c_flowsxc(icmpl, icompe, jcmplg, jdof)
        enddo
      enddo
    enddo
  enddo

  ! Build well flows across procs if required

  if (w_crossproc) then
    call MPI_Allreduce( w_flows   , w_flowsG   , ws_ncompe,                  &
                        MPI_DOUBLE_PRECISION, MPI_SUM, w_comm, ierr          )
    call checkErr(ierr)
    call MPI_Allreduce( w_flowspw , w_flowsGpw , ws_ncompe,                  &
                        MPI_DOUBLE_PRECISION, MPI_SUM, w_comm, ierr          )
    call checkErr(ierr)
    call MPI_Allreduce( w_flowsxc , w_flowsGxc , ws_ncompe*w_ncmplg*ws_ndof, &
                        MPI_DOUBLE_PRECISION, MPI_SUM, w_comm, ierr          )
    call checkErr(ierr)
  else
    w_flowsG    = w_flows
    w_flowsGpw  = w_flowspw
    w_flowsGxc  = w_flowsxc
  endif

end subroutine buildWellFlows2

! *************************************************************************** !

subroutine findWellboreSolution(pw, option)
  !
  ! Find gravity*(wellbore mass density) and well bore composition
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscReal, intent(in) :: pw
  type(option_type), pointer :: option

  PetscInt  :: iterwbs, ir, irw, ixw, icmpl, icompe, jdof, jcmplg, ierr, &
               itry, mtry
  PetscReal :: epsconv,epsdiff,epssoln, Rnorm, so, sg, pb, deti, &
               diff, anal, rat, pwe, sum
  PetscReal :: lam, rNormhold
  PetscBool :: finished, usediff, usediff2, issathold

  ! Initialise

  ierr = 0
  sum = 0.0
  epsconv = 1.0E-6
  epsdiff = 1.0E-6
  epssoln = 1.0E-6
  finished = PETSC_FALSE
  usediff  = PETSC_FALSE
  usediff2 = PETSC_FALSE

  ! Pack wellbore solution into xwbs

  call packWellboreSolution(option)

  ! Iterate to find wellbore solution

  do iterwbs = 1, max_iterWBS

   ! Get residuals for difference checks

    if (usediff .or. usediff2) then
      call getRwbsAndJwbs(pw, option,sum)
      do ir = 1, ws_nxw
        rwbshold(ir) = rwbs(ir)
        xwbshold(ir) = xwbs(ir)
      enddo
    endif

    if (usediff) then
      do ixw = 1, ws_nxw

       xwbs(ixw) = xwbshold(ixw)

        epsdiff = 1.0E-6
        if (ixw == 1 .and. (xwbs(ixw)>0.5) ) epsdiff = -1.0E-6
        if (ixw == 2 ) epsdiff = 1.0E-2

        xwbs(ixw) = xwbs(ixw)+epsdiff
        call getRwbsAndJwbs(pw, option,sum)
        do ir = 1, ws_nxw
          diff = (rwbs(ir)-rwbshold(ir))/epsdiff
          anal = jwbs(ir, ixw)
          rat  = (diff+1.0E-20)/(anal+1.0E-20)
          print *, 'RwbsXw:ir, ix, d, a, rat ' , ir, ixw, diff, anal, rat
        enddo

        xwbs(ixw) = xwbshold(ixw)
      enddo
    endif
    if (usediff2) then
      epsdiff = 1.0
      pwe = pw+epsdiff
      call getRwbsAndJwbs(pwe, option,sum)
      do ir = 1, ws_nxw
        diff = (rwbs(ir)-rwbshold(ir))/epsdiff
        anal =  rwbspw(ir)
        rat  = (diff+1.0E-20)/(anal+1.0E-20)
        print *, 'dRwbsPw:ir, d, a, r ', ir, diff, anal, rat
      enddo
    endif

    ! Get residual and norm

    call getRwbsAndJwbs(pw, option,sum)
    Rnorm = getRnormwbs()

    !  Does this norm indicate we are finished?

    finished = PETSC_FALSE
    if (Rnorm < epsconv .and. iterwbs>1 ) finished = PETSC_TRUE

    !  Check that we all agree

    finished = allTrue(finished)

    ! Exit if finished

    if (finished) then
      call invertJacobian(deti)
      exit
    endif

    ! Invert Jacobian

    call invertJacobian(deti)

    ! Find solution change

    dxwbs = 0.0

    do ixw = 1, ws_nxw
      do ir = 1, ws_nxw
        dxwbs(ixw) = dxwbs(ixw) - jwbsi(ixw, ir)*rwbs(ir)
      enddo
    enddo

    ! Do update

    xwbshold  = xwbs
    issathold = w_issat
    Rnormhold = Rnorm

    ! Loop over attempts to update using a backtracker

    lam=1.0
    mtry = 4
    do itry=1,mtry

      do ixw = 1, ws_nxw
        xwbs(ixw) = xwbs(ixw)+lam*dxwbs(ixw)
      enddo

      if (xwbs(1) < 0.0) xwbs(1) = 0.0
      if (xwbs(1) > 1.0) xwbs(1) = 1.0

      ! Check for state change

      if (ws_is_black_oil) then
        so = xwbs(1)
        if (w_issat) then
          sg = xwbs(2)
          pb = pw
          if (sg < 0.0 .and. so > epssoln) then
            w_issat = PETSC_FALSE
            xwbs(2) = pb - epssoln
          endif
          if (sg>1.0) then
            sg = 1.0
            xwbs(2) = sg
          endif
        else
          sg = 0.0
          pb = xwbs(2)
          if (pb > pw .or. so <= epssoln) then
            w_issat = PETSC_TRUE
            xwbs(2) = epssoln
          endif
        endif
      endif

      ! Limit to positive now change detected

      do ixw = 1, ws_nxw
        if (xwbs(ixw) < 0.0) xwbs(ixw) = 0.0
      enddo

      ! Get residual and norm

      call getRwbsAndJwbs(pw, option,sum)
      Rnorm = getRnormwbs()

      ! Check if gain, leave if so

      if ((Rnorm < epsconv ) .or. (Rnorm<Rnormhold) .or. (itry==mtry)) exit

      lam=0.5*lam

      xwbs    =  xwbshold
      w_issat = issathold

    enddo
  enddo

  ! Check for lack of convergence

  if (.not.finished) then
    ws_unconv_b = ws_unconv_b + 1
  endif

  ! Use implicit function theorem to find derivatives of Xw wrt Pw and Xc

  ! Find dXw/dPw = -(1/(dRwbs/dXw))*(dRbs/dPw)

  w_dxwdpw = 0.0
  do ixw = 1, ws_nxw

    do irw = 1, ws_nxw
      w_dxwdpw(ixw) = w_dxwdpw(ixw) - jwbsi(ixw, irw)*rwbspw(irw)
    enddo

  enddo

  ! Sort out dXw/dXc = -(1/(dRwbs/dXw))*(dRwbs/dXc)

  w_dxwdxc = 0.0

  do ixw = 1, ws_nxw
    do jcmplg = 1, w_ncmplg
      do jdof = 1, ws_ndof

        do irw = 1, ws_nxw
           w_dxwdxc(ixw, jcmplg, jdof) &
        =  w_dxwdxc(ixw, jcmplg, jdof) &
         - jwbsi(ixw, irw)*rwbsxc(irw, jcmplg, jdof)
        enddo

      enddo
    enddo
  enddo

  ! Prepare the output flows by eliminating wellbore solution

  do icompe = 1, ws_ncompe
    do icmpl = 1, w_ncmpl

      sum = 0.0
      do ixw = 1, ws_nxw
        sum = sum + c_flowsxw(icmpl, icompe, ixw)*w_dxwdpw(ixw)
      enddo
         c_flowspw(icmpl, icompe) &
       = c_flowspw(icmpl, icompe) + sum

      do jcmplg = 1, w_ncmplg
        do jdof = 1, ws_ndof

          sum = 0.0
          do ixw = 1, ws_nxw
            sum = sum + c_flowsxw(icmpl, icompe, ixw) &
                       *w_dxwdxc(ixw, jcmplg, jdof)
          enddo

            c_flowsxc(icmpl, icompe, jcmplg, jdof) &
          = c_flowsxc(icmpl, icompe, jcmplg, jdof) + sum
        enddo
      enddo
    enddo
  enddo

  ! Assemble the total derivatives to find well flows

  call buildWellFlows2() ! Includes other procs sums

end subroutine findWellboreSolution

!*****************************************************************************!

subroutine packWellboreSolution(option)
  !
  ! Pack the wellbore solution into the well solution array
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  type(option_type), pointer :: option

  if (ws_oil) then
    xwbs(ws_loc_soil) = w_sp(option%oil_phase)
  endif

  if (ws_gas) then
    if (w_issat) then
      xwbs(ws_loc_sgpb) = w_sp(option%gas_phase)
    else
      xwbs(ws_loc_sgpb) = w_pb
    endif
  endif

  if (.not.ws_isothermal) then
    xwbs(ws_loc_trel) = w_trel
  endif

  if (is_solvent) then
    xwbs(ws_loc_sslv) = w_sp(option%solvent_phase)
  endif

end subroutine packWellboreSolution

! *************************************************************************** !

subroutine unpackWellboreSolution(pw, so, sg, sw, ss, pb, t)
  !
  ! Unpack the wellbore solution from the well solution array
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscReal, intent(in)  :: pw
  PetscReal, intent(out) :: so, sg, sw, ss, pb, t

  PetscReal :: sn, sni

  ! Initialise

  ixwsg = 0
  ixwpb = 0
  ixwss = 0
  ixwt  = 0

  so =  0.0
  sg =  0.0
  ss =  0.0
  sw =  0.0
  t  = 15.0

  ! Set up saturations and pointers

  if (ws_oil) then
    ixwso = ws_loc_soil
    ! Oil in (0 to 1)
    if (xwbs(ixwso)<0.0) xwbs(ixwso) = 0.0
    if (xwbs(ixwso)>1.0) xwbs(ixwso) = 1.0
    so = xwbs(ixwso)
  endif

  if (ws_gas) then
    if (w_issat .or. (.not.ws_is_black_oil)) then
      ixwsg = ws_loc_sgpb
      ! Gas  >0.0 if not black oil and <1.0 all cases
      if ((.not.ws_is_black_oil) .and. (xwbs(ixwsg)<0.0)) xwbs(ixwsg) = 0.0
      if (                             (xwbs(ixwsg)>1.0)) xwbs(ixwsg) = 1.0
      sg    = xwbs(ixwsg)
      pb    = pw
    else
      ixwpb = ws_loc_sgpb
      sg    = 0.0
      ! Pb >0.0 all cases
      if (xwbs(ixwpb)<0.0) xwbs(ixwpb) = 0.0
      pb    = xwbs(ixwpb)
    endif
  endif

  if (.not.ws_isothermal) then
    t    = xwbs(ws_loc_trel)
    ixwt = ws_loc_trel
  endif

  if (is_solvent) then
    ixwss = ws_loc_sslv
    if (xwbs(ixwss)<0.0) xwbs(ixwss) = 0.0
    if (xwbs(ixwss)>1.0) xwbs(ixwss) = 1.0
    ss = xwbs(ixwss)
  endif

  ! Water is always present and dependent

  sn = so+sg+ss
  sw = 1.0-sn

  if (sn>1.0) then
    sni = 1.0/sn
    if (ws_oil) then
      xwbs(ws_loc_soil) = xwbs(ws_loc_soil)*sni
      so                = xwbs(ws_loc_soil)
    endif
    if ( ws_gas .and. (w_issat .or. (.not.ws_is_black_oil))) then
      xwbs(ws_loc_sgpb) = xwbs(ws_loc_sgpb)*sni
      sg                = xwbs(ws_loc_sgpb)
    endif
    if (ws_slv) then
      xwbs(ws_loc_sslv) = xwbs(ws_loc_sslv)*sni
      ss                = xwbs(ws_loc_sslv)
    endif
    sn = so+sg+ss
    sw = 1.0-sn
  endif

end subroutine unpackWellboreSolution

! *************************************************************************** !

subroutine getRwbsAndJwbs(pw, option, sum)
  !
  ! Get the residual and Jacobian for the wellbore solution
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscReal, intent(in) :: pw
  type(option_type), pointer :: option
  PetscReal, intent(out) :: sum

  PetscReal :: so, sg, sw, ss, pb, t
  PetscReal :: sumpw
  PetscInt  :: icomp, icompc, icomph, nic, &
               ixw, ipoil, ipgas, ipwat, ipslv, eid
  PetscInt  :: jcmplg, jdof, ippref

  PetscReal, allocatable :: sumxw(:)
  PetscReal, allocatable :: sumxc(:,:)

  ! Initialise

  allocate(sumxw(ws_nxw ))
  allocate(sumxc(w_ncmplg, ws_ndof))

  nic = ws_nxw
  if (.not.ws_isothermal) nic = ws_nxw-1

  icomph = -1

  ipoil = option%oil_phase
  ipgas = option%gas_phase
  ipwat = option%liquid_phase
  ipslv = option%solvent_phase

  ippref = 1
  if (ws_isproducer) then
    if (ws_oil) then
       ippref = ipoil
     else
       ippref = ipwat
    endif
   endif

  eid = option%energy_id

  ! Zero wellbore residual, derivative wrt Pw and Jacobian

  rwbs   = 0.0
  rwbspw = 0.0
  jwbs   = 0.0
  rwbsxc = 0.0

  ! Unpack the wellbore solution from xwbs

  call unpackWellboreSolution(pw, so, sg, sw, ss, pb, t)

  ! Get wellbore properties

  call findWellborePropertiesAndDerivatives(pw, so, sg, sw, ss, pb, t, option)

  ! Now calculate all the completion flows

  call findAllCompletionFlows(pw, option)

  ! Sum to find well flows

  call buildWellFlows1() ! Includes other procs sums

  ! Find if producing or not

  sum   = 0.0
  sumpw = 0.0
  sumxw = 0.0
  sumxc = 0.0

  do icomp = 1, ws_ncomp

    sum   = sum   + w_flowsG  (icomp)
    sumpw = sumpw + w_flowsGpw(icomp)
    do ixw = 1, ws_nxw
      sumxw(ixw) = sumxw(ixw) + w_flowsGxw(icomp, ixw)
    enddo

    do jcmplg = 1, w_ncmplg
      do jdof = 1, ws_ndof
        sumxc(jcmplg, jdof) =  sumxc(jcmplg, jdof) &
                             + w_flowsGxc(icomp, jcmplg, jdof)
      enddo
    enddo

  enddo

  ! Setup flow part of residual (if isothermal, just ncomps-1 flows)

  do icompc = 1, nic
    rwbs  (icompc) = w_flowsG  (icompc)
    rwbspw(icompc) = w_flowsGpw(icompc)
    do ixw = 1, ws_nxw
      jwbs(icompc, ixw) = w_flowsGxw(icompc, ixw)
    enddo
    do jcmplg = 1, w_ncmplg
      do jdof = 1, ws_ndof
        rwbsxc(icompc, jcmplg, jdof) &
      = rwbsxc(icompc, jcmplg, jdof) + w_flowsGxc(icompc, jcmplg, jdof)
      enddo
    enddo
  enddo

  ! In thermal case, last equation is heat

  if (.not.ws_isothermal) then

    icomph = ws_nxw

    rwbs  (icomph) = w_flowsG  (eid)
    rwbspw(icomph) = w_flowsGpw(eid)

    do ixw = 1, ws_nxw
      jwbs(icomph, ixw) = w_flowsGxw(eid, ixw)
    enddo

    do jcmplg = 1, w_ncmplg
      do jdof = 1, ws_ndof
        rwbsxc(icomph, jcmplg, jdof) &
      = rwbsxc(icomph, jcmplg, jdof) + w_flowsGxc(eid, jcmplg, jdof)
      enddo
    enddo

  endif

  ! Update residuals with flow to/from surface term

  if (sum>0.0) then

  ! Is producing = subtract outflow equal to (sum times wellbore molar density)

    do icompc = 1, nic

      rwbs  (icompc) = rwbs  (icompc) -    sum  *w_zmf  (icompc)
      rwbspw(icompc) = rwbspw(icompc) - (  sumpw*w_zmf  (icompc) &
                                         + sum  *w_zmfpw(icompc) )
      do ixw = 1, ws_nxw
        jwbs(icompc, ixw) =  jwbs(icompc, ixw) &
                           - (  sumxw(ixw)*w_zmf  (icompc     ) &
                              + sum       *w_zmfxw(icompc, ixw) )
      enddo

      do jcmplg = 1, w_ncmplg
        do jdof = 1, ws_ndof
           rwbsxc(icompc, jcmplg, jdof) &
         = rwbsxc(icompc, jcmplg, jdof) - sumxc(jcmplg, jdof)*w_zmf(icompc)
        enddo
      enddo

    enddo

  ! Thermal part

    if (.not.ws_isothermal) then
      rwbs  (icomph) = rwbs  (icomph)-sum  *w_hpm
      rwbspw(icomph) = rwbspw(icomph)-sumpw*w_hpm   &
                                     -sum  *w_hpmdpw
      do ixw = 1, ws_nxw
        jwbs(icomph, ixw) = jwbs(icomph, ixw)-sumxw(ixw)*w_hpm         &
                                             -sum       *w_hpmdxw(ixw)
      enddo

      do jcmplg = 1, w_ncmplg
        do jdof = 1, ws_ndof
          rwbsxc(icomph, jcmplg, jdof) &
        = rwbsxc(icomph, jcmplg, jdof)-sumxc(jcmplg, jdof)*w_hpm
        enddo
      enddo

    endif
  else

    ! Is injecting, subtract outflow (add inflow)
    ! equal to sum times injection component molar density

    do icompc = 1, nic
      if (     (ws_isproducer    .and. icompc == ippref) &
          .or. (ws_isoilinjector .and. icompc == ipoil ) &
          .or. (ws_isgasinjector .and. icompc == ipgas ) &
          .or. (ws_iswatinjector .and. icompc == ipwat ) &
          .or. (ws_isslvinjector .and. icompc == ipslv )) then

        rwbs  (icompc) = rwbs  (icompc)-sum
        rwbspw(icompc) = rwbspw(icompc)-sumpw

        do ixw = 1, ws_nxw
          jwbs(icompc, ixw) = jwbs(icompc, ixw)-sumxw(ixw)
        enddo

        do jcmplg = 1, w_ncmplg
          do jdof = 1, ws_ndof
             rwbsxc(icompc, jcmplg, jdof) &
           = rwbsxc(icompc, jcmplg, jdof) - sumxc(jcmplg, jdof)
          enddo
        enddo

      endif
    enddo

    ! Thermal part

    if (.not.ws_isothermal) then

      rwbs  (icomph) = rwbs  (icomph)-sum  *ws_injection_h
      rwbspw(icomph) = rwbspw(icomph)-sumpw*ws_injection_h

      do ixw = 1, ws_nxw
        jwbs(icomph, ixw) = jwbs(icomph, ixw)-sumxw(ixw)*ws_injection_h
      enddo

      do jcmplg = 1, w_ncmplg
        do jdof = 1, ws_ndof
           rwbsxc(icomph, jcmplg, jdof) &
         = rwbsxc(icomph, jcmplg, jdof) - sumxc(jcmplg, jdof)*ws_injection_h
        enddo
      enddo

    endif
  endif

  ! Release work arrays

  deallocate(sumxw)
  deallocate(sumxc)

end subroutine getRwbsAndJwbs

! *************************************************************************** !

subroutine findWellboreGravityDensityPredictor()
  !
  ! Find predictor value for well gravity density (does not need Pw)
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscReal :: sumgdm, sumgds, summ, sums, sp, gdp, mp, wd
  PetscInt  :: icmpl, iphase, ierr
  PetscReal :: sl(FOUR_INTEGER)
  PetscReal :: sg(FOUR_INTEGER)

  wd = 0.0

  ! Form saturation and mobility weighted sums

  sumgdm = 0.0
  sumgds = 0.0

  summ = 0.0
  sums = 0.0

  do icmpl = 1, w_ncmpl
    do iphase = 1, ws_nphase
      sp  = c_sp  (icmpl, iphase)
      mp  = c_mob (icmpl, iphase)
      gdp = c_kgdp(icmpl, iphase)
      sumgdm = sumgdm+mp*gdp
      sumgds = sumgds+sp*gdp
      summ   = summ  +mp
      sums   = sums  +sp
    enddo
  enddo

  ! In case of cross-proc well continue summation over processors

  if (w_crossproc) then
    ierr = 0
    sg   = 0.0
    sl(1) = sumgdm
    sl(2) = sumgds
    sl(3) = summ
    sl(4) = sums
    call MPI_Allreduce(sl, sg, FOUR_INTEGER, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, w_comm, ierr)
    call checkErr(ierr)
    sumgdm = sg(1)
    sumgds = sg(2)
    summ   = sg(3)
    sums   = sg(4)
  endif

  ! Form density in Kg/rm3

  if (summ>0.0) then
    wd = sumgdm/summ
  else
    wd = sumgds/sums
  endif

  ! Store result with gravity constant included

  w_gd = wd*ws_gravity

end subroutine findWellboreGravityDensityPredictor

!  ************************************************************************** !

subroutine findWellborePropertiesAndDerivatives(pw, so, sg, sw, ss, &
                                                pb, t, option)
  !
  ! Find values and derivs of properties used to obtain well gravity density
  !
  ! Author: Dave Ponting
  ! Date  : 10/15/18

  use EOS_Water_module
  use EOS_Gas_module
  use EOS_Oil_module
  use EOS_Slv_module

  implicit none

  PetscReal, intent(in)  :: pw, so, sg, sw, ss, pb, t
  type(option_type), pointer :: option

  PetscReal :: mdos, mdost, mdospb
  PetscReal :: usf, usft, usfpw, usfpb
  PetscReal :: cr, crt, crpb
  PetscReal :: crusp, cruspt, crusppw, crusppb
  PetscReal :: mdo, mdot, mdopw, mdopb, hot, hop
  PetscReal :: uodum, uotdum, uopdum
  PetscReal :: mdg, mdgt, mdgp, hgt, hgp, ugdum, ugtdum, ugpdum
  PetscReal :: hs, hst, hsp, usdum, ustdum, uspdum
  PetscReal :: mdw, mdwp, mdwt, kgdw, hwp, hwt
  PetscReal :: mds, mdsp, mdst, z
  PetscReal :: ho, hg, hw

  PetscReal :: sp, mdp, hp
  PetscReal :: spX, mdpX, hpX

  PetscReal :: rsm, rsmt, rsmpb
  PetscReal :: den, deni, deni2, dent, denpb, denpw
  PetscReal :: xo, xg, xot, xopw, xopb, xgt, xgpw , xgpb
  PetscReal :: mwo, mwg, mws, mww, sum, sumi, sumpw

  PetscInt :: ierr, iphase, icomp, icoil, icgas
  PetscInt :: ipoil, ipgas, ipwat, ipslv, ixw

  PetscReal, allocatable :: sumxw(:)

  PetscInt, pointer :: table_idx(:)

  ! Initialise local variables

  ipoil = option%oil_phase
  ipgas = option%gas_phase
  ipwat = option%liquid_phase
  ipslv = option%solvent_phase

  ! Molar density of saturated oil

  mdos    = 0.0
  mdost   = 0.0
  mdospb  = 0.0

  ! Undersaturation

  usf     = 0.0
  usft    = 0.0
  usfpw   = 0.0
  usfpb   = 0.0

  ! Compressibility

  cr      = 0.0
  crt     = 0.0
  crpb    = 0.0

  ! Compressibility times undersaration

  crusp   = 0.0
  cruspt  = 0.0
  crusppw = 0.0
  crusppb = 0.0

  ! Molar and enthapy density of oil

  mdo     = 0.0
  mdot    = 0.0
  mdopw   = 0.0
  mdopb   = 0.0

  ho      = 0.0
  hg      = 0.0
  hw      = 0.0
  hs      = 0.0

  ! Molar density and enthalpy of gas

  mdg     = 0.0
  mdgt    = 0.0
  mdgp    = 0.0

  hg      = 0.0
  hgt     = 0.0
  hgp     = 0.0

  ugdum   = 0.0
  ugtdum  = 0.0
  ugpdum  = 0.0

  ! Set up work array

  allocate(sumxw(ws_ncomp))
  sumxw = 0.0

  ! Set up dummy table_idx

  table_idx => null()

  if (option%num_table_indices > 0) then
    allocate(table_idx(option%num_table_indices))
    table_idx = 1
  endif

  ! Get oil and gas molecular weights

  mwo = EOSOilGetFMW()
  mwg = EOSGasGetFMW()
  if (is_solvent) then
    mws = EOSSlvGetFMW()
  else
    mws = 0.0
  endif
  mww = FMWH2O

  ! Default oil composition (as mole fractions)

  xo = 1.0
  xg = 0.0

  xot  = 0.0
  xopb = 0.0
  xopw = 0.0

  xgt  = 0.0
  xgpb = 0.0
  xgpw = 0.0

  den   = 0.0
  dent  = 0.0
  denpb = 0.0
  denpw = 0.0

  ! Initialise the output properties and derivatives

  ! Saturations
  w_sp     = 0.0
  w_spxw   = 0.0

  ! Molar density of each phase
  w_mdp    = 0.0
  w_mdppw  = 0.0
  w_mdpxw  = 0.0

  ! Total molar density of each component
  w_mdc    = 0.0
  w_mdcpw  = 0.0
  w_mdcxw  = 0.0

  ! Total mole fractions
  w_zmf    = 0.0
  w_zmfpw  = 0.0
  w_zmfxw  = 0.0

  ! Mass density of each phase
  w_kgdp   = 0.0
  w_kgdppw = 0.0
  w_kgdpxw = 0.0

  ! Molar enthaply of each phase
  w_hp     = 0.0
  w_hpdpw  = 0.0
  w_hpdxw  = 0.0

  ! Total enthalpy density
  w_h      = 0.0
  w_hdpw   = 0.0
  w_hdxw   = 0.0

  ! Molar enthalpy density
  w_hpm    = 0.0
  w_hpmdpw = 0.0
  w_hpmdxw = 0.0

  ! Total gravity density
  w_gd     = 0.0
  w_gdpw   = 0.0
  w_gdxw   = 0.0

  ! Set up saturations and derivatives

  if (ws_oil) then
    w_sp  (ipoil) = so
    w_spxw(ipoil, ixwso) =  1.0
    w_spxw(ipwat, ixwso) = -1.0
  endif

  if (ws_gas) then
    w_sp(ipgas) = sg
    if (ixwsg>0) then
      w_spxw(ipgas, ixwsg) =  1.0
      w_spxw(ipwat, ixwsg) = -1.0
    endif
  endif

  w_sp(ipwat) = sw

  if (is_solvent) then
    w_sp  (ipslv       ) =  ss
    w_spxw(ipslv, ixwss) =  1.0
    w_spxw(ipwat, ixwss) = -1.0
  endif

  ! Get oil phase molar density and phase molar enthalpy

  if (ws_oil) then

    call EOSOilDensityEnergy(t    , pw    , mdo   , mdot, mdopw, &
                             ho   , hot   , hop   , &
                             uodum, uotdum, uopdum, ierr, table_idx)
    ho  = ho  * 1.d-6 ! J/kmol -> MJ/kmol
    hot = hot * 1.d-6
    hop = hop * 1.d-6

    ! Special case of black oil
    if (ws_is_black_oil) then

      call EOSOilRS(t, pb, rsm, rsmt, rsmpb, ierr, table_idx)
      den = 1.0+rsm
      deni = 0.0
      if (den /= 0.0) then
        deni = 1.0/den
      endif
      deni2 = deni*deni

      dent  = rsmt
      denpb = rsmpb
      denpw = 0.0

      xo =     deni
      xg = rsm*deni

      xot  = -deni2*dent
      xopb = -deni2*denpb
      xopw = 0.0

      xgt  = (rsmt -xg*dent )*deni
      xgpb = (rsmpb-xg*denpb)*deni
      xgpw = 0.0

      if (w_issat) then
        xopw  = xopb
        xgpw  = xgpb
        denpw = denpb
        xopb  = 0.0
        xgpb  = 0.0
        denpb = 0.0
      endif

      call EOSOilDensity        (t, pb, mdos, mdost, mdospb, ierr, table_idx)
      call EOSOilCompressibility(t, pb, cr  , crt  , crpb  , ierr, table_idx)

      if (w_issat) then

        mdo   = mdos  *den
        mdot  = mdost *den+mdos*dent
        mdopw = mdospb*den+mdos*denpw

      else

        crusp   =  cr *(pw-pb)
        cruspt  =  crt*(pw-pb)
        crusppw =  cr
        crusppb = -cr+crpb*(pw-pb)

        usf   =  1.0+crusp*(1.0+0.5*crusp)
        usft  = (1.0+crusp)*cruspt
        usfpw = (1.0+crusp)*crusppw
        usfpb = (1.0+crusp)*crusppb

        mdo   =  mdos  *usf  *den
        mdot  =  mdost *usf  *den   &
               + mdos  *usft *den   &
               + mdos  *usf  *dent
        mdopw =  mdos  *usfpw*den   &
               + mdos  *usf  *denpw
        mdopb =  mdospb*usf  *den   &
               + mdos  *usfpb*den   &
               + mdos  *usf  *denpb

      endif

    endif

    ! Store oil molar density and derivatives

    if (ws_is_black_oil) then

                   w_mdp   (ipoil)        = mdo
                   w_mdppw (ipoil)        = mdopw
      if (ixwt >0) w_mdpxw (ipoil, ixwt ) = mdot
      if (ixwpb>0) w_mdpxw (ipoil, ixwpb) = mdopb

                   w_kgdp  (ipoil)        = mdo  *( xo  *mwo + xg  *mwg)
                   w_kgdppw(ipoil)        = mdopw*( xo  *mwo + xg  *mwg) &
                                           +mdo  *( xopw*mwo + xgpw*mwg)
      if (ixwt >0) w_kgdpxw(ipoil, ixwt ) = mdot *( xo  *mwo + xg  *mwg) &
                                           +mdo  *( xot *mwo + xgt *mwg)
      if (ixwpb>0) w_kgdpxw(ipoil, ixwpb) = mdopb*( xo  *mwo + xg  *mwg) &
                                           +mdo  *( xopb*mwo + xgpb*mwg)
    else
                   w_mdp   (ipoil)        = mdo
                   w_mdppw (ipoil)        = mdopw
      if (ixwt >0) w_mdpxw (ipoil, ixwt ) = mdot

                   w_kgdp  (ipoil)         = mdo  *mwo
                   w_kgdppw(ipoil)         = mdopw*mwo
      if (ixwt >0) w_kgdpxw(ipoil, ixwt )  = mdot *mwo

    endif

                 w_hp   (ipoil)       = ho
                 w_hpdpw(ipoil)       = hop
    if (ixwt >0) w_hpdxw(ipoil, ixwt) = hot

  endif

  ! Gas molar density and derivatives

  if (ws_gas) then

    call EOSGasDensityEnergy( t, pw, mdg  , mdgt  , mdgp &
                                   , hg   , hgt   , hgp  &
                                   , ugdum, ugtdum, ugpdum, ierr )
    hg  = hg  * 1.d-6 ! J/kmol -> MJ/kmol
    hgt = hgt * 1.d-6
    hgp = hgp * 1.d-6

                w_mdp   (ipgas)       = mdg
                w_mdppw (ipgas)       = mdgp
    if (ixwt>0) w_mdpxw (ipgas, ixwt) = mdgt

                w_kgdp  (ipgas)       = mdg *mwg
                w_kgdppw(ipgas)       = mdgp*mwg
    if (ixwt>0) w_kgdpxw(ipgas, ixwt) = mdgt*mwg

                w_hp   (ipgas)        = hg
                w_hpdpw(ipgas)        = hgp
    if (ixwt>0) w_hpdxw(ipgas, ixwt)  = hgt

! Correct the oil phase enphalpy for dissolved gas

    if (ws_is_black_oil) then
                   w_hp   (ipoil)        = ho *xo   + hg *xg
                   w_hpdpw(ipoil)        = hop*xo   + hgp*xg   &
                                          +ho *xopw + hg *xgpw
      if (ixwt>0)  w_hpdxw(ipoil, ixwt)  = hot*xo   + hgt*xg   &
                                          +ho *xot  + hg *xgt
      if (ixwpb>0) w_hpdxw(ipoil, ixwpb) = ho *xopb + hg *xgpb
    endif

  endif

  ! Water molar density and derivatives

  if (ws_wat) then

    call EOSWaterDensity    (t, pw, kgdw, mdw, mdwp, mdwt, ierr)
    call EOSWaterEnthalpy   (t, pw, hw, hwp, hwt         , ierr)
    hw  = hw  * 1.d-6 ! J/kmol -> MJ/kmol
    hwt = hwt * 1.d-6
    hwp = hwp * 1.d-6

                w_mdp   (ipwat)       = mdw
                w_mdppw (ipwat)       = mdwp
    if (ixwt>0) w_mdpxw (ipwat, ixwt) = mdwt

                w_kgdp  (ipwat)       = mdw *mww
                w_kgdppw(ipwat)       = mdwp*mww
    if (ixwt>0) w_kgdpxw(ipwat, ixwt) = mdwt*mww

                w_hp   (ipwat)        = hw
                w_hpdpw(ipwat)        = hwp
    if (ixwt>0) w_hpdxw(ipwat, ixwt)  = hwt

  endif

  ! Solvent molar density and derivatives

  if (is_solvent) then

    call EOSSlvDensity(t, pw, mds, mdst, mdsp, ierr, table_idx)
    call EOSSlvEnergy (t, pw, hs, hst, hsp, usdum, ustdum, uspdum, ierr)

    hs  = hs  * 1.d-6 ! J/kmol -> MJ/kmol
    hst = hst * 1.d-6
    hsp = hsp * 1.d-6

                w_mdp   (ipslv)       = mds
                w_mdppw (ipslv)       = mdsp
    if (ixwt>0) w_mdpxw (ipslv, ixwt) = mdst

                w_kgdp  (ipslv)       = mds *mws
                w_kgdppw(ipslv)       = mdsp*mws
    if (ixwt>0) w_kgdpxw(ipslv, ixwt) = mdst*mws

                w_hp    (ipslv)       = hs
                w_hpdpw (ipslv)       = hsp
    if (ixwt>0) w_hpdxw (ipslv, ixwt) = hst

  endif

  ! Form total well molar densities

  w_mdc   = 0.0
  w_mdcpw = 0.0
  w_mdcxw = 0.0

  if (ws_is_black_oil) then
    ! Black-oil - oil contains both oil and gas compositions
    do iphase = 1, ws_nphase
      sp  = w_sp (iphase)
      mdp = w_mdp(iphase)
      if (iphase == option%oil_phase) then
        icoil = option%oil_phase
        icgas = option%gas_phase
        w_mdc  (icoil) = w_mdc  (icoil)+xo*sp*mdp
        w_mdc  (icgas) = w_mdc  (icgas)+xg*sp*mdp
        w_mdcpw(icoil) = w_mdcpw(icoil)+xo*sp*w_mdppw(iphase)
        w_mdcpw(icgas) = w_mdcpw(icgas)+xg*sp*w_mdppw(iphase)
        do ixw = 1, ws_nxw
          w_mdcxw(icoil, ixw) = w_mdcxw(icoil, ixw) &
                               + xo*(  w_spxw(iphase, ixw)*mdp &
                                     + sp*w_mdpxw(iphase, ixw) )
          w_mdcxw(icgas, ixw) = w_mdcxw(icgas, ixw) &
                               + xg*(  w_spxw(iphase, ixw)*mdp &
                                     + sp*w_mdpxw(iphase, ixw) )
          if (ixw == ixwpb) then
            w_mdcxw(icoil, ixw) = w_mdcxw(icoil, ixw) + xopb*sp*mdp
            w_mdcxw(icgas, ixw) = w_mdcxw(icgas, ixw) + xgpb*sp*mdp
          endif
          if (ixw == ixwt) then
            w_mdcxw(icoil, ixw) = w_mdcxw(icoil, ixw) + xot*w_sp(iphase)*mdp
            w_mdcxw(icgas, ixw) = w_mdcxw(icgas, ixw) + xgt*w_sp(iphase)*mdp
          endif
        enddo
      else
        icomp = iphase
        w_mdc  (icomp) = w_mdc  (icomp)+sp*mdp
        w_mdcpw(icomp) = w_mdcpw(icomp)+sp*w_mdppw(iphase)

        do ixw = 1, ws_nxw
          w_mdcxw(icomp, ixw) = w_mdcxw(icomp, ixw) &
                               + sp*w_mdpxw(iphase, ixw) &
                               + w_spxw(iphase, ixw)*mdp
        enddo
      endif
    enddo
  else
    ! Non black-oil - each phase contains just its own composition
    do iphase = 1, ws_nphase
      icomp = iphase
      w_mdc  (icomp) = w_mdc  (icomp)+w_sp(iphase)*w_mdp  (iphase)
      w_mdcpw(icomp) = w_mdcpw(icomp)+w_sp(iphase)*w_mdppw(iphase)

      do ixw = 1, ws_nxw
        w_mdcxw(icomp, ixw) =   w_mdcxw(icomp, ixw) &
                              + w_spxw(iphase, ixw)*w_mdp  (iphase     ) &
                              + w_sp  (iphase     )*w_mdpxw(iphase, ixw)
      enddo
    enddo
  endif

  ! Build up wellbore mole fractions

  sum   = 0.0
  sumi  = 0.0
  sumpw = 0.0
  sumxw = 0.0

  do icomp = 1, ws_ncomp
    sum   = sum  + w_mdc  (icomp)
    sumpw = sumpw+ w_mdcpw(icomp)
    do ixw = 1, ws_nxw
      sumxw(ixw) = sumxw(ixw)+ w_mdcxw(icomp, ixw)
    enddo
  enddo

  if (sum>0.0) then
    sumi = 1.0/sum
    do icomp = 1, ws_ncomp
      z = w_mdc(icomp)*sumi
      w_zmf  (icomp) =  z
      w_zmfpw(icomp) = (w_mdcpw(icomp)-z*sumpw)*sumi
      do ixw = 1, ws_nxw
        w_zmfxw(icomp, ixw) = (w_mdcxw(icomp, ixw)-z*sumxw(ixw))*sumi
      enddo
    enddo
  endif

  ! Build up wellbore saturation*density

  do iphase = 1, ws_nphase

    w_gd   = w_gd   + w_sp (iphase)*w_kgdp  (iphase)*ws_gravity
    w_gdpw = w_gdpw + w_sp (iphase)*w_kgdppw(iphase)*ws_gravity
    do ixw = 1, ws_nxw
      w_gdxw(ixw) = w_gdxw(ixw) &
                   + ( w_spxw(iphase, ixw)*w_kgdp  (iphase     ) &
                      +w_sp  (iphase     )*w_kgdpxw(iphase, ixw) )*ws_gravity
    enddo

  enddo

  ! Build up enthalpy per unit reservoir volume

  do iphase = 1, ws_nphase

    sp  = w_sp (iphase)
    hp  = w_hp (iphase)
    mdp = w_mdp(iphase)

    w_h    = w_h    + sp*mdp*hp
    w_hdpw = w_hdpw + sp*( w_mdppw(iphase)*hp              &
                          +  mdp          *w_hpdpw(iphase) )

    do ixw = 1, ws_nxw

      spX  = w_spxw (iphase, ixw)
      hpX  = w_hpdxw(iphase, ixw)
      mdpX = w_mdpxw(iphase, ixw)

      w_hdxw(ixw) = w_hdxw(ixw) + spX*mdp*hp + sp*mdpX*hp + sp*mdp*hpX
    enddo

  enddo

  ! Build up enthalpy per mole as (enthalpy/res.vol.)/(moles/res.vol.)

  w_hpm    =  w_h                *sumi
  w_hpmdpw = (w_hdpw-w_hpm*sumpw)*sumi
  do ixw = 1, ws_nxw
    w_hpmdxw(ixw) = (w_hdxw(ixw)-w_hpm*sumxw(ixw))*sumi
  enddo

  ! Free the table index array and sumxw arrays

  if (option%num_table_indices > 0) then
    deallocate(table_idx)
  endif

  deallocate(sumxw)

end subroutine findWellborePropertiesAndDerivatives

! *************************************************************************** !

subroutine wellSolverLoaderTOWG(aux, option)
  !
  ! Load completion solution arrays for TOWG mode
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Auxiliary_module
  use PM_TOWG_Aux_module

  implicit none

  type(auxiliary_type) :: aux
  type(option_type), pointer :: option

  class(pm_towg_aux_type), pointer :: TOWG

  PetscInt :: icmpl, ghosted_id

  ws_xpr = TOWG_OIL_PRESSURE_DOF

  towg => aux%TOWG

  do icmpl = 1, w_ncmpl
    ghosted_id = c_ghosted_id(icmpl)
    call loadCellDataTOWG(towg%auxvars(ZERO_INTEGER, ghosted_id), &
                          option, icmpl)
  enddo

end subroutine wellSolverLoaderTOWG

! *************************************************************************** !

subroutine wellSolverLoaderTOIL(aux, option)
  !
  ! Well solver loader for TOIl mode
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Auxiliary_module
  use PM_TOilIms_Aux_module

  implicit none

  type(auxiliary_type) :: aux
  type(option_type), pointer :: option

  class(pm_toil_ims_aux_type), pointer :: toil

  PetscInt :: icmpl, ghosted_id

  ws_xpr = TOIL_IMS_PRESSURE_DOF

  toil => aux%TOil_ims

  do icmpl = 1, w_ncmpl
    ghosted_id = c_ghosted_id(icmpl)
    call loadCellDataTOIL(toil%auxvars(ZERO_INTEGER, ghosted_id), &
                          option, icmpl)
  enddo

end subroutine wellSolverLoaderTOIL

! *************************************************************************** !

subroutine loadCellDataTOWG(auxvar, option, icmpl)
  !
  ! Load completion cell arrays for TOWG mode (for one cell)
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  use AuxVars_TOWG_module

  implicit none

  class(auxvar_towg_type) :: auxvar
  type(option_type), pointer :: option
  PetscInt, intent(in) :: icmpl
  PetscInt :: iphase, idof

  c_p (icmpl) = auxvar%pres(option%oil_phase)
  c_t (icmpl) = auxvar%temp
  if (ws_is_black_oil) then
    c_pb(icmpl) = auxvar%bo%bubble_point
  endif

  do iphase = 1, ws_nphase

    c_mdp (icmpl, iphase) = auxvar%den     (iphase)
    c_sp  (icmpl, iphase) = auxvar%sat     (iphase)
    c_mob (icmpl, iphase) = auxvar%mobility(iphase)
    c_kgdp(icmpl, iphase) = auxvar%den_kg  (iphase)
    c_hp  (icmpl, iphase) = auxvar%H       (iphase)

    if (.not.option%flow%numerical_derivatives) then
      do idof = 1, ws_ndof
        c_mdpX (icmpl, iphase, idof) = auxvar%D_den     (iphase, idof)
        c_spX  (icmpl, iphase, idof) = auxvar%D_sat     (iphase, idof)
        c_mobX (icmpl, iphase, idof) = auxvar%D_mobility(iphase, idof)
        c_kgdpX(icmpl, iphase, idof) = auxvar%D_den_kg  (iphase, idof)
        c_hpX  (icmpl, iphase, idof) = auxvar%D_H       (iphase, idof)
        if (ws_isothermal .and. (idof == towg_energy_dof)) then
          c_mdpX (icmpl, iphase, idof) = 0.0
          c_spX  (icmpl, iphase, idof) = 0.0
          c_mobX (icmpl, iphase, idof) = 0.0
          c_kgdpX(icmpl, iphase, idof) = 0.0
          c_hpX  (icmpl, iphase, idof) = 0.0
        endif
      enddo
    endif

  enddo

  if (ws_is_black_oil) then

    c_xo(icmpl) = auxvar%bo%xo
    c_xg(icmpl) = auxvar%bo%xg

    if (.not.option%flow%numerical_derivatives) then
      do idof = 1, ws_ndof
        c_xoX(icmpl, idof) = auxvar%bo%D_xo(idof)
        c_xgX(icmpl, idof) = auxvar%bo%D_xg(idof)
        if (ws_isothermal .and. (idof == towg_energy_dof)) then
          c_xoX(icmpl, idof) = 0.0
          c_xgX(icmpl, idof) = 0.0
        endif
      enddo
    endif

  endif

end subroutine loadCellDataTOWG

! *************************************************************************** !

subroutine loadCellDataTOIL(auxvar, option, icmpl)
  !
  ! Load completion cell arrays for TOIL mode (for one cell)
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  use Material_Aux_class
  use AuxVars_TOilIms_module
  use PM_TOilIms_Aux_module

  implicit none

  class(auxvar_toil_ims_type) :: auxvar
  type(option_type), pointer :: option
  PetscInt, intent(in) :: icmpl
  PetscInt :: iphase, idof

  c_p (icmpl) = auxvar%pres(option%oil_phase)
  c_t (icmpl) = auxvar%temp

  do iphase = 1, ws_nphase

    c_mdp (icmpl, iphase) = auxvar%den     (iphase)
    c_sp  (icmpl, iphase) = auxvar%sat     (iphase)
    c_mob (icmpl, iphase) = auxvar%mobility(iphase)
    c_kgdp(icmpl, iphase) = auxvar%den_kg  (iphase)
    c_hp  (icmpl, iphase) = auxvar%H       (iphase)

    if (.not.option%flow%numerical_derivatives) then
      do idof = 1, ws_ndof
        c_mdpX (icmpl, iphase, idof) = auxvar%D_den     (iphase, idof)
        c_spX  (icmpl, iphase, idof) = auxvar%D_sat     (iphase, idof)
        c_mobX (icmpl, iphase, idof) = auxvar%D_mobility(iphase, idof)
        c_kgdpX(icmpl, iphase, idof) = auxvar%D_den_kg  (iphase, idof)
        c_hpX  (icmpl, iphase, idof) = auxvar%D_H       (iphase, idof)
        if (ws_isothermal .and. (idof == TOIL_IMS_ENERGY_DOF)) then
          c_mdpX (icmpl, iphase, idof) = 0.0
          c_spX  (icmpl, iphase, idof) = 0.0
          c_mobX (icmpl, iphase, idof) = 0.0
          c_kgdpX(icmpl, iphase, idof) = 0.0
          c_hpX  (icmpl, iphase, idof) = 0.0
        endif
      enddo
    endif

  enddo

end subroutine loadCellDataTOIL

! *************************************************************************** !

subroutine invertJacobian(deti)
  !
  ! Invert the wellbore composition Jacobian
  !
  ! Author: Dave Ponting
  ! Date  : 09/19/18

  implicit none

  PetscReal, intent(out) :: deti
  PetscReal :: det, a, b, c, d, ai, bi, ci, di
  PetscReal :: cf(3, 3), eps
  PetscInt  :: i, j

  eps = 1.0E-6

  if (ws_nxw == 1) then
    det  = jwbs(1, 1)
    deti = 0.0
    if (abs(det)>0.0) deti = 1.0/det
    jwbsi(1, 1) = deti
  endif

  if (ws_nxw == 2) then

    a = jwbs(1, 1)
    b = jwbs(1, 2)
    c = jwbs(2, 1)
    d = jwbs(2, 2)

    det = a*d-b*c

    deti = 0.0
    if (abs(det)>0.0) deti = 1.0/det

    ai =  d*deti
    bi = -b*deti
    ci = -c*deti
    di =  a*deti

    jwbsi(1, 1) = ai
    jwbsi(1, 2) = bi
    jwbsi(2, 1) = ci
    jwbsi(2, 2) = di

  endif

  if (ws_nxw == 3) then

     cf(1, 1) =  jwbs(2, 2)*jwbs(3, 3)-jwbs(3, 2)*jwbs(2, 3)
     cf(1, 2) = -jwbs(2, 1)*jwbs(3, 3)+jwbs(3, 1)*jwbs(2, 3)
     cf(1, 3) =  jwbs(2, 1)*jwbs(3, 2)-jwbs(3, 1)*jwbs(2, 2)

     cf(2, 1) = -jwbs(1, 2)*jwbs(3, 3)+jwbs(3, 2)*jwbs(1, 3)
     cf(2, 2) =  jwbs(1, 1)*jwbs(3, 3)-jwbs(3, 1)*jwbs(1, 3)
     cf(2, 3) = -jwbs(1, 1)*jwbs(3, 2)+jwbs(3, 1)*jwbs(1, 2)

     cf(3, 1) =  jwbs(1, 2)*jwbs(2, 3)-jwbs(2, 2)*jwbs(1, 3)
     cf(3, 2) = -jwbs(1, 1)*jwbs(2, 3)+jwbs(2, 1)*jwbs(1, 3)
     cf(3, 3) =  jwbs(1, 1)*jwbs(2, 2)-jwbs(2, 1)*jwbs(1, 2)

     det = jwbs(1, 1)*cf(1, 1)+jwbs(1, 2)*cf(1, 2)+jwbs(1, 3)*cf(1, 3)

     deti = 0.0
     if (abs(det) > 0.0) deti = 1.0/det

     do i = 1, 3
       do j = 1, 3
         jwbsi(i, j) = cf(j, i)*deti
       enddo
     enddo

  endif

  if (ws_nxw > 3) then
    det = invertGauss(jwbs, jwbsi, ws_nxw)
  endif

end subroutine invertJacobian

! ************************************************************************** !

function invertGauss(j, jinv, n)

  !
  ! Invert an nxn matrix using Gauss elimination with pivoting
  !
  ! If ei is the unit matrix with one entry at location i, then Ainv.ei=ci
  ! where ci is the ith column of Ainv. So if we solve A.ci=ei for ci
  ! get ith column of Ainv.
  !
  ! In the following, initialise the matrix e to a set of n columns, such that
  ! column ic contains ei, so e is the unit matrix.

  ! We solve A.ci=ei for all the columns, and at the end of the elimination
  ! e contains Ainv.
  ! This can also be viewed as a set of matrix pre-multiplications on both
  ! A and E, such that after the transformation T, A has been transformed
  ! to the unit matrix, and E to Ainv (if T.A=E, T must be Ainv, so T.E=Ainv)
  ! A[ir][ic] is row iw, column ic.
  !
  ! Returns minimum pivot (will be zero if determinant is zero)
  !
  ! Author: Dave Ponting
  ! Date  : 01/25/19

  implicit none

  PetscReal :: invertGauss
  PetscReal, intent(in)  :: j(:,:)
  PetscReal, intent(out) :: jinv(:,:)
  PetscInt , intent(in)  :: n

  PetscReal :: det
  PetscBool :: firstdet
  PetscReal, allocatable :: a(:,:)
  PetscReal, allocatable :: e(:,:)
  PetscInt, allocatable  :: ip(:)

  PetscInt :: ir, ic, iclast, iclastp, irbest &
             , irbestp, irp, irap, irbp, icp, jcp, jc
  PetscReal :: fabsFinal, best, fval, pivot, pivotinv, fabsPivot, v, d, dinv

  allocate(a(n, n))
  allocate(e(n, n))
  allocate(ip(n))

  det = 0.0
  firstdet = PETSC_TRUE
  a = 0.0
  e = 0.0

  ! Store working matrix a and unit rhs matrix, set up pivot matrix.
  ! The actual location of the irth row will be ip(ir)

  do ir = 1, n
    do ic = 1, n
      a(ir, ic) = j(ir, ic)
      if (ir == ic) e(ir, ic) = 1.0
    enddo
    ip(ir) = ir
  enddo

  ! Pivot by column into upper triangular form - no need to do the last column

  do ic = 1, n-1

    ! Find the row with max absolute value in this column
    ! starting from diagonal element

    irbest  = ic
    irbestp = ip(irbest)
    best = abs(a(irbestp, ic))

    ! Loop down the column looking for better pivot rows

    do ir = ic+1, n
      irp = ip(ir)
      fval = abs(a(irp, ic))
      if (fval > best) then
        irbest = ir
        best   = fval
      endif
    enddo

    ! Check if irbest is not ic, pivot if so

    if (irbest /= ic) then
      ! Find current pointers to these rows and swap them
      irap = ip(irbest)
      irbp = ip(ic)
      ip(irbest) = irbp
      ip(ic    ) = irap
    endif

    ! Now do the actual elimination, zeroing
    ! entries in column ic below the diagonal

    icp = ip(ic)

    ! Find the inverse pivot, keeping track of the smallest pivot
    ! to detect zero determinants

    pivot = a(icp, ic)
    fabspivot = abs(pivot)
    if (firstdet) then
      det = fabsPivot
      firstdet = PETSC_FALSE
    else
      det = min(det, fabsPivot)
    endif

    pivotinv = 0.0
    if (fabsPivot > 0.0) pivotinv = 1.0/pivot

    ! Go down the column, subtracting value pivotinv
    ! times equation ic from equation ir

    do ir = ic+1, n
      irp = ip(ir)
      v = a(irp, ic)*pivotinv

      ! Carry out pivoting operation on row of matrix
      ! and all the right hand sides

      do jc = 1, n
        a(irp, jc) = a(irp, jc)-v*a(icp, jc)
        e(irp, jc) = e(irp, jc)-v*e(icp, jc)
      enddo
    enddo
  enddo

  ! Include final equation in the determinant check
  ! (no need to set firstDet now as will not be used again)

  iclast  =    n
  iclastp = ip(n)
  fabsFinal = abs(a(iclastp, iclast))
  if (firstdet) then
    det =          fabsFinal
  else
    det = min(det, fabsFinal)
  endif

  ! Back-substitute each of the right hand sides

  do ic = 1, n

    ! For this rhs, back substitute by row in reverse order

    do ir = n, 1, -1

      irp = ip(ir)

      ! Use previously evaluated solution values

      if (ir<n) then
        do jc = ir+1, n
          jcp = ip(jc)
          e(irp, ic) = e(irp, ic)-a(irp, jc)*e(jcp, ic)
        enddo
      endif

      ! And finally establish value using inverse diagonal

      d    = a(irp, ir)
      dinv = 0.0
      if (abs(d) > 0.0) dinv = 1.0/d
      e(irp, ic) = e(irp, ic)*dinv
    enddo
  enddo

  ! Put inverse matrix back into Jinv

  do ir = 1, n
    do ic = 1, n
      Jinv(ir, ic) = e(ip(ir), ic)
    enddo
  enddo

  invertGauss = det

  deallocate(a)
  deallocate(e)
  deallocate(ip)

end function invertGauss

! *************************************************************************** !

subroutine findInjectionEnthalpy()
  !
  ! Find the molar enthalpy of the injection fluid
  !
  ! Author: Dave Ponting
  ! Date  : 10/16/18

  use EOS_Oil_module
  use EOS_Gas_module
  use EOS_Water_module
  use EOS_Slv_module

  implicit none

  PetscReal :: p, t, u, h
  PetscInt  :: ierr

  ws_injection_h = 0.0

  p = ws_injection_p
  t = ws_injection_t

  ierr = 0
  u    = 0.0
  h    = 0.0

  if (w_type == OIL_INJ_WELL_TYPE) then
    call EOSOilEnthalpy(t , p, h, ierr)
  endif

  if (w_type == GAS_INJ_WELL_TYPE) then
    call EOSGasEnergy(t, p, h, u, ierr)
  endif

  if (w_type == WAT_INJ_WELL_TYPE) then
    call EOSWaterEnthalpy(t, p, h, ierr)
  endif

  if (w_type == SLV_INJ_WELL_TYPE) then
    call EOSSlvEnergy(t, p, h, u, ierr)
  endif

  ws_injection_h = h*1.0d-6

end subroutine findInjectionEnthalpy

! ************************************************************************** !

subroutine incrTactX(ip, c, tact, tactpw)
  !
  ! Build up the actual flow
  !
  ! Author: Dave Ponting
  ! Date  : 10/16/18

  implicit none

  PetscInt , intent(in) :: ip
  PetscReal, intent(in) :: c
  PetscReal, intent(inout) :: tact
  PetscReal, intent(inout) :: tactpw

  PetscInt :: icmplg, idof

  ! Actual flow

  tact = w_flowsG(ip)*c

  ! Derivs wrt Pw (expanding out the Xw dependence)

  ! Basic term

  tactpw = w_flowsGpw(ip)*c

  ! Derivs wrt Xc (expanding out the Xw dependence)

  do icmplg = 1, w_ncmplg
    do idof = 1, ws_ndof
      w_tactxc(icmplg, idof) = &
      w_tactxc(icmplg, idof) + w_flowsGxc(ip, icmplg, idof)*c
    enddo
  enddo

end subroutine incrTactX

! *************************************************************************** !

subroutine findFullFlowDerivatives(jwwi)
  !
  ! Find fully expanded flow derivatives (ie total derivatives in Xp)
  ! based on :
  !
  ! df/dPw, dfw/dXw and df/dXc from flow calculation
  !
  ! dXw/dPw and dXw/dXc from Xw-solution
  !             dPw/dXc from Pw solution (formed in this routine)
  !
  ! Author: Dave Ponting
  ! Date  : 10/16/18

  implicit none

  PetscReal, intent(in) :: jwwi

  PetscInt  :: icomp, icmpl, jcmplg, jdof
  PetscReal :: sum

  ! Step 1: construct dPw/dXc=-(1/(dRw/dPw))*(dRw/dXc)
  !         Note derivatives have expanded into Pw and Xc derivatives

  do jcmplg = 1, w_ncmpl
    do jdof = 1, ws_ndof
      w_dpwdxc(jcmplg, jdof) = -jwwi*w_rwxc(jcmplg, jdof)
    enddo
  enddo

  ! Step 2: expand flow derivatives in Pw and Xw into Xc

  ! Loop over flows
  do icmpl = 1, w_ncmpl
    do icomp = 1, ws_ncomp

      ! Loop over Xc
      do jcmplg = 1, w_ncmplg
        do jdof = 1, ws_ndof

          ! Extend df/dXc with df/fPw.dPw/dXc

          sum = c_flowspw(icmpl, icomp)*w_dpwdxc(jcmplg, jdof)

          ! Add extension to this flow term

            c_flowsxc(icmpl, icomp, jcmplg, jdof) &
          = c_flowsxc(icmpl, icomp, jcmplg, jdof)+sum

        enddo
      enddo
    enddo
  enddo

end subroutine findFullFlowDerivatives

! *************************************************************************** !

subroutine checkSurfaceDensities(option)
  !
  ! Check that a valid surface density calculation exists for the well model
  !
  ! Author: Dave Ponting
  ! Date  : 01/23/19

  use EOS_Oil_module
  use EOS_Gas_module
  use EOS_Water_module
  use EOS_Slv_module

  implicit none

  type (option_type), pointer :: option

  PetscInt  :: iphase, ierr
  PetscReal :: mw, sd

  ierr = 0

  do iphase = 1, option%nphase

    mw = 1.0
    sd = 1.0

    if (iphase == option%oil_phase) then
      mw = EOSOilGetFMW()
      if (mw<0.0) then
        option%io_buffer = 'Oil molecular weight not set';ierr = 1
      endif
      sd = EOSOilGetSurfaceDensity()
      if (sd<0.0) then
        option%io_buffer = 'Oil surface reference density not set';ierr = 1
      endif
    endif
    if (iphase == option%gas_phase) then
      mw = EOSGasGetFMW()
      if (mw<0.0) then
        option%io_buffer = 'Gas molecular weight not set';ierr = 1
      endif
      sd = EOSGasGetSurfaceDensity()
      if (sd<0.0) then
        option%io_buffer = 'Gas surface reference density not set';ierr = 1
      endif
    endif
    if (iphase == option%liquid_phase) then
      mw = FMWH2O
      if (mw<0.0) then
        option%io_buffer = 'Water molecular weight not set';ierr = 1
      endif
      sd = EOSWaterGetSurfaceDensity()
      if (sd<0.0) then
        option%io_buffer = 'Water surface reference density not set';ierr = 1
      endif
     endif
    if (iphase == option%solvent_phase) then
      mw = EOSSlvGetFMW()
      if (mw<0.0) then
        option%io_buffer = 'Solvent molecular weight not set';ierr = 1
      endif
      sd = EOSSlvGetSurfaceDensity()
      if (sd<0.0) then
        option%io_buffer = 'Solvent surface reference density not set';ierr = 1
      endif
    endif

  enddo

  if (ierr == 1) then
    call printErrMsg(option)
  endif

end subroutine checkSurfaceDensities

! *************************************************************************** !

subroutine checkSolverAvailable(option)
  !
  ! Check that a well solver exists for this run
  !
  ! Author: Dave Ponting
  ! Date  : 01/23/19

  use Grid_Grdecl_module, only : GetIsGrdecl

  implicit none

  type (option_type), pointer :: option
  PetscInt  :: ncmp, iphase, nphase, wsnx
  PetscBool :: water_found, mode_ok, is_grdecl

  ncmp   = option%nphase+1
  nphase = option%nphase

  ! Do we have a direct matrix inverter for this system?

  if (ws_isothermal) then
    wsnx = ncmp-2
  else
    wsnx = ncmp-1
  endif
  if (wsnx>3) then
    option%io_buffer = 'Wellsolver needs n>3 direct solver for this system'
    call printErrMsg(option)
  endif

  ! Is water present?

  water_found = PETSC_FALSE
  do iphase = 1, nphase
    if (iphase ==  option%liquid_phase) water_found = PETSC_TRUE
  enddo
  if (.not.water_found) then
    option%io_buffer = 'Well solver assumes water phase is present'
    call printErrMsg(option)
  endif

  ! Is is a supported mode ?

  mode_ok = PETSC_FALSE

  if ((option%iflowmode == TOWG_MODE    ) .or. &
      (option%iflowmode == TOIL_IMS_MODE)) mode_ok = PETSC_TRUE

  if (.not.mode_ok) then
    option%io_buffer = 'This mode not yet supported by well solver'
    call printErrMsg(option)
  endif

  is_grdecl = GetIsGrdecl()
  if (.not.is_grdecl) then
    option%io_buffer = 'Well solver requires grdecl type input'
    call printErrMsg(option)
  endif

end subroutine checkSolverAvailable

! *************************************************************************** !

function getRnormwbs()

  !
  ! Find norm for well bore solution
  !
  ! Author: Dave Ponting
  ! Date  : 04/04/19

  implicit none

  PetscReal :: getRnormwbs,Rnorm
  PetscInt  :: ierr,ixw

  PetscReal :: bl(1), bg(1)

  ierr = 0

  Rnorm = 0.0
  do ixw = 1, ws_nxw
    Rnorm = Rnorm + rwbs(ixw)*rwbs(ixw)
  enddo

  if (w_crossproc) then
    ierr = 0
    bl(1) = Rnorm
    call MPI_Allreduce(bl, bg, ONE_INTEGER, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, w_comm, ierr)
    call checkErr(ierr)
    Rnorm = bg(1)
  endif

  getRnormwbs = sqrt(Rnorm)

end function getRnormwbs

! *************************************************************************** !

function allTrue(vote)

  !
  ! Check that all the procs solving this well agree on a true value
  !
  ! Author: Dave Ponting
  ! Date  : 03/28/19

  implicit none

  PetscBool:: alltrue
  PetscBool, intent(in) :: vote

  PetscInt :: il(1),ig(1),ierr

  il(1) = 1
  ig(1) = 0
  ierr  = 0

  if (w_crossproc) then

    ! Is cross-proc, so only return true iff all votes true

    allTrue = PETSC_FALSE

    ! Increment local to 1 if vote is true

    if (vote) il(1)=1

    !  Get the minimum value (i.e. 0 if any vote is false)

    call MPI_Allreduce(il, ig, ONE_INTEGER, &
                        MPI_INTEGER, MPI_MIN, w_comm, ierr)
    call checkErr(ierr)

    !  Set allTrue to true iff all votes are true

    if (ig(1)==1) allTrue = PETSC_TRUE

else

  ! Not cross-proc, so output follows vote

  allTrue = vote
endif

end function allTrue

! *************************************************************************** !

subroutine checkErr(ierr)

  !
  ! Check for an MPI error
  !
  ! Author: Dave Ponting
  ! Date  : 03/29/19

  implicit none

  PetscInt,intent(in)::ierr

  if (ierr /= 0) ws_MPI_errs = ws_MPI_errs + 1

end subroutine checkErr

end module Well_Solver_module
