module CPR_Preconditioner_module
! Implements a CPR preconditioner using the PCSHELL
! funcitonality of PETSC.
! Daniel Stone and Sebastien Loisel
#include <petsc/finclude/petscsys.h>
#include "petsc/finclude/petscts.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscpc.h"
#include "petsc/finclude/petscviewer.h"
  use petscmat
  use petscksp
  use petscpc
  use petscvec
  use PFLOTRAN_Constants_module
  use Option_module
  implicit none

  private

  type, public :: cpr_pc_type
    Mat :: A, Ap, As
    ! The two stage CPR preconditioner calls two other preconditioners:
    PC :: T1_PC        ! T1 will be a SHELL pc that extracts the pressure system residual from the
                    ! input residual, then approximates the inverse of Ap acting on that residual
                    ! using AMG.
    KSP :: T1_KSP   ! Usually PCONLY, the PC for this KSP is hypre/boomeramg. This is called
                    ! by T1 (above) as the AMG part.
    PC :: T2_PC        ! A regular PC, usually BJACOBI
    PC :: T3_PC
    KSP :: T3_KSP

    Vec :: T1r, T3r, r2, s, z, factors1vec, factors2vec, factors3vec
    PetscInt ::  t2fillin, timescalled, asmoverlap, exrow_offset
    PetscBool :: firstT1Call, firstT2call, firstT3call, asmfactorinplace, &
                 t2shiftinblocks, zeroing, useGAMG, mv_output, t2_zeroing, &
                 skip_T1, amg_report, amg_manual, T1_scale, T3_scale
    character(len=MAXWORDLENGTH) :: T1_type, T2_type, T3_type, extract_type, &
                                    CPR_type
    ! following are buffers/workers for the pressure system extraction.
    ! 1d arrays:
    PetscReal, dimension(:), allocatable :: vals, insert_vals, insert_vals2
    PetscInt, dimension(:), allocatable :: colIdx, colIdx_keep, insert_colIdx
    ! 2d arrays:
    PetscReal, dimension(:,:), allocatable :: all_vals
    ! point at the option object, needed for error outputs
    type(option_type), pointer :: option
  end type cpr_pc_type

  ! interfaces need to make the set and get context routines
  ! work
  interface
    subroutine PCShellSetContext (P_in, ctx_in, ierr_in)
      use petscksp
      Import :: cpr_pc_type
      PC :: P_in
      type(cpr_pc_type) :: ctx_in
      PetscErrorCode :: ierr_in
    end subroutine PCShellSetContext
  end interface

  interface
    subroutine PCShellGetContext (P_in, ctx_in, ierr_in)
      use petscksp
      Import :: cpr_pc_type
      PC :: P_in
      type(cpr_pc_type), pointer :: ctx_in
      PetscErrorCode :: ierr_in
    end subroutine PCShellGetContext
  end interface

public :: DeallocateWorkersInCPRStash, &
          CPRMake

contains

! ************************************************************************** !
!  CPR Apply routines

subroutine CPRApply(p, r, y,ierr)
  ! To be used as a PCSHELL apply routine

  !
  ! Applies the CPR preconditioner to r
  ! (output y)
  !
  ! Author: Sebastien Loisel, Daniel Stone
  ! Modified by: Heeho Park, Jan 2020.
  ! Date: Oct 2017 - March 2018
  !

  implicit none

  ! fixed signature: will approximate
  ! y = inv(A)r using shell preconditioner P:
  PC :: P
  Vec :: r
  Vec :: y
  PetscErrorCode :: ierr

  ! worker vectors:
  Vec :: t1r, t3r, r2
  ! PC context holds workers:
  type(cpr_pc_type), pointer :: ctx
  ! misc variables and parms:
  PetscReal:: one
  Mat :: a

  call PCShellGetContext(p, ctx, ierr); CHKERRQ(ierr)
  one = 1.0
  r2 = ctx%r2
  t1r = ctx%T1r
  t3r = ctx%T3r
  a = ctx%a

  if (ctx%skip_T1) then
    call VecZeroEntries(t1r,ierr); CHKERRQ(ierr)
  else
    call PCApply(ctx%T1_PC,r,t1r,ierr); CHKERRQ(ierr) ! CPR to pres blocks
  endif

  if (ctx%CPR_type == "ADDITIVE") then
    ! solving saturations blocks with amg if diffusion dominated
    call PCApply(ctx%T3_PC,r,t3r,ierr); CHKERRQ(ierr) ! CPR to sat blocks
    call VecAYPX(t1r,one,t3r,ierr); CHKERRQ(ierr) ! t1r = t3r + t1r
  end if

  ! t1r and t3r are pressure and saturation amg solutions, respectively
  ! r2 is a residual after the amg step, r is the original residual.
  call MatMult(a,t1r,r2,ierr); CHKERRQ(ierr)  ! r2 = a*t1r
  call VecAYPX(r2,-one,r,ierr); CHKERRQ(ierr) ! r2 = r - r2

  call PCApply(ctx%T2_PC,r2,y,ierr); CHKERRQ(ierr) ! ILU(0) with Ay=r2
  call VecAYPX(y,one,t1r,ierr); CHKERRQ(ierr) ! y = t1r + y

  !y is the overall changes to pressure and saturation

end subroutine CPRApply

! ************************************************************************** !

subroutine CPRT1Apply(p, x, y,ierr)
  ! To be used as a PCSHELL apply routine

  !
  ! Applies the T1 part of the CPR preconditioner
  ! to r (ourput y)
  !
  ! Author: Sebastien Loisel, Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  implicit none

  ! fixed signature: will approximate
  ! y = inv(A)r using shell preconditioner P:
  PC :: p
  Vec :: x
  Vec :: y
  PetscErrorCode :: ierr

  ! PC context holds workers:
  type(cpr_pc_type), pointer :: ctx
  ! some other KSPs and PCs we'll use:
  KSP :: ksp
  PC :: amg_pc
  ! misc workers, etc:
  PetscInt :: b,start,k,its
  Vec :: s, z

  call PCShellGetContext(p, ctx, ierr); CHKERRQ(ierr)
  s = ctx%s
  z = ctx%z

  call QIRHS(ctx%factors1Vec, ctx%factors2vec, x, s, ierr)

  ksp = ctx%T1_KSP

  if (ctx%T1_type /= "NONE") then
    call KSPSolve(ksp,s,z,ierr); CHKERRQ(ierr)
    call KSPGetIterationNumber(ksp, its, ierr); CHKERRQ(ierr)
  else
    call KSPGetPC(ksp, amg_pc, ierr);CHKERRQ(ierr)
    call PCApply(amg_pc,s,z,ierr); CHKERRQ(ierr)

    if (ctx%amg_report) then
      call PCView(amg_pc, PETSC_VIEWER_STDOUT_WORLD, ierr);CHKERRQ(ierr)
    endif
  endif

  call VecZeroEntries(y,ierr); CHKERRQ(ierr)
  call VecStrideScatter(z,0,y,INSERT_VALUES,ierr); CHKERRQ(ierr)

end subroutine CPRT1Apply

! ************************************************************************** !

subroutine CPRT3Apply(p, x, y,ierr)
  ! To be used as a PCSHELL apply routine

  !
  ! Applies the T3 part of the CPR preconditioner
  ! to r (ourput y)
  ! where saturation block information is running though amg
  !
  ! Author: Heeho Park
  ! January 2020
  !

  implicit none

  ! fixed signature: will approximate
  ! y = inv(A)r using shell preconditioner P:
  PC :: p
  Vec :: x
  Vec :: y
  PetscErrorCode :: ierr

  ! PC context holds workers:
  type(cpr_pc_type), pointer :: ctx
  ! some other KSPs and PCs we'll use:
  KSP :: ksp
  PC :: amg_pc
  ! misc workers, etc:
  PetscInt :: b,start,k,its
  Vec :: s, z

  call PCShellGetContext(p, ctx, ierr); CHKERRQ(ierr)
  s = ctx%s
  z = ctx%z

  call QIRHS(ctx%factors3Vec, ctx%factors2vec, x, s, ierr)

  ksp = ctx%T3_KSP

  if (ctx%T3_type /= "NONE") then
    call KSPSolve(ksp,s,z,ierr); CHKERRQ(ierr)
    call KSPGetIterationNumber(ksp, its, ierr); CHKERRQ(ierr)
  else
    call KSPGetPC(ksp, amg_pc, ierr);CHKERRQ(ierr)
    call PCApply(amg_pc,s,z,ierr); CHKERRQ(ierr)

    if (ctx%amg_report) then
      call PCView(amg_pc, PETSC_VIEWER_STDOUT_WORLD, ierr);CHKERRQ(ierr)
    endif
  endif

  ! supposedly saturation is the 2nd unknown of the block.
  call VecStrideScatter(z,1,y,INSERT_VALUES,ierr); CHKERRQ(ierr)

end subroutine CPRT3Apply

!  end of CPR apply routines
! ************************************************************************** !

! ************************************************************************** !
!  CPR Setup Routines (every time Jaocobian updates)

subroutine CPRSetup(p,ierr)
  ! to be used as a PCSHELL setup routine

  !
  ! sets up the CPR preconditioner, is called
  ! every time the matrix updates
  !
  ! Author:  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  implicit none

  ! fixed signature:
  PC :: p
  PetscErrorCode :: ierr

  ! PC context holds workers:
  type(cpr_pc_type), pointer ::ctx

  call PCShellGetContext(p, ctx, ierr); CHKERRQ(ierr)
  call CPRSetupT1(ctx, ierr)
  if (ctx%CPR_type == "ADDITIVE") then
    call CPRSetupT3(ctx, ierr)
  end if
  call CPRSetupT2(ctx, ierr)

end subroutine CPRSetup

! ************************************************************************** !

subroutine CPRSetupT1(ctx,  ierr)
  !
  ! set up the T1 part of the CPR preconditioner,
  ! mostly by extracting the new pressure system
  !
  ! Author:  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  implicit none

  type(cpr_pc_type) :: ctx
  PetscErrorCode :: ierr

  PetscInt :: b, mx
  Mat :: A
  KSP :: ksp

  A = ctx%A

  call MatGetBlockSize(A,b,ierr); CHKERRQ(ierr)

  if (ctx%firstT1Call) then ! some final initializations that need knowlege
                            ! of the fully assembled Jacobian
    ! let T1 and its inner ksp solver know what their inner operators are
    ksp = ctx%T1_KSP
    call PCSetOperators(ctx%T1_PC,A,A,ierr); CHKERRQ(ierr)
    call KSPSetOperators(ksp,ctx%Ap,ctx%Ap,ierr); CHKERRQ(ierr)

    ! now sparsity pattern of the Jacobian is defined,
    ! allocate worker arrays that are big enough to
    ! be used in the extraction routines coming up.
    call MatGetMaxRowCount(A, mx, ierr)
    call AllocateWorkersInCPRStash(ctx, mx, b)

    ctx%firstT1Call = .FALSE.
  end if

  call MatZeroEntries(ctx%Ap, ierr); CHKERRQ(ierr)
  call MatZeroEntries(ctx%As, ierr); CHKERRQ(ierr)

  select case(ctx%extract_type)
    case('ABF')
      if (b == 2) then
        call MatGetSubABFImmiscible(A, ctx%Ap, ctx%As, ctx%factors1vec, &
                                    ctx%factors3vec, ierr, ctx)
      else
        ctx%option%io_buffer = 'ABF not available for more than 3 unknowns per cell'
        call PrintErrMsg(ctx%option)
      end if
    case('QIMPES_TWO_UNKNOWNS')
      call MatGetSubQIMPESImmiscible(A, ctx%Ap, ctx%As, ctx%factors1vec, &
                                     ctx%factors3vec, ierr, ctx)
    case('QIMPES_THREE_UNKNOWNS')
      ! we have a more efficient version for 3x3 blocks so do that if we can instead
      if (b == 3) then
        call MatGetSubQIMPES(A, ctx%Ap, ctx%factors1vec,  ierr, ctx)
      else  ! talk to Daniel or Paolo to remove this statement - Heeho
        call MatGetSubQIMPES_var(A, ctx%Ap, ctx%factors1vec,  ierr,   b,  ctx)
      end if
    case('QIMPES')
      if (b == 2) then
        ! more efficient for 2x2 blocks
        call MatGetSubQIMPESImmiscible(A, ctx%Ap, ctx%As, ctx%factors1vec, &
                                       ctx%factors3vec, ierr, ctx)
      else if (b == 3) then
        call MatGetSubQIMPES(A, ctx%Ap, ctx%factors1vec,  ierr, ctx)
      else
        call MatGetSubQIMPES_var(A, ctx%Ap, ctx%factors1vec,  ierr,   b,  ctx)
      end if
    case('QIMPES_VARIABLE_FORCE')
      ! force to use the variables block size implementation even for 3x3 blocks
      call MatGetSubQIMPES_var(A, ctx%Ap, ctx%factors1vec,  ierr,   b,  ctx)
    case default
      ctx%option%io_buffer = 'CPRSetupT1, extraction type not defined'
      call PrintErrMsg(ctx%option)
  end select

end subroutine CPRSetupT1

! ************************************************************************** !

subroutine CPRSetupT3(ctx,  ierr)

  ! Decoupling of saturation block is done while the pressure block
  ! was decoupled in CPRSetupT1.
  !
  ! Only the simple declaration is required here.
  !
  ! Author:  Heeho Park
  ! January 2020


  implicit none

  type(cpr_pc_type) :: ctx
  PetscErrorCode :: ierr

  KSP :: ksp

  if (ctx%firstT3Call) then
    ksp = ctx%T3_KSP
    call PCSetOperators(ctx%T3_PC,ctx%A,ctx%A,ierr); CHKERRQ(ierr)
    call KSPSetOperators(ksp,ctx%As,ctx%As,ierr); CHKERRQ(ierr)
    ctx%firstT3Call = .FALSE.
  end if

end subroutine CPRSetupT3

! ************************************************************************** !

subroutine CPRSetupT2(ctx, ierr)
  !
  ! set up the T2 part of the CPR preconditioner,
  ! if anything needs to be done
  !
  ! Author:  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  use String_module

  implicit none

  type(cpr_pc_type) :: ctx
  PetscErrorCode :: ierr

  PC :: T2, pc_inner
  Mat :: A
  PetscInt :: nsub_ksp, first_sub_ksp, i
  KSP, pointer :: sub_ksps(:)

  ! set operator on first call.

  ! also if first call then time is right
  ! to adjust things like fill in, etc.

  ! if BJACOBI, this is the best place to try
  ! modifying sub-ksps
  if (ctx%firstT2Call) then

  T2 = ctx%T2_PC

    if (StringCompare(ctx%T2_type, 'PCASM')) then
      call PCASMSetOverlap(T2, ctx%asmoverlap, ierr);CHKERRQ(ierr)
    endif

    A = ctx%A
    call PCSetOperators(T2,A,A,ierr); CHKERRQ(ierr)
    call PCSetFromOptions(T2, ierr);CHKERRQ(ierr)
    call PCSetUp(T2, ierr); CHKERRQ(ierr)

    if (StringCompare(ctx%T2_type, 'PCASM')) then

      ! default should be preonly
      call PCASMGetSubKSP(T2,   nsub_ksp, first_sub_ksp, &
                               PETSC_NULL_KSP, ierr); CHKERRQ(ierr)
      !                         ksp array
      ! allocate ksp array now number known:
      allocate(sub_ksps(nsub_ksp))
      ! call again:
      call PCASMGetSubKSP(T2,   nsub_ksp, first_sub_ksp, &
                               sub_ksps, ierr); CHKERRQ(ierr)

      do i = 1,nsub_ksp
        call KSPGetPC(sub_ksps(i), pc_inner,ierr);CHKERRQ(ierr)
        if (ctx%t2shiftinblocks) then
          call PCFactorSetShiftType(pc_inner,MAT_SHIFT_INBLOCKS,&
                                    ierr);CHKERRQ(ierr)
        endif
        call PCFactorSetLevels(pc_inner, ctx%t2fillin, ierr); CHKERRQ(ierr)
        if (ctx%asmfactorinplace) then
          call PCFactorSetUseInPlace(pc_inner, PETSC_TRUE,  ierr);CHKERRQ(ierr)
        endif
      enddo
      deallocate(sub_ksps)
      nullify(sub_ksps)
    ! specifically for block jacobi:
    elseif (StringCompare(ctx%T2_type, 'PCBJACOBI')) then

      ! default should be preonly
      call PCBJacobiGetSubKSP(T2,   nsub_ksp, first_sub_ksp, &
                               PETSC_NULL_KSP, ierr); CHKERRQ(ierr)
      !                         ksp array
      ! allocate ksp array now number known:
      allocate(sub_ksps(nsub_ksp))
      ! call again:
      call PCBJacobiGetSubKSP(T2,   nsub_ksp, first_sub_ksp, &
                               sub_ksps, ierr); CHKERRQ(ierr)

      do i = 1,nsub_ksp
        ! default should be ilu with 0 fill
        call KSPGetPC(sub_ksps(i), pc_inner, ierr); CHKERRQ(ierr)

        ! MUST DO THIS: the pc is otherwise not setup by this
        ! point and will not accept changes like factor levels.
        call PCSetType(pc_inner, PCILU, ierr); CHKERRQ(ierr)
        call PCSetUp(PC_inner, ierr); CHKERRQ(ierr)

        call PCFactorSetLevels(pc_inner, ctx%t2fillin, ierr); CHKERRQ(ierr)

      enddo
      deallocate(sub_ksps)
      nullify(sub_ksps)
    endif

    ctx%firstT2call = PETSC_FALSE
  endif

end subroutine CPRSetupT2

!  end of CPR Setup Routines
! ************************************************************************** !

! ************************************************************************** !
!  CPR Creation Routines

subroutine CPRMake(p, ctx, c, ierr, option)
  !
  ! make the CPR preconditioner
  ! create all necessary sub PC/KSP objects
  ! and set up as much as can be done at this
  ! point
  !
  ! Author: Sebastien Loisel,  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !
  use String_module

  implicit none

  PC :: p
  PetscErrorCode :: ierr
  type(cpr_pc_type) :: ctx
  MPI_Comm :: c
  type(option_type), target :: option

  ctx%option => option

  call PCSetType(p,PCSHELL,ierr); CHKERRQ(ierr)
  call PCShellSetApply(p,CPRapply,ierr); CHKERRQ(ierr)
  call PCShellSetSetUp(p,CPRSetup,ierr); CHKERRQ(ierr)

  call PCShellSetContext(p, ctx, ierr); CHKERRQ(ierr)

#if PETSC_VERSION_GE(3,11,99)
#if !defined(PETSC_HAVE_HYPRE)
  if (.not. ctx%useGAMG .or. &
      (StringCompare(ctx%T2_type,'SAILS') .or. &
       StringCompare(ctx%T2_type,'PILUT') .or. &
       StringCompare(ctx%T2_type,'EUCLID'))) then
    option%io_buffer = 'CPR solver settings require that PETSc be &
      &configured with hypre (--download-hypre=yes).'
    call PrintErrMsg(option)
  endif
#endif
#else
#if !defined(PETSC_HAVE_LIBHYPRE)
  if (.not. ctx%useGAMG .or. &
      (StringCompare(ctx%T2_type,'SAILS') .or. &
       StringCompare(ctx%T2_type,'PILUT') .or. &
       StringCompare(ctx%T2_type,'EUCLID'))) then
    option%io_buffer = 'CPR solver settings require that PETSc be &
      &configured with hypre (--download-hypre=yes).'
    call PrintErrMsg(option)
  endif
#endif
#endif


  call CPRCreateT1(c, ctx,  ierr); CHKERRQ(ierr)
  if (ctx%CPR_type == "ADDITIVE") then
    call CPRCreateT3(c, ctx,  ierr); CHKERRQ(ierr)
  end if
  call CPRCreateT2(c, ctx,  ierr); CHKERRQ(ierr)

end subroutine CPRMake

! ************************************************************************** !

subroutine CPRCreateT1(c,  ctx,   ierr)
  !
  ! This will perform bare minimum
  ! creation of objects for T1 that do not
  ! need an existing system (i.e. matrix A).
  ! The system dependent setup must be repeated
  ! every step and therefore is seperated out.
  !
  ! Author: Sebastien Loisel,  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  implicit none

  MPI_Comm :: c
  type(cpr_pc_type) :: ctx
  PetscErrorCode :: ierr

  KSP :: ksp
  PC :: prec

  ! nice default options for boomeramg:
  if (.NOT. ctx%amg_manual) then
    call SetCPRDefaults(ierr)
  endif

  call KSPCreate(c,ksp,ierr); CHKERRQ(ierr)

  select case(ctx%T1_type)
    case('RICHARDSON')
      call KSPSetType(ksp,KSPRICHARDSON,ierr); CHKERRQ(ierr)
    case('FGMRES')
      call KSPSetType(ksp,KSPFGMRES,ierr); CHKERRQ(ierr)
      call KSPSetTolerances(ksp, 1.0d-3, 1.d0-3, PETSC_DEFAULT_REAL, &
                            PETSC_DEFAULT_INTEGER, ierr);CHKERRQ(ierr)
    case('GMRES')
      call KSPSetType(ksp,KSPGMRES,ierr); CHKERRQ(ierr)
    case default
      ! really we will just skip over the ksp entirely
      ! in this case, but for completeness..
      call KSPSetType(ksp,KSPPREONLY,ierr); CHKERRQ(ierr)
  end select

  call KSPGetPC(ksp,prec,ierr); CHKERRQ(ierr)

  if (ctx%useGAMG) then
    call PCSetType(prec,PCGAMG,ierr); CHKERRQ(ierr)
  else
    call PCSetType(prec,PCHYPRE,ierr); CHKERRQ(ierr)
#if PETSC_VERSION_GE(3,11,99)
#if defined(PETSC_HAVE_HYPRE)
    call PCHYPRESetType(prec,"boomeramg",ierr); CHKERRQ(ierr)
#endif
#else
#if defined(PETSC_HAVE_LIBHYPRE)
    call PCHYPRESetType(prec,"boomeramg",ierr); CHKERRQ(ierr)
#endif
#endif

  endif

  call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
  call PetscObjectSetName(ksp,"T1",ierr); CHKERRQ(ierr)

  ctx%T1_KSP = ksp
  call PCCreate(C,ctx%T1_PC,ierr); CHKERRQ(ierr)
  call PCSetType(ctx%T1_PC,PCSHELL,ierr); CHKERRQ(ierr)
  call PCShellSetApply(ctx%T1_PC,CPRT1apply,ierr); CHKERRQ(ierr)

  call PCShellSetContext(ctx%T1_PC, ctx, ierr); CHKERRQ(ierr)

end subroutine CPRCreateT1

! ************************************************************************** !

subroutine CPRCreateT3(c,  ctx,   ierr)
  !
  ! The parameters and options of the AMG solver are consistent with
  ! T1 and cannot create separtion options for this; for now.
  ! if some solver, does well for parabolic equations then we'll separate
  ! this from T1.
  !
  ! Author: Heeho Park
  ! January 2020
  !

  implicit none

  MPI_Comm :: c
  type(cpr_pc_type) :: ctx
  PetscErrorCode :: ierr

  KSP :: ksp
  PC :: prec

  call KSPCreate(c,ksp,ierr); CHKERRQ(ierr)

  select case(ctx%T3_type)
    case('RICHARDSON')
      call KSPSetType(ksp,KSPRICHARDSON,ierr); CHKERRQ(ierr)
    case('FGMRES')
      call KSPSetType(ksp,KSPFGMRES,ierr); CHKERRQ(ierr)
      call KSPSetTolerances(ksp, 1.0d-3, 1.d0-3, PETSC_DEFAULT_REAL, &
                            PETSC_DEFAULT_INTEGER, ierr);CHKERRQ(ierr)
    case('GMRES')
      call KSPSetType(ksp,KSPGMRES,ierr); CHKERRQ(ierr)
    case default
      ! really we will just skip over the ksp entirely
      ! in this case, but for completeness..
      call KSPSetType(ksp,KSPPREONLY,ierr); CHKERRQ(ierr)
  end select

  call KSPGetPC(ksp,prec,ierr); CHKERRQ(ierr)

  if (ctx%useGAMG) then
    call PCSetType(prec,PCGAMG,ierr); CHKERRQ(ierr)
  else
    call PCSetType(prec,PCHYPRE,ierr); CHKERRQ(ierr)
#if PETSC_VERSION_GE(3,11,99)
#if defined(PETSC_HAVE_HYPRE)
    call PCHYPRESetType(prec,"boomeramg",ierr); CHKERRQ(ierr)
#endif
#else
#if defined(PETSC_HAVE_LIBHYPRE)
    call PCHYPRESetType(prec,"boomeramg",ierr); CHKERRQ(ierr)
#endif
#endif

  endif

  call KSPSetFromOptions(ksp,ierr); CHKERRQ(ierr)
  call PetscObjectSetName(ksp,"T3",ierr); CHKERRQ(ierr)

  ctx%T3_KSP = ksp
  call PCCreate(C,ctx%T3_PC,ierr); CHKERRQ(ierr)
  call PCSetType(ctx%T3_PC,PCSHELL,ierr); CHKERRQ(ierr)
  call PCShellSetApply(ctx%T3_PC,CPRT3apply,ierr); CHKERRQ(ierr)

  call PCShellSetContext(ctx%T3_PC, ctx, ierr); CHKERRQ(ierr)

end subroutine CPRCreateT3

! ************************************************************************** !

subroutine CPRCreateT2(c, ctx, ierr)
  !
  ! create the T2 preconditioner for the CPR PC
  !
  ! Author: Sebastien Loisel,  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  implicit none

  MPI_Comm :: c
  type(cpr_pc_type) :: ctx
  PetscErrorCode :: ierr

  PC :: t2

  call PCCreate(c,t2,ierr); CHKERRQ(ierr)

  !select case(ctx%T2_type)
  select case(trim(ctx%T2_type))
    case('SAILS')
      call PCSetType(T2,PCHYPRE,ierr); CHKERRQ(ierr)
#if PETSC_VERSION_GE(3,11,99)
#if defined(PETSC_HAVE_HYPRE)
      call PCHYPRESetType(t2,"parasails",ierr); CHKERRQ(ierr)
#endif
#else
#if defined(PETSC_HAVE_LIBHYPRE)
      call PCHYPRESetType(t2,"parasails",ierr); CHKERRQ(ierr)
#endif
#endif
    case('PBJ')
      call PCSetType(t2,PCPBJACOBI,ierr); CHKERRQ(ierr)
    case('NONE')
      call PCSetType(t2,PCNONE,ierr); CHKERRQ(ierr)
    case('PCASM')
      call PCSetType(t2,PCASM,ierr); CHKERRQ(ierr)
    case('PCGASM')
      call PCSetType(t2,PCGASM,ierr); CHKERRQ(ierr)
    case('PILUT')
      call PCSetType(t2,PCHYPRE,ierr); CHKERRQ(ierr)
#if PETSC_VERSION_GE(3,11,99)
#if defined(PETSC_HAVE_HYPRE)
      call PCHYPRESetType(t2,"pilut",ierr); CHKERRQ(ierr)
#endif
#else
#if defined(PETSC_HAVE_LIBHYPRE)
      call PCHYPRESetType(t2,"pilut",ierr); CHKERRQ(ierr)
#endif
#endif
    case('EUCLID')
      call PCSetType(t2,PCHYPRE,ierr); CHKERRQ(ierr)
#if PETSC_VERSION_GE(3,11,99)
#if defined(PETSC_HAVE_HYPRE)
      call PCHYPRESetType(t2,"euclid",ierr); CHKERRQ(ierr)
#endif
#else
#if defined(PETSC_HAVE_LIBHYPRE)
      call PCHYPRESetType(t2,"euclid",ierr); CHKERRQ(ierr)
#endif
#endif
    case('ILU')
      call PCSetType(t2,PCILU,ierr); CHKERRQ(ierr)
    case default
      call PCSetType(t2,PCBJACOBI,ierr); CHKERRQ(ierr)
  end select

  ctx%T2_PC = t2

end subroutine CPRCreateT2

!  end of CPR Creation Routines
! ************************************************************************** !

! ************************************************************************** !
! suplementary setup/init/deinit/routines

subroutine SetCPRDefaults(ierr)
  !
  ! set in Petsc's options tables some good
  ! default options for th HYPRE BOOMERAMG
  ! AMG preconditioner.
  !
  ! Author:  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  implicit none

  PetscErrorCode :: ierr

  character(len=MAXSTRINGLENGTH) :: string, word

  ! set sensible defualts for the boomeramg pc here,
  ! the defaults it comes with are rarely preferable to
  ! us.

  ! strong threshold, 0.5 is ok for 3D
  string =  '-pc_hypre_boomeramg_strong_threshold'
  word   =  '0.5'
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            trim(string),trim(word), &
                            ierr);CHKERRQ(ierr)

  ! coarsen type, PMIS is generally more efficient
  string =  '-pc_hypre_boomeramg_coarsen_type'
  word   =  'PMIS'
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            trim(string),trim(word), &
                            ierr);CHKERRQ(ierr)

  ! interpolation type, ext+i is reccomended
  string =  '-pc_hypre_boomeramg_interp_type'
  word   =  'ext+i'
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            trim(string),trim(word), &
                            ierr);CHKERRQ(ierr)

  ! relaxer, Jacobi should be super weak but cheap
  string =  '-pc_hypre_boomeramg_relax_type_all'
  word   =  'Jacobi'
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            trim(string),trim(word), &
                            ierr);CHKERRQ(ierr)

end subroutine SetCPRDefaults

! ************************************************************************** !

subroutine AllocateWorkersInCPRStash(ctx, n, b)
  !
  ! allocate the various arrays that are used to hold data
  ! during the pressure extraction phase of the CPR pc.
  !
  ! Author:  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  implicit none


  type(cpr_pc_type) :: ctx
  PetscInt :: entriesInARow, entriesInAReducedRow

  PetscInt :: n, b

  ! n: the maxiumum number of nonzero columns in any row
  !    of the matrix we will work with
  ! b: the block size of the matrix we will work with

  ! entriesInARow: the maximum number of entries in a row in
  !               the matrix we will work with, plus a
  !               safety buffer in case it was wrong
  entriesInARow = n + 10

  ! entriesInAReducedRow: the maximum number of entries in a row
  !                       of the extracted matrix we will work with.
  !
  entriesInAReducedRow = n/b
  !                       plus a safety buffer in case it was wrong
  entriesInAReducedRow = entriesInAReducedRow + 10

  ! vals: a real buffer to hold the matrix values output from
  !       matgetrows
  allocate(ctx%vals (0:entriesInARow))
  ! insert_vals: a real buffer to hold the values we will
  !               insert to the reduced matrix
  allocate(ctx%insert_vals (0:entriesInAReducedRow))
  allocate(ctx%insert_vals2 (0:entriesInAReducedRow))
  ! colIdx: an integer buffer to hold the column indexes
  !         output from matgetvals
  allocate(ctx%colIdx (0:entriesInARow))
  ! colIdx_keep: will copy colIdx into here, since matrestorerows
  !              resets colIdx
  allocate(ctx%colIdx_keep (0:entriesInARow))
  ! insert_colIdx: an integer buffer of the column indexes
  !                to be input to the reduced matrix we will
  !                work with
  allocate(ctx%insert_colIdx (0:entriesInAReducedRow))
  ! all_vals: will just copy every value from the current row block
  !           (i.e. set of b rows) into here to work on more
  !           easilly. Can very much be improved.
  allocate(ctx%all_vals (0:b-1, 0:entriesInARow))

end subroutine AllocateWorkersInCPRStash

! ************************************************************************** !

subroutine DeallocateWorkersInCPRStash(ctx)
  !
  ! DEallocate the various arrays that are used to hold data
  ! during the pressure extraction phase of the CPR pc.
  !
  ! Author:  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  implicit none


  type(cpr_pc_type) :: ctx

  if (allocated(ctx%vals))deallocate(ctx%vals)
  if (allocated(ctx%insert_vals))deallocate(ctx%insert_vals)
  if (allocated(ctx%insert_vals2))deallocate(ctx%insert_vals2)
  if (allocated(ctx%colIdx))deallocate(ctx%colIdx)
  if (allocated(ctx%colIdx_keep))deallocate(ctx%colIdx_keep)
  if (allocated(ctx%insert_colIdx))deallocate(ctx%insert_colIdx)
  if (allocated(ctx%all_vals))deallocate(ctx%all_vals)

  !! also the pointers:
  nullify(ctx%option)

end subroutine DeallocateWorkersInCPRStash

! end of suplementary setup/init/deinit/routines
! ************************************************************************** !


! ************************************************************************** !
subroutine MatGetSubABFImmiscible(A, App, Ass, factors1Vec,  &
                                  factors3Vec, ierr, ctx)
  !
  ! extraction of the pressure system matrix of for the
  ! CPR preconditioner with
  ! Alternate Block Factorization (ABF)
  ! Decoupling method for 2 unknowns (pressure and saturation)
  ! in the case of ADDITIVE option,
  ! this subroutine extracts both pressure and saturation system of the matrix

  ! Author:  Heeho Park
  ! Date: January 2020.

  ! applies to 2x2 blocks ONLY (composed of Jacobian element, j)
  ! [j_pp j_ps]
  ! [j_sp j_ss]


  implicit none


  Mat :: A, App, Ass
  Vec :: factors1Vec, factors3Vec
  PetscErrorCode :: ierr
  type(cpr_pc_type) :: ctx

  PetscInt, dimension(0:0) :: insert_rows
  ! misc workers:
  PetscReal :: j_pp, j_ps, j_sp, j_ss, lambda_inv, d_ss_abf, d_ps_abf, &
               fac0, fac1, fac3, fac4, scaling_factor, scaling_factor2
  PetscInt :: block_size, rows, cols, num_blocks, num_blocks_local, &
              first_row, cur_col_index, num_col_blocks, &
              diag_row_index, loop_index, i, j, num_cols
  PetscMPIInt :: rank, row_start, row_end


  ctx%vals = 0.d0
  ctx%insert_vals = 0.d0
  ctx%insert_vals2 = 0.d0
  ctx%all_vals = 0.d0

  ctx%colIdx = 0
  ctx%colIdx_keep = 0
  ctx%insert_colIdx = 0
  num_cols = 0

  call MPI_Comm_Rank(PETSC_COMM_WORLD, rank, ierr)
  call MatGetOwnershipRange(A, row_start, row_end, ierr); CHKERRQ(ierr)
  call MatGetBlockSize(A,block_size,ierr); CHKERRQ(ierr)
  call MatGetSize(A, rows, cols, ierr); CHKERRQ(ierr)
  if (rows /= cols) then
     ctx%option%io_buffer = 'MatGetSubQIMPES, given a nonsquare matrix'
    call PrintErrMsg(ctx%option)
  endif
  num_blocks = rows/block_size
  num_blocks_local = (row_end-row_start)/block_size


  ! loop over the row blocks
  do i = 0,num_blocks_local-1

    first_row = i*block_size + row_start

    ! a) extract [j_pp j_ps] of the diagonal block
    !    and all the values of the first row
    call MatGetRow(A, first_row, num_cols, ctx%colIdx, ctx%vals, ierr); CHKERRQ(ierr)
    do j = 0,num_cols-1
      ! a.1) store all the values in the first row and their indices
      ctx%all_vals(0, j) = ctx%vals(j)
      ctx%colIdx_keep(j) = ctx%colIdx(j)
    end do

    diag_row_index = -1
    do loop_index = 0,num_cols-1,block_size
      if (ctx%colIdx(loop_index) == first_row) then
        diag_row_index = loop_index
      endif
    enddo

    if (diag_row_index == -1) then
      ctx%option%io_buffer = 'MatGetSubQIMPESImmiscible, cannot find &
                              &diagonal entry, check matrix.'
      call PrintErrMsg(ctx%option)
    endif

    ! a.2) extract j_pp j_ps
    j_pp = ctx%vals(diag_row_index)
    j_ps = ctx%vals(diag_row_index+1)

    call MatRestoreRow(A, first_row, num_cols, ctx%colIdx, ctx%vals, ierr)
    CHKERRQ(ierr)

    ! b) extract second row
    call MatGetRow(A, first_row+1, num_cols, PETSC_NULL_INTEGER, ctx%vals, &
                  ierr); CHKERRQ(ierr)
    do j = 0,num_cols-1
      ctx%all_vals(1, j) = ctx%vals(j)
    end do

    ! b.1) extract [j_sp j_ss] and invert the value
    j_sp = ctx%vals(diag_row_index)
    j_ss = ctx%vals(diag_row_index+1)
    lambda_inv = 1.0d0/(j_pp*j_ss-j_ps*j_sp)

    ! c) storing factors to later multiply to the RHS vector, b, in QIRHS
    ! inv(Lamda)*D_ss*r_p - inv(Labmda)*D_ps*r_s -> fac0*r_p + fac1*r_s
    ! scaling by the first column shown by Daniel's implementation
    ! this comes from numerical experiments done by David Ponting.
    ! The scaling factor keeps the shape of long waves of diffusion
    ! which AMG is really good at resolving; hence, gaining large speed-up.
    if (ctx%T1_scale) then
      scaling_factor = abs(j_pp) + abs(j_sp)
    else
      scaling_factor = 1.d0
    end if
    fac0 = lambda_inv*j_ss*scaling_factor
    fac1 = -lambda_inv*j_ps*scaling_factor
    call VecSetValue(factors1Vec, first_row, fac0, INSERT_VALUES, ierr)
    CHKERRQ(ierr)
    call VecSetValue(factors1Vec, first_row+1, fac1, INSERT_VALUES, ierr)
    CHKERRQ(ierr)

    if (ctx%CPR_type == "ADDITIVE") then
      ! this is effective when you expect saturation to be
      ! diffusion-dominant
      if (ctx%T3_scale) then
        scaling_factor2 = abs(j_ps) + abs(j_ss)
      else
        scaling_factor2 = 1.d0
      end if
      fac3 = -lambda_inv*j_sp*scaling_factor2
      fac4 = lambda_inv*j_pp*scaling_factor2
      call VecSetValue(factors3Vec, first_row, fac3, INSERT_VALUES, ierr)
      CHKERRQ(ierr)
      call VecSetValue(factors3Vec, first_row+1, fac4, INSERT_VALUES, ierr)
      CHKERRQ(ierr)
    end if

    num_col_blocks = num_cols/block_size
    call MatRestoreRow(A, first_row+1, num_cols, PETSC_NULL_INTEGER, &
                       ctx%vals, ierr); CHKERRQ(ierr)

    ! d) prepare to set values
    insert_rows = i + row_start/block_size
    do j = 0,num_col_blocks-1
      cur_col_index = j*block_size
      ctx%insert_colIdx(j) = ctx%colIdx_keep(cur_col_index)/block_size

      ! lam_inv*Dss*App-lam_inv*Dps*Asp
      ! apply decoupling to the left hand side
      ctx%insert_vals(j) = fac0*ctx%all_vals(0,cur_col_index) &
                           + fac1*ctx%all_vals(1,cur_col_index)
      if (ctx%CPR_type == "ADDITIVE") then
        ! -lam_inv*Dsp*Aps+lam_inv*Dpp*Ass
        ctx%insert_vals2(j) = fac3*ctx%all_vals(0,cur_col_index+1) &
                              + fac4*ctx%all_vals(1,cur_col_index+1)
      end if
    end do

    ! e) set values
    call MatSetValues(App, 1, insert_rows, num_col_blocks, &
                      ctx%insert_colIdx(0:num_col_blocks-1), &
                      ctx%insert_vals(0:num_col_blocks-1), INSERT_VALUES, ierr)
    CHKERRQ(ierr)

    if (ctx%CPR_type == "ADDITIVE") then
      call MatSetValues(Ass, 1, insert_rows, num_col_blocks, &
                        ctx%insert_colIdx(0:num_col_blocks-1), &
                        ctx%insert_vals2(0:num_col_blocks-1), INSERT_VALUES, ierr)
      CHKERRQ(ierr)
    end if

  end do

  call MatAssemblyBegin(App,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
  call MatAssemblyEnd(App,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)

  call VecAssemblyBegin(factors1Vec, ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(factors1vec, ierr);CHKERRQ(ierr)

  if (ctx%CPR_type == "ADDITIVE") then
    call MatAssemblyBegin(Ass,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(Ass,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    call VecAssemblyBegin(factors3Vec, ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(factors3vec, ierr);CHKERRQ(ierr)
  end if

end subroutine MatGetSubABFImmiscible


! ************************************************************************** !

subroutine MatGetSubQIMPESImmiscible(A, App, Ass, factors1Vec, &
                                     factors3Vec, ierr, ctx)
  !
  ! extraction of the pressure system matrix of for the
  ! CPR preconditioner with
  ! Quasi Implicit Pressure Explicit Saturation (IMPES)
  ! Decoupling method for 2 unknowns (pressure and saturation)

  ! Author:  Heeho Park
  ! Date: January 2020.

  ! applies to 2x2 blocks ONLY (composed of Jacobian element, j)
  ! [j_pp j_ps]
  ! [j_sp j_ss]


  implicit none


  Mat :: A, App, Ass
  Vec :: factors1Vec, factors3Vec
  PetscErrorCode :: ierr
  type(cpr_pc_type) :: ctx

  PetscInt, dimension(0:0) :: insert_rows
  ! misc workers:
  PetscReal :: j_ps, j_ss, j_ps_j_ss_inv, fac0, fac1, &
               j_pp, j_sp, fac3, fac4, j_sp_j_pp_inv
  PetscInt :: block_size, rows, cols, num_blocks, num_blocks_local, &
              first_row, cur_col_index, num_col_blocks, &
              diag_row_index, loop_index, i, j, num_cols
  PetscMPIInt :: rank, row_start, row_end


  ctx%vals = 0.d0
  ctx%insert_vals = 0.d0
  ctx%all_vals = 0.d0

  ctx%colIdx = 0
  ctx%colIdx_keep = 0
  ctx%insert_colIdx = 0
  num_cols = 0

  call MPI_Comm_Rank(PETSC_COMM_WORLD, rank, ierr)
  call MatGetOwnershipRange(A, row_start, row_end, ierr); CHKERRQ(ierr)
  call MatGetBlockSize(A,block_size,ierr); CHKERRQ(ierr)
  call MatGetSize(A, rows, cols, ierr); CHKERRQ(ierr)
  if (rows /= cols) then
     ctx%option%io_buffer = 'MatGetSubQIMPES, given a nonsquare matrix'
    call PrintErrMsg(ctx%option)
  endif
  num_blocks = rows/block_size
  num_blocks_local = (row_end-row_start)/block_size


  ! loop over the row blocks
  do i = 0,num_blocks_local-1

    first_row = i*block_size + row_start

    ! a) extract j_ps of the diagonal block and all the values of the first row
    call MatGetRow(A, first_row, num_cols, ctx%colIdx, ctx%vals, ierr); CHKERRQ(ierr)
    do j = 0,num_cols-1
      ! a.1) store all the values in the first row and their indices
      ctx%all_vals(0, j) = ctx%vals(j)
      ctx%colIdx_keep(j) = ctx%colIdx(j)
    end do

    diag_row_index = -1
    do loop_index = 0,num_cols-1,block_size
      if (ctx%colIdx(loop_index) == first_row) then
        diag_row_index = loop_index
      endif
    enddo

    if (diag_row_index == -1) then
      ctx%option%io_buffer = 'MatGetSubQIMPESImmiscible, cannot find &
                              &diagonal entry, check matrix.'
      call PrintErrMsg(ctx%option)
    endif

    ! a.2) extract j_pp,j_ps
    j_pp = ctx%vals(diag_row_index)
    j_ps = ctx%vals(diag_row_index+1)
    call MatRestoreRow(A, first_row, num_cols, ctx%colIdx, ctx%vals, ierr)
    CHKERRQ(ierr)

    ! b) extract second row
    call MatGetRow(A, first_row+1, num_cols, PETSC_NULL_INTEGER, ctx%vals, &
                  ierr); CHKERRQ(ierr)
    do j = 0,num_cols-1
      ctx%all_vals(1, j) = ctx%vals(j)
    end do

    ! b.1) extract j_sp, j_ss and invert the value
    j_sp = ctx%vals(diag_row_index)
    j_ss = ctx%vals(diag_row_index+1)
    if (j_ss == 0.d0) then
      ! to avoid NaNs
      j_ps_j_ss_inv = 0.d0
    else
      j_ps_j_ss_inv = j_ps*1.d0/j_ss
    end if

    num_col_blocks = num_cols/block_size
    call MatRestoreRow(a, first_row+1, num_cols, PETSC_NULL_INTEGER, &
                       ctx%vals, ierr); CHKERRQ(ierr)

    ! c) storing factors to later multiply to the RHS vector, b, in QIRHS
    ! r_p - D_ps*inv(D_ss)*r_s -> fac0*r_p + fac1*r_s
    fac0 = 1.d0
    fac1 = -j_ps_j_ss_inv
    call VecSetValue(factors1Vec, first_row, fac0, INSERT_VALUES, ierr)
    CHKERRQ(ierr)
    call VecSetValue(factors1Vec, first_row+1, fac1, INSERT_VALUES, ierr)
    CHKERRQ(ierr)
    if (ctx%CPR_type == "ADDITIVE") then
      ! -D_sp*inv(D_pp)*r_p + r_s -> fac3*r_p + fac4*r_s
      ! this is effective when you expect saturation to be
      ! diffusion-dominant
      if (j_pp == 0.d0) then
        ! to avoid NaNs
        j_sp_j_pp_inv = 0.d0
      else
        j_sp_j_pp_inv = j_sp*1.d0/j_pp
      end if
      fac3 = -j_sp_j_pp_inv
      fac4 = 1.d0
      call VecSetValue(factors3Vec, first_row, fac3, INSERT_VALUES, ierr)
      CHKERRQ(ierr)
      call VecSetValue(factors3Vec, first_row+1, fac4, INSERT_VALUES, ierr)
      CHKERRQ(ierr)
    end if

    ! d) prepare to set values
    insert_rows = i + row_start/block_size
    do j = 0,num_col_blocks-1
      cur_col_index = j*block_size
      ctx%insert_colIdx(j) = ctx%colIdx_keep(cur_col_index)/block_size

      ! Apply decoupling to the left hand side matrix
      ! A_pp - D_ps*inv(D_ss)*A_sp
      ctx%insert_vals(j) = fac0*ctx%all_vals(0,cur_col_index) &
                           + fac1*ctx%all_vals(1,cur_col_index)
      if (ctx%CPR_type == "ADDITIVE") then
        ! -D_sp*inv(D_pp)*A_ps + A_ss
        ctx%insert_vals2(j) = fac3*ctx%all_vals(0,cur_col_index+1) &
                              + fac4*ctx%all_vals(1,cur_col_index+1)
      end if
    end do

    ! e) set values
    call MatSetValues(App, 1, insert_rows, num_col_blocks, &
                      ctx%insert_colIdx(0:num_col_blocks-1), &
                      ctx%insert_vals(0:num_col_blocks-1), INSERT_VALUES, ierr)
    CHKERRQ(ierr)
    if (ctx%CPR_type == "ADDITIVE") then
      call MatSetValues(Ass, 1, insert_rows, num_col_blocks, &
                        ctx%insert_colIdx(0:num_col_blocks-1), &
                        ctx%insert_vals2(0:num_col_blocks-1), INSERT_VALUES, ierr)
      CHKERRQ(ierr)
    end if
  end do

  call MatAssemblyBegin(App,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
  call MatAssemblyEnd(App,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
  call VecAssemblyBegin(factors1Vec, ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(factors1vec, ierr);CHKERRQ(ierr)

  if (ctx%CPR_type == "ADDITIVE") then
    call MatAssemblyBegin(Ass,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(Ass,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
    call VecAssemblyBegin(factors3Vec, ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(factors3vec, ierr);CHKERRQ(ierr)
  end if

end subroutine MatGetSubQIMPESImmiscible


! ************************************************************************** !
!        pressure system extraction routines

subroutine MatGetSubQIMPES(a, ap, factors1Vec,  ierr, ctx)
  !
  ! extraction of the pressure system matrix of for the
  ! CPR preconditioner, and store the pivoting factors
  !
  ! Author:  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  ! 3x3 blocks ONLY

  implicit none


  Mat :: a, ap
  Vec :: factors1Vec
  PetscErrorCode :: ierr
  type(cpr_pc_type) :: ctx

  PetscInt, dimension(0:0) :: insert_rows
  ! the elements of the diagonal block:
  PetscReal :: aa,bb,cc,dd,ee,ff,gg,hh,ii
  ! misc workers:
  PetscReal :: det, fac0, fac1, fac2, sm
  PetscInt :: b, rws, cls, nblks, nblks_l, firstRow, cur_coldex, ncolblks, &
              firstrowdex, loopdex, i, j, numcols, numcols_keep
  PetscMPIInt :: rnk, r_st, r_nd



  ctx%vals = 0.d0
  ctx%insert_vals = 0.d0
  ctx%all_vals = 0.d0

  ctx%colIdx = 0
  ctx%colIdx_keep = 0
  ctx%insert_colIdx = 0

  call MPI_Comm_Rank(PETSC_COMM_WORLD, rnk, ierr)
  call MatGetOwnershipRange(a, r_st, r_nd, ierr); CHKERRQ(ierr)
  call MatGetBlockSize(a,b,ierr); CHKERRQ(ierr)
  call MatGetSize(A, rws, cls, ierr); CHKERRQ(ierr)
  if (rws /= cls) then
    ctx%option%io_buffer = 'MatGetSubQIMPES, given a nonsquare matrix'
    call PrintErrMsg(ctx%option)
  endif
  nblks = rws/b
  nblks_l = (r_nd-r_st)/b

  sm = 0.d0

  ! loop over the row blocks
  do i = 0,nblks_l-1

    firstRow = i*b + r_st

    ! a) extract first row
    call MatGetRow(a, firstRow, numcols, ctx%colIdx, ctx%vals, ierr); CHKERRQ(ierr)
    ! store vals since we have to put them back
    ! store colIdx this time as well
    do j = 0,numcols-1
      ctx%all_vals(0, j) = ctx%vals(j)
      ctx%colIdx_keep(j) = ctx%colIdx(j)
    end do
    ! b) we can get index of diagonal block here
    firstrowdex = -1
    do loopdex = 0,numcols-1,3
      if (ctx%colIdx(loopdex) == firstrow) then
        firstrowdex = loopdex
      endif
    enddo
    if (firstrowdex == -1) then
      ctx%option%io_buffer = 'MatGetSubQIMPES, cannot find diagonal entry, check matrix'
      call PrintErrMsg(ctx%option)
    endif
    aa = ctx%vals(firstrowdex)
    bb = ctx%vals(firstrowdex+1)
    cc = ctx%vals(firstrowdex+2)
    ! restore
    call MatRestoreRow(a, firstRow, numcols, ctx%colIdx, ctx%vals, ierr)
    CHKERRQ(ierr)

    ! c) second row
    call MatGetRow(a, firstRow+1, numcols, PETSC_NULL_INTEGER, ctx%vals, &
                  ierr); CHKERRQ(ierr)
    do j = 0,numcols-1
      ctx%all_vals(1, j) = ctx%vals(j)
    end do
    dd = ctx%vals(firstrowdex)
    ee = ctx%vals(firstrowdex+1)
    ff = ctx%vals(firstrowdex+2)
    call MatRestoreRow(a, firstRow+1, numcols, PETSC_NULL_INTEGER, &
                       ctx%vals, ierr); CHKERRQ(ierr)

    ! d) third row
    call MatGetRow(a, firstRow+2, numcols, PETSC_NULL_INTEGER, ctx%vals, &
                 ierr); CHKERRQ(ierr)
    do j = 0,numcols-1
      ctx%all_vals(2, j) =ctx%vals(j)
    end do
    gg = ctx%vals(firstrowdex)
    hh = ctx%vals(firstrowdex+1)
    ii = ctx%vals(firstrowdex+2)
    numcols_keep = numcols
    call MatRestoreRow(a, firstRow+2, numcols, PETSC_NULL_INTEGER, &
                       ctx%vals, ierr); CHKERRQ(ierr)

    ! e) factors
    sm = abs(aa)+abs(dd)+abs(gg)
    det = aa*(ee*ii - ff*hh) - bb*(dd*ii-ff*gg) + cc*(dd*hh-ee*gg)
    fac0 = sm*(ee*ii-ff*hh)/det
    fac1 = sm*(cc*hh-bb*ii)/det
    fac2 = sm*(bb*ff-cc*ee)/det

    ! f) store vectors
    call VecSetValue(factors1Vec, firstRow, fac0, INSERT_VALUES, ierr)
    CHKERRQ(ierr)
    call VecSetValue(factors1Vec, firstRow+1, fac1, INSERT_VALUES, ierr)
    CHKERRQ(ierr)
    call VecSetValue(factors1Vec, firstRow+2, fac2, INSERT_VALUES, ierr)
    CHKERRQ(ierr)

    ! g) prepare to set values
    insert_rows(0) = i + r_st/b
    ncolblks = numcols_keep/b
    do j = 0,ncolblks-1
      cur_coldex = j*b
      ctx%insert_colIdx(j) = ctx%colIdx_keep(cur_coldex)/b


      ctx%insert_vals(j) = fac0*ctx%all_vals(0, cur_coldex)
      ctx%insert_vals(j) = ctx%insert_vals(j) + fac1*ctx%all_vals(1,cur_coldex)
      ctx%insert_vals(j) = ctx%insert_vals(j) + fac2*ctx%all_vals(2,cur_coldex)
    end do

    ! h) set values
    call MatSetValues(ap, 1, insert_rows, ncolblks, &
                      ctx%insert_colIdx(0:ncolblks-1), &
                      ctx%insert_vals(0:ncolblks-1), INSERT_VALUES, ierr)
    CHKERRQ(ierr)


  end do
  call MatAssemblyBegin(ap,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
  call MatAssemblyEnd(ap,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)

  call VecAssemblyBegin(factors1Vec, ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(factors1vec, ierr);CHKERRQ(ierr)

end subroutine MatGetSubQIMPES

! ************************************************************************** !

subroutine MatGetSubQIMPES_var(a, ap, factors1Vec,  ierr, &
                              b, ctx                       )
  !
  ! extraction of the pressure system matrix of for the
  ! CPR preconditioner, and store the pivoting factors
  !
  ! Author:  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  ! for arbitrary block size b

  implicit none

  Mat :: a, ap
  Vec :: factors1vec
  PetscErrorCode :: ierr
  PetscInt :: b
  type(cpr_pc_type) :: ctx

  PetscInt, dimension(0:0) :: insert_rows
  PetscReal, dimension(0:b-1,0:b-1) :: diag_block
  PetscReal, dimension(0:b-1) :: local_factors
  PetscMPIInt :: rnk, r_st, r_nd
  PetscInt :: rws, cls, nblks, nblks_l, firstRow, cur_coldex, ncolblks, &
              firstrowdex, loopdex, i, j, k, numcols, numcols_keep
  PetscReal :: sm, offdiagsum, diagpart
  ! for lapack inversion routine:
  integer, dimension(1:b) :: ipiv
  PetscReal, dimension(1:b) :: work
  PetscInt :: lwork, invinfo, luinfo

  lwork = b

  ctx%vals = 0.d0
  ctx%insert_vals = 0.d0
  ctx%all_vals = 0.d0

  ctx%colIdx = 0
  ctx%colIdx_keep = 0
  ctx%insert_colIdx = 0

  call MPI_Comm_Rank(PETSC_COMM_WORLD, rnk, ierr)
  call MatGetOwnershipRange(a, r_st, r_nd, ierr); CHKERRQ(ierr)
  call MatGetSize(a, rws, cls, ierr); CHKERRQ(ierr)
  if (rws /= cls) then
    ctx%option%io_buffer = 'MatGetSubQIMPES, given a nonsquare matrix'
    call PrintErrMsg(ctx%option)
  endif

  nblks = rws/b
  nblks_l = (r_nd-r_st)/b

  ! loop over the row blocks
  do i = 0,nblks_l-1

    firstRow = i*b + r_st

    ! first row special treatment
    call MatGetRow(a, firstRow, numcols, ctx%colIdx, ctx%vals, ierr); CHKERRQ(ierr)

    ! store both values of row and the col indexs
    do k = 0,numcols-1
      ctx%all_vals(0, k) = ctx% vals(k)
      ctx%colIdx_keep(k) = ctx%colIdx(k)
    end do
    ! get index of diagonal block
    firstrowdex = -1
    do loopdex = 0,numcols-1,b
      if (ctx%colIdx(loopdex) == firstrow) then
        firstrowdex = loopdex
      endif
    enddo
    if (firstrowdex == -1) then
      ctx%option%io_buffer = 'MatGetSubQIMPES_var, cannot find diagonal entry, check matrix'
      call PrintErrMsg(ctx%option)
    endif
    numcols_keep = numcols
    ! restore first row
    call MatRestoreRow(a, firstRow, numcols, ctx%colIdx, ctx%vals, ierr)
    CHKERRQ(ierr)

    ! loop over remaining rows
    do j = 1,b-1
      call MatGetRow(a, firstRow+j, numcols, PETSC_NULL_INTEGER, ctx%vals, &
                     ierr); CHKERRQ(ierr)
      ! harvest values
      do k = 0,numcols-1
        ctx%all_vals(j, k) = ctx%vals(k)
      end do
      call MatRestoreRow(a, firstRow+j, numcols, PETSC_NULL_INTEGER, &
                         ctx%vals, ierr); CHKERRQ(ierr)
    enddo

    ! get inverse of block
    diag_block = ctx%all_vals(0:b-1, firstrowdex:firstrowdex+b-1)
    ! invert
    ! NOTE: this can fail and seemingly just return permutation
    call DGETRF(b, b, diag_block, b, IPIV, luinfo ) ! factorize first
    call DGETRI(b,   diag_block,       b,   IPIV,   WORK,     LWORK, invinfo)

    if (invinfo > 0) then
      ctx%option%io_buffer = 'MatGetSubQIMPES_var, singular diagonal block'
      call PrintErrMsg(ctx%option)
    endif
    ! scaling factor: the sum of abs of the first column of
    ! diagonal
    sm = 0.d0
    do j = 0,b-1
      !sm = sm + abs(ctx%all_vals(j, firstrowdex))
      sm = sm + abs(ctx%all_vals(j, firstrowdex+ctx%exrow_offset))
    end do
    ! factors: take the top row of the inverse block
    ! and scale
    do j = 0,b-1
        local_factors(j) =  sm*diag_block(0, j) ! diag block has been
                                                ! replaced with inverse by this point

    end do

    offdiagsum = 0.d0
    diagpart   = 0.d0

    ! prepare to set values
    insert_rows(0) = i + r_st/b
    ncolblks = numcols_keep/b

    do j = 0,ncolblks-1
      cur_coldex = j*b
      ctx%insert_colIdx(j) = ctx%colIdx_keep(cur_coldex)/b

      ctx%insert_vals(j) = 0.d0
      do k = 0,b-1
        ctx%insert_vals(j) = ctx%insert_vals(j) + local_factors(k)*ctx%all_vals(k, cur_coldex)
      enddo

      if (ctx%insert_colIdx(j) == insert_rows(0)) then
        diagpart = abs(ctx%insert_vals(j))
      else
        offdiagsum = offdiagsum + abs(ctx%insert_vals(j))
      endif

    end do


    do j = 0,b-1
       call VecSetValue(factors1Vec, firstRow+j, local_factors(j), &
                         INSERT_VALUES, ierr); CHKERRQ(ierr)
    end do


    ! set values
    call MatSetValues(ap, 1, insert_rows, ncolblks, &
                      ctx%insert_colIdx(0:ncolblks-1), &
                      ctx%insert_vals(0:ncolblks-1), INSERT_VALUES, ierr)
     CHKERRQ(ierr)

  end do
  call MatAssemblyBegin(ap,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)
  call MatAssemblyEnd(ap,MAT_FINAL_ASSEMBLY,ierr); CHKERRQ(ierr)

  call VecAssemblyBegin(factors1Vec, ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(factors1vec, ierr);CHKERRQ(ierr)

end subroutine MatGetSubQIMPES_var

! ************************************************************************** !

subroutine QIRHS(factors, worker, r, rhat, ierr)
  !
  ! extract the RHS of the pressure system for the CPR
  ! preconditioner, given pivoting factors (factors)
  ! and full system rhs (r).
  ! output is rhat.
  !
  ! Author:  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  implicit none

  Vec :: factors, worker, r, rhat
  PetscErrorCode :: ierr

  PetscInt :: b, k

  ! Explanation:
  ! residual = [p1,s1,t1,p2,s2,t2, ... pn,sn,tn]
  ! as in pressure saturation temperature
  ! factors = [fac0_1,fac1_1,fac2_1,fac0_2,fac1_2,fac2_3 ... ]
  ! fac could be j_ps*inv(j_ss)
  ! multipily factor to residual element by element
  ! worker =  [p1*fac0_1,s1*fac1_1,t1*fac2_1 ... ]
  ! then VecStrideGather collapses the results applying addition
  ! rhat[1] = p1*fac0_1 + s1*fac1_1 + t1*fac2_1
  ! rhat[2] = p2*fac0_2 + s2*fac1_2 + t2*fac2_2
  ! so that size(rhat) = size(r)/b
  ! - Heeho Park

  call VecPointwiseMult(worker, factors, r, ierr);CHKERRQ(ierr)
  call VecGetBlockSize(worker,b,ierr); CHKERRQ(ierr)
  k = 0
  call VecStrideGather(worker, k, rhat, INSERT_VALUES, ierr);CHKERRQ(ierr)
  do k = 1,b-1
   call VecStrideGather(worker, k, rhat, ADD_VALUES, ierr);CHKERRQ(ierr)
  end do

end subroutine QIRHS

! end of pressure system extraction routines
! ************************************************************************** !

! ************************************************************************** !
!       misc routines

subroutine MatGetMaxRowCount(a, mx, ierr)
  !
  ! loop over rows of matrix, return the greatest number
  ! of nonzeros in a row
  !
  ! Author:  Daniel Stone
  ! Date: Oct 2017 - March 2018
  !

  implicit none

  Mat :: a
  PetscInt :: mx
  PetscErrorCode :: ierr

  PetscInt :: mx_loc, r_st, r_nd, i, numcols

  call MatGetOwnershipRange(a, r_st, r_nd, ierr); CHKERRQ(ierr)

  mx = 0
  mx_loc = 0

  do i = r_st,r_nd-1
    call MatGetRow(a, i, numcols, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, &
                   ierr); CHKERRQ(ierr)
    if (numcols > mx_loc) then
      mx_loc = numcols
    endif
    call MatRestoreRow(a, i, numcols, PETSC_NULL_INTEGER, &
                       PETSC_NULL_SCALAR, ierr); CHKERRQ(ierr)
  end do

  call MPI_Allreduce(mx_loc, mx, ONE_INTEGER_MPI, MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)


end subroutine MatGetMaxRowCount

! end of misc routines
! ************************************************************************** !

end module CPR_Preconditioner_module
