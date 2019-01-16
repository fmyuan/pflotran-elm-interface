module Derivative_tests_module

! Here are some simple routines that are never meant to be used in production runs
! but are very useful in debugging and development.

! The idea is to take the existing structure available for computing numerical
! derivatives (i.e. nflowdof + 1 auxvar objects, wih perturbations in each of the 
! solution variables), and use that to compute numerical derivatives for the 
! auxiliary variables.

! These numerical derivatives are not used anywhere, but are compared to 
! corresponding analytical derivatives of the aux variables that have
! already been computed and stored.

! This allows sanity checks and/or validation of the analytical derivatives
! if they have been computed and stored on the auxvar variable.


#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use AuxVars_Flow_module
  use AuxVars_TOilIms_module
  use AuxVars_TOWG_module

  implicit none

  private

  public :: NumCompare_toil, &
            NumCompare_towg_bo, &
            NumCompare_tl3p, &
            NumCompare_tl4p, &
            Num_as_alyt_tl4p

contains

! ************************************************************************** !

subroutine NumCompare_toil(nphase,ndof,auxvars,option)
  ! intention: given an array of auxvar objects, do differencing 
  ! and compare to stored analytical derivatives

  ! provides some validation/sanity control for analytical 
  ! derivative computation of auxvars

  ! Daniel Stone, July 2018
  implicit none
  PetscInt :: nphase,ndof
  type(option_type) :: option
  type(auxvar_toil_ims_type) :: auxvars(0:)

  PetscInt :: iphase,idof
  PetscReal :: pert,p_unpert,p_pert,nderiv,aderiv,diff,rdiff
  PetscReal :: atol,rtol
  PetscInt :: probs

  atol = flow_aux_debug_tol
  rtol = flow_aux_debug_reltol

  print *
  print *, "NumCompare TOil"
  print *, "Comparing numerical to analytical derivatives of aux variables in one cell."
  print *, "Displaying potential problems if difference is > ", atol, ","
  print *, "or if relative difference is > ", rtol 
  print *

  probs = 0

  !! ********* from auxvars flow *********

  ! pres:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%pres(iphase)
     ! get perturbed value
      p_pert = auxvars(idof)%pres(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_pres(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "pres:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        probs = probs + 1
        print *
      endif
    enddo
  enddo


  ! sat:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%sat(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%sat(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_sat(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "sat:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        probs = probs + 1
        print *
      endif
    enddo
  enddo


  ! pc - needs special treatment because not full (nphase) size:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase-1
      ! get unperturbed value
      p_unpert = auxvars(0)%pc(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%pc(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_pc(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "pc:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        probs = probs + 1
        print *
      endif
    enddo
  enddo

  ! den:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%den(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%den(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_den(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "den:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        probs = probs + 1
        print *
      endif
    enddo
  enddo


  ! den kg:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%den_kg(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%den_kg(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_den_kg(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "den_kg:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        probs = probs + 1
        print *
      endif
    enddo
  enddo


  ! mob:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%mobility(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%mobility(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_mobility(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
      !if (diff>atol .OR. rdiff>rtol .OR. auxvars(0)%sat(1) == 0 .OR. auxvars(0)%sat(1) == 1 ) then
        print *, "mob:"
        print *,  "sats: ", auxvars(0)%sat
        print *, "pert is ", pert
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        probs = probs + 1
        print *
      endif
    enddo
  enddo


  ! poro:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,1
      ! get unperturbed value
      p_unpert = auxvars(0)%effective_porosity
      ! get perturbed value
      p_pert = auxvars(idof)%effective_porosity

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_por(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "poro (note phase is meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        probs = probs + 1
        print *
      endif
    enddo
  enddo

  !! ********* end of from auxvars flow *********

  !! ********* from auxvars energy flow *********

  ! H:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%H(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%H(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_H(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "H:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        probs = probs + 1
        print *
      endif
    enddo
  enddo

  ! U:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%U(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%U(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_U(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "U:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        probs = probs + 1
        print *
      endif
    enddo
  enddo
  !! ********* end of from auxvars energy flow *********

  print *, "NumCompare TOil, total possible problems here: ", probs

end subroutine NumCompare_toil

! ************************************************************************** !

subroutine NumCompare_towg_bo(nphase,ndof,auxvars,option,&
                              dof_op,dof_osat,dof_gsat,dof_temp,&
                              isSat)
  ! intention: given an array of auxvar objects, do differencing 
  ! to compare analytical derivatives

  use PFLOTRAN_Constants_module
  implicit none
  PetscInt :: nphase,ndof
  type(auxvar_towg_type) :: auxvars(0:)
  type(option_type) :: option

  PetscInt :: iphase,idof
  PetscReal :: pert,p_unpert,p_pert,nderiv,aderiv,diff,rdiff
  PetscReal :: atol,rtol
  PetscInt :: probs

  PetscInt :: dof_op,dof_osat,dof_gsat,dof_temp
  PetscBool :: isSat


  atol = flow_aux_debug_tol
  rtol = flow_aux_debug_reltol

#if 0
  print *
  print *, "NumCompare TOWG BO"
  print *, "Comparing numerical to analytical derivatives of aux variables in one cell."
  print *, "Displaying potential problems if difference is > ", atol, ","
  print *, "or if relative difference is > ", rtol 
  if (isSat) then
    print *, "This cell is in saturated state."
  else
    print *, "This cell is in UNsaturated state."
  endif
  print *
#endif

  probs = 0

  !! *********  from auxvars BO*********
  ! xo
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%bo%xo
      ! get perturbed value
      p_pert = auxvars(idof)%bo%xo

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%bo%D_xo(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "xo (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo

  ! xg
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%bo%xg
      ! get perturbed value
      p_pert = auxvars(idof)%bo%xg

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%bo%D_xg(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "xg (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo

  !! *********  end offrom auxvars BO*********


  !! ********* from auxvars flow *********

  ! pres:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%pres(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%pres(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_pres(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "pres:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo


  ! sat:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%sat(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%sat(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_sat(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "sat:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo


  ! pc - needs special treatment because not full (nphase) size:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase-1
      ! get unperturbed value
      p_unpert = auxvars(0)%pc(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%pc(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_pc(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "pc:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo

  ! den:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%den(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%den(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_den(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "den:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo


  ! den kg:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%den_kg(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%den_kg(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_den_kg(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "den_kg:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo


  ! mob:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%mobility(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%mobility(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_mobility(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "mob:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,  "sats: ", auxvars(0)%sat
        print *,  "mob0: ", p_unpert
        print *,  "mob pert: ", p_pert
        print *,  "mob diff: ", p_pert - p_unpert
        print *,  "pert in dof is: ", pert
        print *
        probs = probs + 1
      endif
    enddo
  enddo

  ! poro:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,1
      ! get unperturbed value
      p_unpert = auxvars(0)%effective_porosity
      ! get perturbed value
      p_pert = auxvars(idof)%effective_porosity

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_por(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "poro:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo

  !! ********* end of from auxvars flow *********

  !! ********* from auxvars energy flow *********

  ! H:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%H(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%H(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_H(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "H:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo

  ! U:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%U(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%U(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%D_U(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "U:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo
  !! ********* end of from auxvars energy flow *********

  !print *, "NumCompare TOWG BO, total possible problems here: ", probs
  !print *, 

  if (probs > 0) then
    print *, "NumCompare TOWG BO (atol: ", atol, ", rtol: ", rtol
    print *, ", saturated? ", isSat, "), probs: ", probs
    print *, "more than 0 problems"
  endif


end subroutine NumCompare_towg_bo

! ************************************************************************** !

subroutine NumCompare_tl3p(nphase,ndof,auxvars,option,&
                              dof_op,dof_osat,dof_gsat,dof_temp)

  ! intention: given an array of auxvar objects, do differencing 
  ! to compare analytical derivatives

  !!!! HACK FOR DEVELOPMENT/TESTING  - we will in fact OVERWRITE
  !!!! the analytical derivative arrays, until the implementation
  !!!! to fill them up has been implemented

  use PFLOTRAN_Constants_module
  implicit none
  PetscInt :: nphase,ndof
  type(auxvar_towg_type),intent(INOUT) :: auxvars(0:)
  type(option_type) :: option

  PetscInt :: iphase,idof
  PetscReal :: pert,p_unpert,p_pert,nderiv,aderiv,diff,rdiff
  PetscReal :: atol,rtol
  PetscInt :: probs

  PetscInt :: dof_op,dof_osat,dof_gsat,dof_temp
  PetscBool :: isSat


  atol = flow_aux_debug_tol
  rtol = flow_aux_debug_reltol

  probs = 0

  !! *********  TL INTERMEDIATES************
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tl3TEST%krh
      ! get perturbed value
      p_pert = auxvars(idof)%tl3TEST%krh

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%tl3TEST%D_krh(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krh, note phase meaningless:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,  "sats: ", auxvars(0)%sat
        print *,  "krh0: ", p_unpert
        print *,  "krh pert: ", p_pert
        print *,  "krh diff: ", p_pert - p_unpert
        print *,  "pert in dof is: ", pert
        print *
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tl3TEST%denotl
      ! get perturbed value
      p_pert = auxvars(idof)%tl3TEST%denotl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%tl3TEST%D_denotl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "denotl, note phase meaningless:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tl3TEST%dengtl
      ! get perturbed value
      p_pert = auxvars(idof)%tl3TEST%dengtl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%tl3TEST%D_dengtl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "dengtl, note phase meaningless:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tl3TEST%viscotl
      ! get perturbed value
      p_pert = auxvars(idof)%tl3TEST%viscotl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%tl3TEST%D_viscotl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "viscotl, note phase meaningless:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tl3TEST%viscgtl
      ! get perturbed value
      p_pert = auxvars(idof)%tl3TEST%viscgtl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%tl3TEST%D_viscgtl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "viscgtl, note phase meaningless:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo


  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tl3TEST%krotl
      ! get perturbed value
      p_pert = auxvars(idof)%tl3TEST%krotl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%tl3TEST%D_krotl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krotl, note phase meaningless:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tl3TEST%krgtl
      ! get perturbed value
      p_pert = auxvars(idof)%tl3TEST%krgtl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%tl3TEST%D_krgtl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krgtl, note phase meaningless:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo


  !! ********* END  TL INTERMEDIATES*******

  !! *********  from TL specifically*********
  ! den_oil_eff_kg:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tl%den_oil_eff_kg
      ! get perturbed value
      p_pert = auxvars(idof)%tl%den_oil_eff_kg

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%tl%D_den_oil_eff_kg(idof) = nderiv
#endif
      ! analytical derivative
      aderiv = auxvars(0)%tl%D_den_oil_eff_kg(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "den_oil_eff_kg, note phase meaningless:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo

  ! den_gas_eff_kg:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tl%den_gas_eff_kg
      ! get perturbed value
      p_pert = auxvars(idof)%tl%den_gas_eff_kg

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%tl%D_den_gas_eff_kg(idof) = nderiv
#endif
      ! analytical derivative
      aderiv = auxvars(0)%tl%D_den_gas_eff_kg(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "den_gas_eff_kg, note phase meaningless:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo
  !! *********  end offrom TL specifically*********


  !! ********* from auxvars flow *********

  ! pres:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%pres(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%pres(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%D_pres(iphase,idof) = nderiv
#endif

      ! analytical derivative
      aderiv = auxvars(0)%D_pres(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "pres:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif

    enddo
  enddo


  ! sat:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%sat(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%sat(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%D_sat(iphase,idof) = nderiv
#endif

      ! analytical derivative
      aderiv = auxvars(0)%D_sat(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "sat:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif

    enddo
  enddo


  ! pc - needs special treatment because not full (nphase) size:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase-1
      ! get unperturbed value
      p_unpert = auxvars(0)%pc(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%pc(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%D_pc(iphase,idof) = nderiv
#endif

      ! analytical derivative
      aderiv = auxvars(0)%D_pc(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "pc:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo

  ! den:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%den(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%den(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%D_den(iphase,idof) = nderiv
#endif

      ! analytical derivative
      aderiv = auxvars(0)%D_den(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "den:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo


  ! den kg:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%den_kg(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%den_kg(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%D_den_kg(iphase,idof) = nderiv
#endif
      ! analytical derivative
      aderiv = auxvars(0)%D_den_kg(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "den_kg:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo


  ! mob:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%mobility(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%mobility(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%D_mobility(iphase,idof) = nderiv
#endif
      ! analytical derivative
      aderiv = auxvars(0)%D_mobility(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "mob:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,  "sats: ", auxvars(0)%sat
        print *,  "mob0: ", p_unpert
        print *,  "mob pert: ", p_pert
        print *,  "mob diff: ", p_pert - p_unpert
        print *,  "pert in dof is: ", pert
        print *
        probs = probs + 1
      endif
    enddo
  enddo

  ! poro:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,1
      ! get unperturbed value
      p_unpert = auxvars(0)%effective_porosity
      ! get perturbed value
      p_pert = auxvars(idof)%effective_porosity

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%D_por(idof) = nderiv
#endif
      ! analytical derivative
      aderiv = auxvars(0)%D_por(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "poro:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo

  !! ********* end of from auxvars flow *********

  !! ********* from auxvars energy flow *********

  ! H:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%H(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%H(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%D_H(iphase,idof) = nderiv
#endif
      ! analytical derivative
      aderiv = auxvars(0)%D_H(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "H:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo

  ! U:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%U(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%U(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

#if 0
      !!! assign
      auxvars(0)%D_U(iphase,idof) = nderiv
#endif
      ! analytical derivative
      aderiv = auxvars(0)%D_U(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "U:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo
  !! ********* end of from auxvars energy flow *********


  if (probs > 0) then
    print *, "NumCompare TL3P (atol: ", atol, ", rtol: ", rtol, ")"
    print *, "more than 0 problems"
  endif


end subroutine NumCompare_tl3p

! ************************************************************************** !

subroutine Num_as_alyt_tl4p(nphase,ndof,auxvars,option,&
                              dof_op,dof_osat,dof_gsat,dof_temp,&
                              isSat)
  ! intention: given an array of auxvar objects, do differencing 
  ! to populate auxvars derivative members with numerical derivatives

  use PFLOTRAN_Constants_module
  implicit none
  PetscInt :: nphase,ndof
  type(auxvar_towg_type) :: auxvars(0:)
  type(option_type) :: option

  PetscInt :: iphase,idof
  PetscReal :: pert,p_unpert,p_pert,nderiv,aderiv,diff,rdiff

  PetscInt :: dof_op,dof_osat,dof_gsat,dof_temp
  PetscBool :: isSat


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if 0
  ! pres:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%pres(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%pres(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      !!! assign
      auxvars(0)%D_pres(iphase,idof) = nderiv

    enddo
  enddo



  ! pc - needs special treatment because not full (nphase) size:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase-1
      ! get unperturbed value
      p_unpert = auxvars(0)%pc(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%pc(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      !!! assign
      auxvars(0)%D_pc(iphase,idof) = nderiv

    enddo
  enddo
#endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! den:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%den(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%den(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      !!! assign
      auxvars(0)%D_den(iphase,idof) = nderiv

#if 0
      ! analytical derivative
      aderiv = auxvars(0)%D_den(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "den:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
#endif
    enddo
  enddo


  ! den kg:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%den_kg(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%den_kg(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      !!! assign
      auxvars(0)%D_den_kg(iphase,idof) = nderiv
      

#if 0
      ! analytical derivative
      aderiv = auxvars(0)%D_den_kg(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "den_kg:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
#endif
    enddo
  enddo


  ! mob:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%mobility(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%mobility(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      !!! assign
      auxvars(0)%D_mobility(iphase,idof) = nderiv

#if 0
      ! analytical derivative
      aderiv = auxvars(0)%D_mobility(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "mob:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,  "sats: ", auxvars(0)%sat
        print *,  "mob0: ", p_unpert
        print *,  "mob pert: ", p_pert
        print *,  "mob diff: ", p_pert - p_unpert
        print *,  "pert in dof is: ", pert
        print *
        probs = probs + 1
      endif
#endif
    enddo
  enddo

  ! poro:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,1
      ! get unperturbed value
      p_unpert = auxvars(0)%effective_porosity
      ! get perturbed value
      p_pert = auxvars(idof)%effective_porosity

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      !!! assign
      auxvars(0)%D_por(idof) = nderiv
#if 0

      ! analytical derivative
      aderiv = auxvars(0)%D_por(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "poro:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
#endif
    enddo
  enddo

  !! ********* end of from auxvars flow *********

  !! ********* from auxvars energy flow *********

  ! H:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%H(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%H(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      !!! assign
      auxvars(0)%D_H(iphase,idof) = nderiv

#if 0
      ! analytical derivative
      aderiv = auxvars(0)%D_H(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "H:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
#endif
    enddo
  enddo

  ! U:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%U(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%U(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      !!! assign
      auxvars(0)%D_U(iphase,idof) = nderiv
#if 0

      ! analytical derivative
      aderiv = auxvars(0)%D_U(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "U:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
#endif
    enddo
  enddo
  !! ********* end of from auxvars energy flow *********

end subroutine Num_as_alyt_tl4p

! ************************************************************************** !

subroutine NumCompare_tl4p(nphase,ndof,auxvars,option,&
                              dof_op,dof_osat,dof_gsat,dof_temp,&
                              isSat)
  ! intention: given an array of auxvar objects, do differencing 
  ! to compare analytical derivatives

  use PFLOTRAN_Constants_module
  implicit none
  PetscInt :: nphase,ndof
  type(auxvar_towg_type) :: auxvars(0:)
  type(option_type) :: option

  PetscInt :: iphase,idof
  PetscReal :: pert,p_unpert,p_pert,nderiv,aderiv,diff,rdiff
  PetscReal :: atol,rtol
  PetscInt :: probs

  PetscInt :: dof_op,dof_osat,dof_gsat,dof_temp
  PetscBool :: isSat


  atol = flow_aux_debug_tol
  rtol = flow_aux_debug_reltol

  probs = 0

  !! *********  from auxvars BO*********
  ! xo
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%bo%xo
      ! get perturbed value
      p_pert = auxvars(idof)%bo%xo

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%bo%D_xo(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "xo (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo

  ! xg
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%bo%xg
      ! get perturbed value
      p_pert = auxvars(idof)%bo%xg

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%bo%D_xg(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "xg (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    !enddo
  enddo

  !! *********  end offrom auxvars BO*********


  !! ********* from auxvars flow *********

  ! pres:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%pres(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%pres(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%D_pres(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "pres:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo


  ! sat:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%sat(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%sat(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%D_sat(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "sat:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *
        probs = probs + 1
      endif
    enddo
  enddo


  ! pc - needs special treatment because not full (nphase) size:
  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase-1
      ! get unperturbed value
      p_unpert = auxvars(0)%pc(iphase)
      ! get perturbed value
      p_pert = auxvars(idof)%pc(iphase)

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert

      ! analytical derivative
      aderiv = auxvars(0)%D_pc(iphase,idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "pc:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *, "sats: "
        print *,auxvars(0)%sat
        print *, "fm: ", auxvars(0)%tlt%fm
        print *, auxvars(0)%tlT%D_fm
        print *
        probs = probs + 1
      endif
    enddo
  enddo

  !! ***** end of from auxvars flow *********


  !! ***************  from tl4p intermediates  *********
  if (auxvars(0)%has_TL_test_object) then
      do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    do iphase=1,nphase
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%viso
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%viso

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_viso(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "viso:"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlt%fm
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%fm

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_fm(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "fm (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%krom
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%krom

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_krom(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krom (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *, "sats:"
        print *, auxvars(0)%sat
        print *, "fm: ", auxvars(0)%tlT%fm
        print *, "dfm: "
        print *, auxvars(0)%tlT%D_fm(:)
#if 0
        if (idof == 2) then
          print *, "here"
        endif
#endif
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%krvm
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%krvm

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_krvm(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krvm (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%krgm
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%krgm

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_krgm(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krgm (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%krsm
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%krsm

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_krsm(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krsm (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%krotl
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%krotl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_krotl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krotl (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
#if 0
        if (idof ==  2 .AND. abs(nderiv) > 0.d0) then
          print *, "here"
        endif
#endif
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%kroi
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%kroi

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_kroi(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "kroi (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%krog
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%krog

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_krog(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krog (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%krow
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%krow

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_krow(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krow (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%krgtl
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%krgtl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_krgtl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krgtl (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%krstl
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%krstl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_krstl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "krstl (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%viscotl
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%viscotl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_viscotl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "viscotl (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%viscgtl
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%viscgtl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_viscgtl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "viscgtl (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%viscstl
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%viscstl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_viscstl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "viscstl (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%denotl
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%denotl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_denotl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "denotl (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%dengtl
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%dengtl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_dengtl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "dengtl (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%denstl
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%denstl

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_denstl(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "denstl (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%uoil
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%uoil

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_uoil(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "uoil (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo

  do idof=1,ndof
    ! get perturbation for this dof variable
    pert = auxvars(idof)%pert
    !do iphase=1,nphase
    iphase = 1
      ! get unperturbed value
      p_unpert = auxvars(0)%tlT%uvap
      ! get perturbed value
      p_pert = auxvars(idof)%tlT%uvap

      ! numerical derivative
      nderiv = (p_pert-p_unpert)/pert
      ! analytical derivative
      aderiv = auxvars(0)%tlT%D_uvap(idof)

      ! difference:
      diff = abs(aderiv-nderiv)
      rdiff = diff/abs(nderiv)

      if (diff>atol .OR. rdiff>rtol) then
        print *, "uvap (note phase meaningless):"
        call NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)
        print *,
        probs = probs + 1
      endif
    !enddo
  enddo
  endif

  !! ********* end of from tl4p intermediates  *********

  !print *, "NumCompare TOWG BO, total possible problems here: ", probs
  !print *, 

  if (probs > 0) then
    print *, "NumCompare tl4p (atol: ", atol, ", rtol: ", rtol
    print *, ", saturated? ", isSat, "), probs: ", probs
    print *, "more than 0 problems"
  endif


end subroutine NumCompare_tl4p

! ************************************************************************** !

subroutine NumCompareOutput(idof,iphase,nderiv,aderiv,diff,rdiff)

  implicit none
  PetscInt :: idof,iphase
  PetscReal :: nderiv,aderiv,diff,rdiff

  print *, "phase, ", iphase, ", deriv w.r.t. variable number ", idof
  print *, "alyt: ", aderiv, "; num: ", nderiv
  print *, "; DIFFERENCE (abs,rel): ", diff, rdiff

end subroutine NumCompareOutput

! ************************************************************************** !

end module Derivative_tests_module
