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
            NumCompare_towg_bo

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
