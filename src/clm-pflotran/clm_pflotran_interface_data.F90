!#include <preproc.h>

module clm_pflotran_interface_data

  ! !USES:
  !  use shr_kind_mod, only     : r8 => shr_kind_r8
  !!  use abortutils, only       : endrun
  !!  use clm_varctl      , only : iulog
  !!  use spmdMod         , only : masterproc
  !
  ! !PUBLIC TYPES:
  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

  type, public :: clm_pflotran_data

     !
     ! Time invariant data:
     !
     ! (i) Soil properties -
     ! Global vectors
     Vec :: hksat_x, hksat_y, hksat_z
     Vec :: succat
     Vec :: watsat
     Vec :: bsw

     ! Local for CLM
     Vec :: hksat_x_clmloc, hksat_y_clmloc, hksat_z_clmloc
     Vec :: succat_clmloc
     Vec :: watsat_clmloc
     Vec :: bsw_clmloc

     ! Local for PFLOTRAN
     Vec :: hksat_x_pfloc, hksat_y_pfloc, hksat_z_pfloc
     Vec :: succat_pfloc
     Vec :: watsat_pfloc
     Vec :: bsw_pfloc

     ! (ii) Initial conditions -
     Vec :: zwt_2d, zwt_2d_clmloc, zwt_2d_pfloc

     !
     ! Time variant data
     !
     ! (i) Sink/Source of water
     Vec :: qflx, qflx_clmloc, qflx_pfloc

     ! (ii) Saturation
     Vec :: sat, sat_clmloc, sat_pfloc

     !
     PetscInt :: clm_num_local
     !PetscInt :: clm_num_ghost
     PetscInt :: pf_num_local
     PetscInt :: clm_num_docells
     PetscInt :: pf_num_docells
	 
	 
	 ! 
	 Vec :: sat_clm, sat_pf


  end type clm_pflotran_data

  public :: clm_pf_data_create, &
       clm_pf_data_allocate_memory

contains
  function clm_pf_data_create()

    implicit none

    type(clm_pflotran_data), pointer :: clm_pf_data_create
    type(clm_pflotran_data), pointer :: clm_pf_data

    allocate(clm_pf_data)

    clm_pf_data%hksat_x = 0
    clm_pf_data%hksat_y = 0
    clm_pf_data%hksat_z = 0
    clm_pf_data%succat  = 0
    clm_pf_data%watsat  = 0
    clm_pf_data%bsw     = 0
    clm_pf_data%zwt_2d  = 0
    clm_pf_data%qflx    = 0
    clm_pf_data%sat     = 0

    clm_pf_data%hksat_x_clmloc = 0
    clm_pf_data%hksat_y_clmloc = 0
    clm_pf_data%hksat_z_clmloc = 0
    clm_pf_data%succat_clmloc  = 0
    clm_pf_data%watsat_clmloc  = 0
    clm_pf_data%bsw_clmloc     = 0
    clm_pf_data%zwt_2d_clmloc  = 0
    clm_pf_data%qflx_clmloc    = 0
    clm_pf_data%sat_clmloc     = 0

    clm_pf_data%hksat_x_pfloc = 0
    clm_pf_data%hksat_y_pfloc = 0
    clm_pf_data%hksat_z_pfloc = 0
    clm_pf_data%succat_pfloc  = 0
    clm_pf_data%watsat_pfloc  = 0
    clm_pf_data%bsw_pfloc     = 0
    clm_pf_data%zwt_2d_pfloc  = 0
    clm_pf_data%qflx_pfloc    = 0
    clm_pf_data%sat_pfloc     = 0

    clm_pf_data%clm_num_local   = -1
    !clm_pf_data%clm_num_ghost = -1
    clm_pf_data%pf_num_local    = -1
    clm_pf_data%clm_num_docells = -1
    clm_pf_data%pf_num_docells  = -1

    clm_pf_data_create => clm_pf_data

  end function clm_pf_data_create

  subroutine clm_pf_data_allocate_memory(data, mycomm)

    implicit none

    PetscErrorCode :: ierr
    PetscMPIInt    :: mycomm, rank
    PetscViewer    :: viewer
    type(clm_pflotran_data), pointer :: data

    call MPI_Comm_rank(MPI_COMM_WORLD,rank, ierr)

    !if (data%clm_num_local.lt.0) then
    !   if (masterproc) write(iulog,*) 'create_clm_pf_data_allocate memory:', &
    !        ' clm_num_local = ',data%clm_num_local
    !   call endrun()
    !endif

    !if (data%pf_num_docells.lt.0) then
    !   if (masterproc) write(iulog,*) 'create_clm_pf_data_allocate memory:', &
    !        ' pf_num_docells = ',data%pf_num_docells
    !   call endrun()
    !endif


    !write(*,*), '>>> rank=',rank, 'clm_num_local=',data%clm_num_local, &
    !     'pf_num_local   =',data%pf_num_local, &
    !     'clm_num_docells=',data%clm_num_docells, &
    !     'pf_num_docells =',data%pf_num_docells

    ! Create MPI Vectors
    call VecCreateMPI(mycomm,data%clm_num_local,PETSC_DETERMINE,data%hksat_x,ierr)
    call VecSetBlockSize(data%hksat_x, 1, ierr)
    call VecSet(data%hksat_x, 0.0d0)
    call VecDuplicate(data%hksat_x, data%hksat_y, ierr)
    call VecDuplicate(data%hksat_x, data%hksat_z, ierr)
    call VecDuplicate(data%hksat_x, data%succat,  ierr)
    call VecDuplicate(data%hksat_x, data%watsat,  ierr)
    call VecDuplicate(data%hksat_x, data%bsw,     ierr)
    call VecDuplicate(data%hksat_x, data%zwt_2d,  ierr)
    call VecDuplicate(data%hksat_x, data%qflx,    ierr)

    call VecCreateMPI(mycomm,data%pf_num_local,PETSC_DETERMINE,data%sat,ierr)
    call VecSetBlockSize(data%sat, 1, ierr)


    ! Create Sequential Vectors
    call VecCreateSeq(PETSC_COMM_SELF,data%clm_num_local,data%hksat_x_clmloc,ierr)
    call VecSetBlockSize(data%hksat_x_clmloc, 1, ierr)

    call VecDuplicate(data%hksat_x_clmloc, data%hksat_y_clmloc, ierr)
    call VecDuplicate(data%hksat_x_clmloc, data%hksat_z_clmloc, ierr)
    call VecDuplicate(data%hksat_x_clmloc, data%succat_clmloc,  ierr)
    call VecDuplicate(data%hksat_x_clmloc, data%watsat_clmloc,  ierr)
    call VecDuplicate(data%hksat_x_clmloc, data%bsw_clmloc,     ierr)
    call VecDuplicate(data%hksat_x_clmloc, data%zwt_2d_clmloc,  ierr)
    call VecDuplicate(data%hksat_x_clmloc, data%qflx_clmloc,    ierr)

    call VecCreateSeq(PETSC_COMM_SELF,data%clm_num_docells,data%sat_clmloc,ierr)
    call VecSetBlockSize(data%sat_clmloc, 1,ierr)

    call VecCreateSeq(PETSC_COMM_SELF,data%pf_num_docells,data%hksat_x_pfloc,ierr)
    call VecSetBlockSize(data%hksat_x_pfloc, 1, ierr)

    call VecDuplicate(data%hksat_x_pfloc, data%hksat_y_pfloc, ierr)
    call VecDuplicate(data%hksat_x_pfloc, data%hksat_z_pfloc, ierr)
    call VecDuplicate(data%hksat_x_pfloc, data%succat_pfloc,  ierr)
    call VecDuplicate(data%hksat_x_pfloc, data%watsat_pfloc,  ierr)
    call VecDuplicate(data%hksat_x_pfloc, data%bsw_pfloc,     ierr)
    call VecDuplicate(data%hksat_x_pfloc, data%zwt_2d_pfloc,  ierr)
    call VecDuplicate(data%hksat_x_pfloc, data%qflx_pfloc,    ierr)

    call VecCreateSeq(PETSC_COMM_SELF,data%pf_num_local,data%sat_pfloc,ierr)
    call VecSetBlockSize(data%sat_pfloc,1,ierr)


    write(*,*), 'size(hksat_x_clmloc) = ',data%clm_num_local
    write(*,*),	'size(hksat_x       ) = ',data%clm_num_local
	  write(*,*), 'size(hksat_x_pfloc ) = ',data%pf_num_docells
    write(*,*), 'size(sat_pfloc     ) = ',data%pf_num_local
    write(*,*),	'size(sat           ) = ',data%pf_num_local
    write(*,*), 'size(sat_clmloc    ) = ',data%clm_num_docells


    call VecCreateMPI(mycomm,data%clm_num_local, PETSC_DETERMINE,data%sat_clm,ierr)
	call VecCreateMPI(mycomm,data%pf_num_local , PETSC_DETERMINE,data%sat_pf ,ierr)

  end subroutine clm_pf_data_allocate_memory

end module clm_pflotran_interface_data
