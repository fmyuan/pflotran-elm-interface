module clm_pflotran_interface_data

  implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
  private

  type, public :: clm_pflotran_interface_data_type

  ! Time invariant data:
  
  ! (i) Soil properties -
  ! Local for CLM  - mpi vectors
  Vec :: hksat_x_clm
  Vec :: hksat_y_clm
  Vec :: hksat_z_clm
  Vec :: sucsat_clm
  Vec :: watsat_clm
  Vec :: bsw_clm
  Vec :: press_clm

  ! Local for PFLOTRAN - seq. vec
  Vec :: hksat_x_pf
  Vec :: hksat_y_pf
  Vec :: hksat_z_pf
  Vec :: sucsat_pf
  Vec :: watsat_pf
  Vec :: bsw_pf
  Vec :: press_pf

  ! Time variant data
  
  ! (i) Sink/Source of water
  !!Vec :: qflx
  Vec :: qflx_clm   ! mpi vec
  Vec :: qflx_pf    ! seq vec

  ! (ii) Saturation
  Vec :: sat_clm    ! seq vec
  Vec :: sat_pf     ! mpi vec

  !
  PetscInt :: nlclm  ! num of local clm cells
  PetscInt :: ngclm  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf   ! num of local pflotran cells
  PetscInt :: ngpf   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  end type clm_pflotran_interface_data_type

  type(clm_pflotran_interface_data_type) , public, target , save :: clm_pf_idata
  
  public :: clm_pflotran_interface_data_allocate_memory
  
contains

  subroutine clm_pflotran_interface_data_allocate_memory(mycomm)
  
    implicit none
    
    PetscErrorCode :: ierr
    PetscMPIInt    :: mycomm, rank
    PetscReal      :: zero = 0.0d0
    Vec :: vec_test

    call MPI_Comm_rank(MPI_COMM_WORLD,rank, ierr)

    ! Create MPI Vectors for CLM
    call VecCreateMPI(MPI_COMM_WORLD,clm_pf_idata%nlclm,PETSC_DECIDE,clm_pf_idata%hksat_x_clm,ierr) 

    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%hksat_y_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%hksat_z_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%sucsat_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%watsat_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%bsw_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%press_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%qflx_clm,ierr)

    ! Create Seq. Vectors for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf,clm_pf_idata%hksat_x_pf,ierr)

    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%hksat_y_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%hksat_z_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%sucsat_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%watsat_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%bsw_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%press_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%qflx_pf,ierr)


    ! Create MPI Vectors for PFLOTRAN
    call VecCreateMPI(MPI_COMM_WORLD,clm_pf_idata%nlpf,PETSC_DECIDE,clm_pf_idata%sat_pf,ierr)
 
    ! Create Seq. Vectors for CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm,clm_pf_idata%sat_clm,ierr)

  end subroutine clm_pflotran_interface_data_allocate_memory

end module clm_pflotran_interface_data
