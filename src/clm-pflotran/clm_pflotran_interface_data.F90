module clm_pflotran_interface_data

  implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  private

  type, public :: clm_pflotran_idata_type

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

  ! (ii) Mesh property

  ! Area of top face
  Vec :: area_top_face_clm ! seq vec
  Vec :: area_top_face_pf  ! mpi vec

  ! Time variant data
  
  ! (i) Sink/Source of water for PFLOTRAN's 3D subsurface domain
  Vec :: qflx_clm   ! mpi vec
  Vec :: qflx_pf    ! seq vec
  
  ! (ii) Source of water and temperature of rain for PFLOTRAN's 2D surface domain
  Vec :: rain_clm   ! mpi vec
  Vec :: rain_pf    ! seq vec
  Vec :: rain_temp_clm ! seq vec
  Vec :: rain_temp_pf  ! mpi vec
  
  ! (iii) Ground heat flux BC for PFLOTRAN's subsurface domain
  !       This BC is applied on top surface of the subsurface domain
  Vec :: gflux_subsurf_clm  ! mpi vec
  Vec :: gflux_subsurf_pf   ! seq vec
  !       When running PFLOTRAN surface-subsurface simulation, ground heat flux
  !       is a SS for PFLOTRAN's surface domain.
  !
  !       Note: CLM decomposes the domain across processors in a horizontal.
  !       Thus, nlclm_2dsub = nlclm_srf across all processors. Thus, there is
  !       no need for 'gflux_surf_clm'
  !
  Vec :: gflux_surf_pf   ! seq vec

  ! (iv) Saturation
  Vec :: sat_clm    ! seq vec
  Vec :: sat_pf     ! mpi vec

  ! (v) Subsurface temperature
  Vec :: temp_clm   ! seq vec
  Vec :: temp_pf    ! mpi vec

  ! (vi) Ice saturation
  Vec :: sat_ice_clm ! seq vec
  Vec :: sat_ice_pf  ! mpi vec

  ! (vii) Stand water head
  Vec :: h2osfc_clm ! seq vec
  Vec :: h2osfc_pf  ! mpi vec

  ! Number of cells for the 3D subsurface domain
  PetscInt :: nlclm_sub ! num of local clm cells
  PetscInt :: ngclm_sub ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_sub  ! num of local pflotran cells
  PetscInt :: ngpf_sub   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the surface of the 3D subsurface domain
  PetscInt :: nlclm_2dsub  ! num of local clm cells
  PetscInt :: ngclm_2dsub  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_2dsub   ! num of local pflotran cells
  PetscInt :: ngpf_2dsub   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  ! Number of cells for the 2D surface domain
  PetscInt :: nlclm_srf  ! num of local clm cells
  PetscInt :: ngclm_srf  ! num of ghosted clm cells (ghosted = local+ghosts)
  PetscInt :: nlpf_srf   ! num of local pflotran cells
  PetscInt :: ngpf_srf   ! num of ghosted pflotran cells (ghosted = local+ghosts)

  end type clm_pflotran_idata_type

  type(clm_pflotran_idata_type) , public, target , save :: clm_pf_idata
  
  public :: CLMPFLOTRANIDataInit, &
            CLMPFLOTRANIDataCreateVec, &
            CLMPFLOTRANIDataDestroy
  
contains

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataInit()
  ! 
  ! This routine initialized the data transfer type.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  ! 
  
    implicit none

    clm_pf_idata%nlclm_sub = 0
    clm_pf_idata%ngclm_sub = 0
    clm_pf_idata%nlpf_sub = 0
    clm_pf_idata%ngpf_sub = 0

    clm_pf_idata%nlclm_2dsub = 0
    clm_pf_idata%ngclm_2dsub = 0
    clm_pf_idata%nlpf_2dsub = 0
    clm_pf_idata%ngpf_2dsub = 0

    clm_pf_idata%nlclm_srf = 0
    clm_pf_idata%ngclm_srf = 0
    clm_pf_idata%nlpf_srf = 0
    clm_pf_idata%ngpf_srf = 0

    clm_pf_idata%hksat_x_clm = 0
    clm_pf_idata%hksat_y_clm = 0
    clm_pf_idata%hksat_z_clm = 0
    clm_pf_idata%sucsat_clm = 0
    clm_pf_idata%watsat_clm = 0
    clm_pf_idata%bsw_clm = 0
    clm_pf_idata%press_clm = 0

    clm_pf_idata%hksat_x_pf = 0
    clm_pf_idata%hksat_y_pf = 0
    clm_pf_idata%hksat_z_pf = 0
    clm_pf_idata%sucsat_pf = 0
    clm_pf_idata%watsat_pf = 0
    clm_pf_idata%bsw_pf = 0
    clm_pf_idata%press_pf = 0

    clm_pf_idata%qflx_clm = 0
    clm_pf_idata%qflx_pf = 0
    
    clm_pf_idata%rain_clm = 0
    clm_pf_idata%rain_pf = 0
    clm_pf_idata%rain_temp_clm = 0
    clm_pf_idata%rain_temp_pf = 0
    
    clm_pf_idata%gflux_subsurf_clm = 0
    clm_pf_idata%gflux_subsurf_pf = 0
    clm_pf_idata%gflux_surf_pf = 0

    clm_pf_idata%sat_clm = 0
    clm_pf_idata%sat_pf = 0

    clm_pf_idata%temp_clm = 0
    clm_pf_idata%temp_pf = 0

    clm_pf_idata%sat_ice_clm = 0
    clm_pf_idata%sat_ice_pf = 0

    clm_pf_idata%h2osfc_clm = 0
    clm_pf_idata%h2osfc_pf = 0

  end subroutine CLMPFLOTRANIDataInit

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataCreateVec(mycomm)
  ! 
  ! This routine creates PETSc vectors required for data transfer between
  ! CLM and PFLOTRAN.
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 2011
  ! 
  
    implicit none
    
    PetscErrorCode :: ierr
    PetscMPIInt    :: mycomm, rank
    PetscReal      :: zero = 0.0d0
    Vec :: vec_test

    call MPI_Comm_rank(mycomm,rank, ierr)

    !
    ! For data transfer from CLM to PFLOTRAN
    !

    ! Create MPI Vectors for CLM
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_sub,PETSC_DECIDE,clm_pf_idata%hksat_x_clm,ierr)
    call VecSet(clm_pf_idata%hksat_x_clm,0.d0,ierr)

    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%hksat_y_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%hksat_z_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%sucsat_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%watsat_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%bsw_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%press_clm,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_clm,clm_pf_idata%qflx_clm,ierr)

    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_2dsub,PETSC_DECIDE,clm_pf_idata%gflux_subsurf_clm,ierr)
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_srf,PETSC_DECIDE,clm_pf_idata%rain_clm,ierr)
    call VecCreateMPI(mycomm,clm_pf_idata%nlclm_srf,PETSC_DECIDE,clm_pf_idata%rain_temp_clm,ierr)
    call VecSet(clm_pf_idata%gflux_subsurf_clm,0.d0,ierr)
    call VecSet(clm_pf_idata%rain_clm,0.d0,ierr)
    call VecSet(clm_pf_idata%rain_temp_clm,0.d0,ierr)

    ! Create Seq. Vectors for PFLOTRAN
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_sub,clm_pf_idata%hksat_x_pf,ierr)
    call VecSet(clm_pf_idata%hksat_x_pf,0.d0,ierr)

    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%hksat_y_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%hksat_z_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%sucsat_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%watsat_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%bsw_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%press_pf,ierr)
    call VecDuplicate(clm_pf_idata%hksat_x_pf,clm_pf_idata%qflx_pf,ierr)

    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_2dsub,clm_pf_idata%gflux_subsurf_pf,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_srf,clm_pf_idata%gflux_surf_pf,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_srf,clm_pf_idata%rain_pf,ierr)
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngpf_srf,clm_pf_idata%rain_temp_pf,ierr)
    call VecSet(clm_pf_idata%gflux_subsurf_pf,0.d0,ierr)
    call VecSet(clm_pf_idata%gflux_surf_pf,0.d0,ierr)
    call VecSet(clm_pf_idata%rain_pf,0.d0,ierr)
    call VecSet(clm_pf_idata%rain_temp_pf,0.d0,ierr)

    !
    ! For data transfer from PFLOTRAN to CLM
    !

    ! Create MPI Vectors for PFLOTRAN
    ! 3D Subsurface PFLOTRAN ---to--- 3D Subsurface CLM
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_sub,PETSC_DECIDE,clm_pf_idata%sat_pf,ierr)
    call VecSet(clm_pf_idata%sat_pf,0.d0,ierr)

    call VecDuplicate(clm_pf_idata%sat_pf,clm_pf_idata%temp_pf,ierr)
    call VecDuplicate(clm_pf_idata%sat_pf,clm_pf_idata%sat_ice_pf,ierr)
    call VecDuplicate(clm_pf_idata%sat_pf,clm_pf_idata%area_top_face_pf,ierr)

    ! 2D Surface PFLOTRAN ---to--- 2D Surface CLM
    call VecCreateMPI(mycomm,clm_pf_idata%nlpf_srf,PETSC_DECIDE,clm_pf_idata%h2osfc_pf,ierr)
    call VecSet(clm_pf_idata%h2osfc_pf,0.d0,ierr)

    ! Create Seq. Vectors for CLM
    ! 3D Subsurface PFLOTRAN ---to--- 3D Subsurface CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%ngclm_sub,clm_pf_idata%sat_clm,ierr)
    call VecSet(clm_pf_idata%sat_clm,0.d0,ierr)

    call VecDuplicate(clm_pf_idata%sat_clm,clm_pf_idata%temp_clm,ierr)
    call VecDuplicate(clm_pf_idata%sat_clm,clm_pf_idata%sat_ice_clm,ierr)
    call VecDuplicate(clm_pf_idata%sat_clm,clm_pf_idata%area_top_face_clm,ierr)

    ! 2D Surface PFLOTRAN ---to--- 2D Surface CLM
    call VecCreateSeq(PETSC_COMM_SELF,clm_pf_idata%nlclm_2dsub,clm_pf_idata%h2osfc_clm,ierr)
    call VecSet(clm_pf_idata%h2osfc_clm,0.d0,ierr)

  end subroutine CLMPFLOTRANIDataCreateVec

! ************************************************************************** !

  subroutine CLMPFLOTRANIDataDestroy()
  ! 
  ! This routine destroys PETSc vectors that were created for data transfer.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 4/10/2013
  ! 
  
    implicit none
    
    PetscErrorCode :: ierr

    if(clm_pf_idata%hksat_x_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_x_clm,ierr)
    if(clm_pf_idata%hksat_y_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_y_clm,ierr)
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)
    if(clm_pf_idata%hksat_z_clm  /= 0) call VecDestroy(clm_pf_idata%hksat_z_clm,ierr)

    if(clm_pf_idata%hksat_x_pf  /= 0) call VecDestroy(clm_pf_idata%hksat_x_pf,ierr)
    if(clm_pf_idata%hksat_y_pf  /= 0) call VecDestroy(clm_pf_idata%hksat_y_pf,ierr)
    if(clm_pf_idata%hksat_z_pf  /= 0) call VecDestroy(clm_pf_idata%hksat_z_pf,ierr)
    if(clm_pf_idata%sucsat_pf  /= 0) call VecDestroy(clm_pf_idata%sucsat_pf,ierr)
    if(clm_pf_idata%watsat_pf  /= 0) call VecDestroy(clm_pf_idata%watsat_pf,ierr)
    if(clm_pf_idata%bsw_pf  /= 0) call VecDestroy(clm_pf_idata%bsw_pf,ierr)
    if(clm_pf_idata%press_pf  /= 0) call VecDestroy(clm_pf_idata%press_pf,ierr)

    if(clm_pf_idata%qflx_clm  /= 0) call VecDestroy(clm_pf_idata%qflx_clm,ierr)
    if(clm_pf_idata%qflx_pf  /= 0) call VecDestroy(clm_pf_idata%qflx_pf,ierr)
    
    if(clm_pf_idata%rain_clm  /= 0) call VecDestroy(clm_pf_idata%rain_clm,ierr)
    if(clm_pf_idata%rain_pf  /= 0) call VecDestroy(clm_pf_idata%rain_pf,ierr)
    if(clm_pf_idata%rain_temp_clm  /= 0) call VecDestroy(clm_pf_idata%rain_temp_clm,ierr)
    if(clm_pf_idata%rain_temp_pf  /= 0) call VecDestroy(clm_pf_idata%rain_temp_pf,ierr)
    
    if(clm_pf_idata%gflux_subsurf_clm  /= 0) call VecDestroy(clm_pf_idata%gflux_subsurf_clm,ierr)
    if(clm_pf_idata%gflux_subsurf_pf  /= 0) call VecDestroy(clm_pf_idata%gflux_subsurf_pf,ierr)
    if(clm_pf_idata%gflux_surf_pf  /= 0) call VecDestroy(clm_pf_idata%gflux_surf_pf,ierr)

    if(clm_pf_idata%sat_clm  /= 0) call VecDestroy(clm_pf_idata%sat_clm,ierr)
    if(clm_pf_idata%sat_pf  /= 0) call VecDestroy(clm_pf_idata%sat_pf,ierr)

    if(clm_pf_idata%temp_clm  /= 0) call VecDestroy(clm_pf_idata%temp_clm,ierr)
    if(clm_pf_idata%temp_pf  /= 0) call VecDestroy(clm_pf_idata%temp_pf,ierr)

    if(clm_pf_idata%sat_ice_clm  /= 0) call VecDestroy(clm_pf_idata%sat_ice_clm,ierr)
    if(clm_pf_idata%sat_ice_pf  /= 0) call VecDestroy(clm_pf_idata%sat_ice_pf,ierr)

    if(clm_pf_idata%h2osfc_clm  /= 0) call VecDestroy(clm_pf_idata%h2osfc_clm,ierr)
    if(clm_pf_idata%h2osfc_pf  /= 0) call VecDestroy(clm_pf_idata%h2osfc_pf,ierr)

    if(clm_pf_idata%area_top_face_clm  /= 0) &
      call VecDestroy(clm_pf_idata%area_top_face_clm,ierr)
    if(clm_pf_idata%area_top_face_pf  /= 0) &
      call VecDestroy(clm_pf_idata%area_top_face_pf,ierr)

  end subroutine CLMPFLOTRANIDataDestroy

end module clm_pflotran_interface_data
