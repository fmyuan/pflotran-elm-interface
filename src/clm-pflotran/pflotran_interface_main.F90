program pflotran_interface_main

  use pflotran_model_module
  use clm_pflotran_interface_data
  use Mapping_module
  use Input_module
  use HDF5_aux_module
  use hdf5
  use Grid_module
  use Field_module
  use Richards_Aux_module
  use Richards_module
  use Discretization_module
  
  implicit none

  include "piof.h"
#include "finclude/petsclog.h"
#include "finclude/petscsysdef.h"
#include "finclude/petscviewer.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

  PetscInt :: time

  type(pflotran_model_type),pointer  :: pflotran_m
  type(clm_pflotran_data),pointer    :: clmpf_data
  PetscInt, pointer                  :: clm_cell_ids(:)
  PetscInt                           :: clm_npts, ii,fileid,num_u_a,jj
  PetscReal, pointer                 :: v_loc(:),u_a(:)
  type(input_type), pointer          :: input
  PetscViewer                        :: viewer
  PetscInt, pointer                  :: tmp_id(:)
  PetscInt                           :: npts
  Vec :: sat_clmloc
  Vec :: sat_pfloc

  ! To read HDF5 soil properties  
  character(len=MAXSTRINGLENGTH)     :: filename
  character(len=MAXSTRINGLENGTH)     :: group_name
  character(len=MAXSTRINGLENGTH)     :: dataset_name
  character(len=MAXWORDLENGTH)       :: card
  character(len=MAXSTRINGLENGTH)     :: alp_dname = '/Material/alpha'//CHAR(0)
  character(len=MAXSTRINGLENGTH)     :: lam_dname = '/Material/lambda'//CHAR(0)
  character(len=MAXSTRINGLENGTH)     :: por_dname = '/Material/porosity'//CHAR(0)
  character(len=MAXSTRINGLENGTH)     :: pxx_dname = '/Material/perm_x'//CHAR(0)
  character(len=MAXSTRINGLENGTH)     :: pyy_dname = '/Material/perm_y'//CHAR(0)
  character(len=MAXSTRINGLENGTH)     :: pzz_dname = '/Material/perm_z'//CHAR(0)
  character(len=MAXSTRINGLENGTH)     :: rch_dname = '/Recharge'//CHAR(0)
  character(len=MAXSTRINGLENGTH)     :: ics_dname = '/Pressure'//CHAR(0)
  
  !
  PetscReal,pointer :: alp_p(:)
  PetscReal,pointer :: lam_p(:)
  PetscReal,pointer :: por_p(:)
  PetscReal,pointer :: pxx_p(:)
  PetscReal,pointer :: pyy_p(:)
  PetscReal,pointer :: pzz_p(:)
  PetscInt          :: data_dims(1)
  PetscInt          :: dataset_dims(1)
  PetscReal,pointer :: rch_2d_p(:,:), rch_1d_p(:)
  PetscInt          :: rch_data_dims(2), rch_dataset_dims(2)
  PetscReal,pointer :: ics_p(:)
  PetscInt          :: ics_data_dims(1), ics_dataset_dims(1)

  !
  Vec :: lam_nat_v,lam_loc_v
  Vec :: por_nat_v,por_loc_v
  Vec :: alp_nat_v,alp_loc_v
  Vec :: pxx_nat_v,pxx_loc_v
  Vec :: pyy_nat_v,pyy_loc_v
  Vec :: pzz_nat_v,pzz_loc_v
  Vec :: rch_nat_v,rch_loc_v
  Vec :: ics_nat_v,ics_loc_v
  Vec :: por_pfloc

  !
  PetscScalar,pointer :: v_loc_1(:),v_loc_2(:),v_loc_3(:),v_loc_4(:)
  PetscScalar,pointer :: v_loc_5(:),v_loc_6(:),v_loc_7(:),v_loc_8(:)
  PetscInt, pointer   :: tmp_int_array(:)
  PetscReal,pointer   :: tmp_real_array(:)
  IS                  :: is_from, is_to
  VecScatter          :: vec_scat
  PetscInt            :: istart,iend
  
  type(grid_type),pointer              :: grid
  type(field_type),pointer             :: field
  type(richards_auxvar_type), pointer  :: rich_aux_vars(:)
  type(richards_auxvar_type), pointer  :: aux_var

  PetscInt :: PRINT_RANK    
  PRINT_RANK = 0

  ! A
  allocate(pflotran_m)
  allocate(clmpf_data)    

  ! Create the model and CLM data
  pflotran_m => pflotranModelCreate()
  clmpf_data => clm_pf_data_create()

  clm_npts = 2*2*10/pflotran_m%option%mycommsize
  allocate(clm_cell_ids(clm_npts))
  do ii = 1,clm_npts
    clm_cell_ids(ii) = ii-1 + clm_npts*pflotran_m%option%myrank
    if(pflotran_m%option%myrank.eq.PRINT_RANK) write(*,*), ii, clm_cell_ids(ii)
  enddo

  filename = 'idealize_mapping_clm2pf_flux_2x2_clm_matching.dat'//CHAR(0)
  filename = 'idealize_mapping_clm2pf_flux_2x2.dat'//CHAR(0)
  call pflotranModelInitMapping3(pflotran_m, &
    clm_cell_ids,clm_npts,1,1)

  filename = 'idealize_mapping_clm2pf_soil_2x2.dat'//CHAR(0)
  call pflotranModelInitMapping3(pflotran_m, &
    clm_cell_ids,clm_npts,2,1)

  !if(pflotran_m%option%myrank.eq.PRINT_RANK) write(*,*), 'pflotranModelStepperRunInit'
  !call pflotranModelStepperRunInit(pflotran_m)
  !do time = 1,0
  !   call pflotranModelStepperRunTillPauseTime(pflotran_m,time * 3600.0d0)
  !enddo
  !call pflotranModelStepperRunFinalize(pflotran_m)

  ! Read soil properties data  
  filename = 'soil_prop_2x2pt.h5'
  call HDF5ReadDatasetReal1D(filename,alp_dname,NONUNIFORM_CONTIGUOUS_READ,&
    pflotran_m%option,alp_p,data_dims,dataset_dims)
  !write(*,*),'alpha - data_dims: ',data_dims(:)

  call HDF5ReadDatasetReal1D(filename,lam_dname,NONUNIFORM_CONTIGUOUS_READ,&
    pflotran_m%option,lam_p,data_dims,dataset_dims)
  !write(*,*),'lambda - data_dims: ',data_dims(:)

  call HDF5ReadDatasetReal1D(filename,por_dname,NONUNIFORM_CONTIGUOUS_READ,&
    pflotran_m%option,por_p,data_dims,dataset_dims)
  !write(*,*),'porosity - data_dims: ',data_dims(:)

  call HDF5ReadDatasetReal1D(filename,pxx_dname,NONUNIFORM_CONTIGUOUS_READ,&
    pflotran_m%option,pxx_p,data_dims,dataset_dims)
  !write(*,*),'px data_dims: ',data_dims(:)

  call HDF5ReadDatasetReal1D(filename,pyy_dname,NONUNIFORM_CONTIGUOUS_READ,&
    pflotran_m%option,pyy_p,data_dims,dataset_dims)
  !write(*,*),'py data_dims: ',data_dims(:)

  call HDF5ReadDatasetReal1D(filename,pzz_dname,NONUNIFORM_CONTIGUOUS_READ,&
    pflotran_m%option,pzz_p,data_dims,dataset_dims)
  !write(*,*),'pz data_dims: ',data_dims(:)
  
  ! Read initial conditions
  !filename = 'init_cond_conus_10min_mesh.h5'
  !call HDF5ReadDatasetReal1D(filename,ics_dname,NONUNIFORM_CONTIGUOUS_READ,&
  !  pflotran_m%option,ics_p,data_dims,dataset_dims)
  !write(*,*),'ics_data_dims: ',data_dims(:),dataset_dims(:)

  ! Create vectors to save soil properties
  call VecCreateMPI(pflotran_m%option%mycomm, PETSC_DECIDE, dataset_dims(1), alp_nat_v, ierr)
  call VecCreateMPI(pflotran_m%option%mycomm, PETSC_DECIDE, dataset_dims(1), lam_nat_v, ierr)
  call VecCreateMPI(pflotran_m%option%mycomm, PETSC_DECIDE, dataset_dims(1), por_nat_v, ierr)
  call VecCreateMPI(pflotran_m%option%mycomm, PETSC_DECIDE, dataset_dims(1), pxx_nat_v, ierr)
  call VecCreateMPI(pflotran_m%option%mycomm, PETSC_DECIDE, dataset_dims(1), pyy_nat_v, ierr)
  call VecCreateMPI(pflotran_m%option%mycomm, PETSC_DECIDE, dataset_dims(1), pzz_nat_v, ierr)
  
  ! Create vector to save initial conditions
  !call VecCreateMPI(pflotran_m%option%mycomm, PETSC_DECIDE, dataset_dims(1), ics_nat_v, ierr)
  
  ! Save the data into vectors in natural index  
  call VecGetArrayF90(alp_nat_v,v_loc_1,ierr)
  call VecGetArrayF90(lam_nat_v,v_loc_2,ierr)
  call VecGetArrayF90(por_nat_v,v_loc_3,ierr)
  call VecGetArrayF90(pxx_nat_v,v_loc_4,ierr)
  call VecGetArrayF90(pyy_nat_v,v_loc_5,ierr)
  call VecGetArrayF90(pzz_nat_v,v_loc_6,ierr)
  !call VecGetArrayF90(ics_nat_v,v_loc_7,ierr)
  
  do ii = 1,data_dims(1)
    v_loc_1(ii) = alp_p(ii)
    v_loc_2(ii) = lam_p(ii)
    v_loc_3(ii) = por_p(ii)
    v_loc_4(ii) = pxx_p(ii)
    v_loc_5(ii) = pyy_p(ii)
    v_loc_6(ii) = pzz_p(ii)
    !v_loc_7(ii) = ics_p(ii)
  enddo
      
  call VecRestoreArrayF90(alp_nat_v,v_loc_1,ierr)
  call VecRestoreArrayF90(lam_nat_v,v_loc_2,ierr)
  call VecRestoreArrayF90(por_nat_v,v_loc_3,ierr)
  call VecRestoreArrayF90(pxx_nat_v,v_loc_4,ierr)
  call VecRestoreArrayF90(pyy_nat_v,v_loc_5,ierr)
  call VecRestoreArrayF90(pzz_nat_v,v_loc_6,ierr)
  !call VecRestoreArrayF90(ics_nat_v,v_loc_7,ierr)

  ! Free memory
  deallocate(alp_p)
  deallocate(lam_p)
  deallocate(por_p)
  deallocate(pxx_p)
  deallocate(pyy_p)
  deallocate(pzz_p)
  !deallocate(ics_p)
  
  call VecCreateMPI(pflotran_m%option%mycomm, &
                    pflotran_m%map_clm2pf%d_ncells_ghosted, &
                    PETSC_DECIDE, por_pfloc, ierr)

  call MappingSourceToDestination( pflotran_m%map_clm2pf_soils, pflotran_m%option, por_nat_v, &
       por_pfloc)

  call PetscViewerASCIIOpen(pflotran_m%option%mycomm, 'por_pfloc.out', viewer, ierr)
  call VecView(por_pfloc, viewer,ierr)
  call PetscViewerDestroy(viewer, ierr)

  call PetscViewerASCIIOpen(pflotran_m%option%mycomm, 'por.out', viewer, ierr)
  call VecView(por_nat_v, viewer,ierr)
  call PetscViewerDestroy(viewer, ierr)

  if(pflotran_m%option%myrank.eq.PRINT_RANK) write(*,*), 'pflotranModelDestroy'
  call pflotranModelDestroy(pflotran_m)



#if 0 
  grid      => pflotran_m%realization%patch%grid
  field     => pflotran_m%realization%field
  
  ! Create index set - Scattering from global vector
  allocate(tmp_int_array(grid%ngmax))
  do ii = 1,grid%ngmax
    tmp_int_array(ii) = grid%nG2A(ii)
  enddo

  ! Create vectors to save soil properties and initial conditions corresponding 
  ! to ghosted cells present on each proc
  call VecCreateMPI(pflotran_m%option%mycomm, grid%ngmax, PETSC_DECIDE, alp_loc_v, ierr)
  call VecCreateMPI(pflotran_m%option%mycomm, grid%ngmax, PETSC_DECIDE, lam_loc_v, ierr)
  call VecCreateMPI(pflotran_m%option%mycomm, grid%ngmax, PETSC_DECIDE, por_loc_v, ierr)
  call VecCreateMPI(pflotran_m%option%mycomm, grid%ngmax, PETSC_DECIDE, pxx_loc_v, ierr)
  call VecCreateMPI(pflotran_m%option%mycomm, grid%ngmax, PETSC_DECIDE, pyy_loc_v, ierr)
  call VecCreateMPI(pflotran_m%option%mycomm, grid%ngmax, PETSC_DECIDE, pzz_loc_v, ierr)

  tmp_int_array = tmp_int_array - 1
  call ISCreateBlock(pflotran_m%option%mycomm, 1, grid%ngmax, tmp_int_array, PETSC_COPY_VALUES, &
         is_from, ierr)
  deallocate(tmp_int_array)
  !call PetscViewerASCIIOpen(pflotran_m%option%mycomm, 'is_from.out', viewer, ierr)
  !call ISView(is_from, viewer,ierr)
  !call PetscViewerDestroy(viewer, ierr)
  
  ! Create index set - Scattering to 
  call VecGetOwnershipRange(alp_loc_v,istart,iend,ierr)
  allocate(tmp_int_array(grid%ngmax))
  do ii = 1,grid%ngmax
    tmp_int_array(ii) = ii-1+istart
  enddo
  call ISCreateBlock(pflotran_m%option%mycomm, 1, grid%ngmax, tmp_int_array, PETSC_COPY_VALUES, &
         is_to, ierr)
  !call PetscViewerASCIIOpen(pflotran_m%option%mycomm, 'is_to.out', viewer, ierr)
  !call ISView(is_to, viewer,ierr)
  !call PetscViewerDestroy(viewer, ierr)

  ! Create vector scatter
  call VecScatterCreate(alp_nat_v, is_from, alp_loc_v, is_to, vec_scat, ierr)
  !call PetscViewerASCIIOpen(pflotran_m%option%mycomm, 'vec_scat.out', viewer, ierr)
  !call VecScatterView(vec_scat, viewer,ierr)
  !call PetscViewerDestroy(viewer, ierr)
  !call ISDestroy(is_from)
  !call ISDestroy(is_to)
  
  ! Scatter vectors
  call VecScatterBegin(vec_scat, alp_nat_v, alp_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(  vec_scat, alp_nat_v, alp_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)

  call VecScatterBegin(vec_scat, lam_nat_v, lam_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(  vec_scat, lam_nat_v, lam_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)
  
  call VecScatterBegin(vec_scat, por_nat_v, por_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(  vec_scat, por_nat_v, por_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)

  call VecScatterBegin(vec_scat, pxx_nat_v, pxx_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(  vec_scat, pxx_nat_v, pxx_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)

  call VecScatterBegin(vec_scat, pyy_nat_v, pyy_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(  vec_scat, pyy_nat_v, pyy_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)

  call VecScatterBegin(vec_scat, pzz_nat_v, pzz_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(  vec_scat, pzz_nat_v, pzz_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)
  
  call VecScatterDestroy(vec_scat,ierr)
#if 0
  ! ========================================================================
  ! Save the porosity and permeability data
  ! ========================================================================
  call GridVecGetArrayF90(grid,field%porosity_loc, v_loc_1, ierr)
  call GridVecGetArrayF90(grid,field%perm_xx_loc,  v_loc_2, ierr)
  call GridVecGetArrayF90(grid,field%perm_yy_loc,  v_loc_3, ierr)
  call GridVecGetArrayF90(grid,field%perm_zz_loc,  v_loc_4, ierr)

  call VecGetArrayF90(por_loc_v, v_loc_5, ierr)
  call VecGetArrayF90(pxx_loc_v, v_loc_6, ierr)
  call VecGetArrayF90(pyy_loc_v, v_loc_7, ierr)
  call VecGetArrayF90(pzz_loc_v, v_loc_8, ierr)

  do ii = 1,grid%ngmax
    !if(ii.eq.1) write (*,*), v_loc_1(ii),v_loc_2(ii),v_loc_3(ii),v_loc_4(ii)
    !if(ii.eq.1) write (*,*), v_loc_5(ii),v_loc_6(ii),v_loc_7(ii),v_loc_8(ii)
    v_loc_1(ii) = v_loc_5(ii)
    v_loc_2(ii) = v_loc_6(ii)
    v_loc_3(ii) = v_loc_7(ii)
    v_loc_4(ii) = v_loc_8(ii)
    !if(ii.eq.1) write (*,*), v_loc_1(ii),v_loc_2(ii),v_loc_3(ii),v_loc_4(ii)
  enddo
  
  call GridVecRestoreArrayF90(grid,field%porosity_loc, v_loc_1, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_xx_loc,  v_loc_2, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_yy_loc,  v_loc_3, ierr)
  call GridVecRestoreArrayF90(grid,field%perm_zz_loc,  v_loc_4, ierr)

  call VecRestoreArrayF90(por_loc_v, v_loc_5, ierr)
  call VecRestoreArrayF90(pxx_loc_v, v_loc_6, ierr)
  call VecRestoreArrayF90(pyy_loc_v, v_loc_7, ierr)
  call VecRestoreArrayF90(pzz_loc_v, v_loc_8, ierr)

  ! ========================================================================
  ! Save lambda and alpha
  ! ========================================================================
  rich_aux_vars   => pflotran_m%realization%patch%aux%Richards%aux_vars
  call VecGetArrayF90(alp_loc_v, v_loc_1, ierr)
  call VecGetArrayF90(lam_loc_v, v_loc_2, ierr)

  do ii = 1,grid%ngmax
    aux_var => rich_aux_vars(ii)
    aux_var%bc_alpha  = v_loc_1(ii)
    aux_var%bc_lambda = v_loc_2(ii)
  enddo
  
  call VecRestoreArrayF90(alp_loc_v, v_loc_1, ierr)
  call VecRestoreArrayF90(lam_loc_v, v_loc_2, ierr)
  
  ! ========================================================================
  ! Save initial conditions
  ! ========================================================================
  call VecCreateMPI(pflotran_m%option%mycomm, grid%nlmax, PETSC_DECIDE, ics_loc_v, ierr)

  ! Create index set - Scattering from global vector
  allocate(tmp_int_array(grid%ngmax))
  do ii = 1,grid%nlmax
    tmp_int_array(ii) = grid%nL2A(ii)
  enddo

  !tmp_int_array = tmp_int_array - 1
  call ISCreateBlock(pflotran_m%option%mycomm, 1, grid%nlmax, tmp_int_array, PETSC_COPY_VALUES, &
         is_from, ierr)
  deallocate(tmp_int_array)
  !call PetscViewerASCIIOpen(pflotran_m%option%mycomm, 'is_fromm.out', viewer, ierr)
  !call ISView(is_from, viewer,ierr)
  !call PetscViewerDestroy(viewer, ierr)
  
  ! Create index set - Scattering to 
  call VecGetOwnershipRange(ics_loc_v,istart,iend,ierr)
  allocate(tmp_int_array(grid%nlmax))
  do ii = 1,grid%nlmax
    tmp_int_array(ii) = ii-1+istart
  enddo
  call ISCreateBlock(pflotran_m%option%mycomm, 1, grid%nlmax, tmp_int_array, PETSC_COPY_VALUES, &
         is_to, ierr)
  !call PetscViewerASCIIOpen(pflotran_m%option%mycomm, 'is_too.out', viewer, ierr)
  !call ISView(is_to, viewer,ierr)
  !call PetscViewerDestroy(viewer, ierr)

  ! Create vector scatter
  call VecScatterCreate(ics_nat_v, is_from, ics_loc_v, is_to, vec_scat, ierr)

  call VecScatterBegin(vec_scat, ics_nat_v, ics_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(  vec_scat, ics_nat_v, ics_loc_v, INSERT_VALUES, SCATTER_FORWARD, ierr)

  call VecScatterDestroy(vec_scat, ierr)

  call GridVecGetArrayF90(grid,field%flow_xx,tmp_real_array,ierr)
  call VecGetArrayF90(ics_loc_v,v_loc_2,ierr)

  do ii = 1, grid%nlmax
       tmp_real_array(ii) = v_loc_2(ii)
       !tmp_real_array(ii) = REAL(v_loc_2(ii))
  enddo

  call GridVecRestoreArrayF90(grid,field%flow_xx,tmp_real_array,ierr)
  call VecRestoreArrayF90(ics_loc_v,v_loc_2,ierr)

  call DiscretizationGlobalToLocal(pflotran_m%realization%discretization,&
    field%flow_xx,field%flow_xx_loc,NFLOWDOF)
  call VecCopy(field%flow_xx, field%flow_yy, ierr)
  call RichardsUpdateAuxVars(pflotran_m%realization)
#endif
  ! ========================================================================
  !                             Read forcing data
  ! ========================================================================
  !filename = '01_recharge.h5'
  !call HDF5ReadDatasetReal2D(filename,rch_dname,NONUNIFORM_CONTIGUOUS_READ,&
  !  pflotran_m%option,rch_2d_p,rch_data_dims,rch_dataset_dims)
  !write(*,*),'rch_data_dims: ',rch_data_dims(:),rch_dataset_dims(:)

  ! ========================================================================
  !                             Read forcing data
  ! ========================================================================

  ! ========================================================================
  !                             Mapping
  ! ========================================================================
#if 0
  clm_npts = 57*98*10/pflotran_m%option%mycommsize
  allocate(clm_cell_ids(clm_npts))
  do ii = 1,clm_npts
    clm_cell_ids(ii) = ii-1 + clm_npts*pflotran_m%option%myrank
  enddo

  filename = 'conus_10min_from_clm_subset_wts_matrix.txt'//CHAR(0)
  call pflotranModelInitMapping3(pflotran_m, &
    clm_cell_ids,clm_npts,1,1)
#endif


  call pflotranModelStepperRunInit(pflotran_m)
  do time = 1,10
     call pflotranModelStepperRunTillPauseTime(pflotran_m,time * 3600.0d0)
  enddo
  call pflotranModelStepperRunFinalize(pflotran_m)

  call pflotranModelDestroy(pflotran_m)


#if 0

#if 1
  !clm_npts = 8192/pflotran_m%option%mycommsize
  clm_npts = 5586/pflotran_m%option%mycommsize
  allocate(clm_cell_ids(clm_npts))
  do ii = 1,clm_npts
    clm_cell_ids(ii) = ii-1 + clm_npts*pflotran_m%option%myrank
  enddo

#else

  filename = 'proc_id_4.txt'
  fileid   = 10
  input => InputCreate(fileid,filename)
  call InputReadFlotranString(input,pflotran_m%option)
  call InputReadInt(input,pflotran_m%option,npts)
  write(*,*), 'npts = ', npts

  allocate(tmp_id(npts))
  
  clm_npts = 0
  do ii = 1,npts
    call InputReadFlotranString(input,pflotran_m%option)
	call InputReadInt(input,pflotran_m%option,tmp_id(ii))
	if (pflotran_m%option%myrank.eq.tmp_id(ii)) clm_npts = clm_npts + 1
  enddo
  write(*,*), 'clm_npts = ',pflotran_m%option%myrank,clm_npts
  allocate(clm_cell_ids(clm_npts))
  
  clm_npts = 0
  do ii = 1,npts
    if (pflotran_m%option%myrank.eq.tmp_id(ii)) then
	  clm_npts = clm_npts + 1
	  clm_cell_ids(clm_npts) = ii-1
	endif
  enddo

  deallocate(tmp_id)

#endif

  !call pflotranModelInitMapping2(pflotran_m,clmpf_data,clm_cell_ids,clm_npts)
  !call pflotranModelInitMapping3(pflotran_m,clmpf_data,clm_cell_ids,clm_npts)
  filename = 'conus_10min_smallmesh_from_clm_subset_wts_matrix.txt'//CHAR(0)
  call pflotranModelInitMapping3(pflotran_m,filename,&
    clm_cell_ids,clm_npts,1,1)

#if 0
  filename = 'u_a.txt'
  fileid   = 10
  input => InputCreate(fileid,filename)
  call InputReadFlotranString(input,pflotran_m%option)
  call InputErrorMsg(input,pflotran_m%option,'number of cells',card)

  call InputReadInt(input,pflotran_m%option,num_u_a)

  allocate(u_a(num_u_a))
  
  do ii = 1,num_u_a
    call InputReadFlotranString(input,pflotran_m%option)
	call InputReadInt(input,pflotran_m%option,jj)
    call InputReadDouble(input,pflotran_m%option,u_a(ii))
  enddo

  call VecCreateMPI(pflotran_m%option%mycomm, pflotran_m%map_clm2pf%s_ncells_local, &
       PETSC_DECIDE, sat_clmloc, ierr)
  call VecGetArrayF90(sat_clmloc,v_loc,ierr)
    
  do ii = 1,pflotran_m%map_clm2pf%s_ncells_local
    !v_loc(ii) = u_a(ii + pflotran_m%option%myrank*pflotran_m%map_clm2pf%s_ncells_local)
	v_loc(ii) = u_a(clm_cell_ids(ii)+1)
  enddo 
  call VecRestoreArrayF90(sat_clmloc,v_loc,ierr)
  !call VecView(sat_clmloc,PETSC_VIEWER_STDOUT_WORLD)
  !call PetscViewerASCIIOpen(pflotran_m%option%mycomm, 'sat_clmloc.out',viewer,ierr)
  !call VecView(sat_clmloc, viewer,ierr)
  !call PetscViewerDestroy(viewer, ierr)
  

  call VecCreateMPI(pflotran_m%option%mycomm, pflotran_m%map_clm2pf%d_ncells_ghosted, &
       PETSC_DECIDE, sat_pfloc, ierr)

  !write(*,*), 's_ncells_local = ',pflotran_m%map_clm2pf%s_ncells_local

  call MappingSourceToDestination( pflotran_m%map_clm2pf, pflotran_m%option, sat_clmloc, &
       sat_pfloc)
#endif

#endif ! #if 0 @ line no 64

#endif !

end program pflotran_interface_main
