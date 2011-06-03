program pflotran_interface_main

  use pflotran_model_module
  use clm_pflotran_interface_data
  use Mapping_module
  use Input_module
  
  implicit none

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
  character(len=MAXSTRINGLENGTH)     :: filename
  character(len=MAXWORDLENGTH)       :: card
  PetscViewer                        :: viewer
  PetscInt, pointer                  :: tmp_id(:)
  PetscInt                           :: npts
  
  
  Vec :: sat_clmloc
  Vec :: sat_pfloc


  allocate(pflotran_m)
  allocate(clmpf_data)

  ! Create the model and CLM data
  pflotran_m => pflotranModelCreate()
  clmpf_data => clm_pf_data_create()


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
  call PetscViewerASCIIOpen(pflotran_m%option%mycomm, 'sat_clmloc.out',viewer,ierr)
  call VecView(sat_clmloc, viewer,ierr)
  call PetscViewerDestroy(viewer, ierr)
  

  call VecCreateMPI(pflotran_m%option%mycomm, pflotran_m%map_clm2pf%d_ncells_ghosted, &
       PETSC_DECIDE, sat_pfloc, ierr)

  !write(*,*), 's_ncells_local = ',pflotran_m%map_clm2pf%s_ncells_local

  call MappingSourceToDestination( pflotran_m%map_clm2pf, pflotran_m%option, sat_clmloc, &
       sat_pfloc)
#endif

  !call pflotranModelStepperRunInit(pflotran_m)
  !do time = 1,0
  !   call pflotranModelStepperRunTillPauseTime(pflotran_m,time * 1800.0d0)
  !enddo
  !call pflotranModelStepperRunFinalize(pflotran_m)

  !call pflotranModelDestroy(pflotran_m)

end program pflotran_interface_main
