module Output_MassBalance_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Logging_module 
  use Output_Aux_module
  use Output_Common_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  ! flags signifying the first time a routine is called during a given
  ! simulation
  PetscBool :: mass_balance_first

  public :: OutputMassBalance
            
contains

! ************************************************************************** !

subroutine OutputMassBlanceInit(num_steps)
  ! 
  ! Initializes module variables for observation
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/16/13
  ! 

  use Option_module

  implicit none
  
  PetscInt :: num_steps
  
  if (num_steps == 0) then
    mass_balance_first = PETSC_TRUE
  else
    mass_balance_first = PETSC_FALSE
  endif

end subroutine OutputMassBlanceInit

! ************************************************************************** !

subroutine OutputMassBalance(realization_base)
  ! 
  ! Print to Tecplot POINT format, for MassBalance of all forms of species by species
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/18/08
  ! 

  use Realization_Subsurface_class, only : realization_subsurface_type
  use Realization_Base_class, only : realization_base_type
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Utility_module
  use Output_Aux_module
  
  use MpFlow_module

  use MpFlow_Aux_module
  use Material_Aux_class
  use String_module

  implicit none

  class(realization_base_type), target :: realization_base

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(output_option_type), pointer :: output_option
  type(coupler_type), pointer :: coupler
  type(mass_balance_region_type), pointer :: cur_mbr
  type(mpflow_auxvar_type), pointer :: mpflow_auxvars_bc_or_ss(:)

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: units
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: fid = 86
  PetscInt :: icol
  PetscInt :: iconn
  PetscInt :: offset

  PetscInt :: iphase, ispec
  PetscReal,pointer :: mass_balance(:,:)
  PetscReal :: sum_local(realization_base%option%flow%nspecflow)
  PetscReal :: sum_global(realization_base%option%flow%nspecflow)
  
  PetscReal :: global_water_mass
  PetscReal :: sum_allwater  ! sum of global water mass and fluxes
  
  PetscMPIInt :: int_mpi
  PetscBool :: bcs_done
  PetscErrorCode :: ierr

  !-----------------------------------------------------------------------------

  patch => realization_base%patch
  grid => patch%grid
  option => realization_base%option

  allocate(mass_balance(option%flow%nfluid+1, &
                        option%flow%nspecflow))

  output_option => realization_base%output_option

  if (option%flow%nspecflow<=0) return

  do ispec = 1, option%flow%nspecflow  ! not yet output those for energy (TODO)

    if (len_trim(output_option%plot_name) > 2) then
      filename = trim(output_option%plot_name) // &
               '_' // trim(StringFormatInt(ispec)) // &
               '-mas.dat'
    else
      filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '_' // trim(StringFormatInt(ispec)) // &
               '-mas.dat'
    endif

    !-------------------------
    ! open file and write headers if first time
    if (option%myrank == option%io_rank) then

      if (output_option%print_column_ids) then
        icol = 1
      else
        icol = -1
      endif
  
      if (mass_balance_first .or. .not.FileExists(filename)) then
        open(unit=fid,file=filename,action="write",status="replace")

        ! write header
        write(fid,'(a)',advance="no") ' "Time [' // trim(output_option%tunit) // &
        ']"'  
      
        if (option%iflowmode > 0) then
          call OutputWriteToHeader(fid,'dt_flow',output_option%tunit,'',icol)
        endif
      
        ! total mass
        select case(option%iflowmode)
          case(MPFLOW_MODE)
            call OutputWriteToHeader(fid,'Global: Water(L)', &
                                    '-kg','',icol)
            call OutputWriteToHeader(fid,'Global: Water(G)', &
                                    '-kg','',icol)
            call OutputWriteToHeader(fid,'Global: Water(S)', &
                                    '-kg','',icol)
          case default
            ! do nothing
        end select

        ! BC/SS
        coupler => patch%boundary_condition_list%first
        bcs_done = PETSC_FALSE
        do
          if (.not.associated(coupler)) then
            if (bcs_done) then
              exit
            else
              bcs_done = PETSC_TRUE
              if (associated(patch%source_sink_list)) then
                coupler => patch%source_sink_list%first
                if (.not.associated(coupler)) exit
              else
                exit
              endif
            endif
          endif
        
          string = trim(coupler%name) // ': Water(L)'
          call OutputWriteToHeader(fid,string,'-kg','',icol)
          string = trim(coupler%name) // ': Water(G)'
          call OutputWriteToHeader(fid,string,'-kg','',icol)
          string = trim(coupler%name) // ': Water(S)'
          call OutputWriteToHeader(fid,string,'-kg','',icol)

          units = 'kg/' // trim(output_option%tunit) // ''
          string = trim(coupler%name) // ': Water(L)'
          call OutputWriteToHeader(fid,string,units,'',icol)
          string = trim(coupler%name) // ': Water(G)'
          call OutputWriteToHeader(fid,string,units,'',icol)
          string = trim(coupler%name) // ': Water(S)'
          call OutputWriteToHeader(fid,string,units,'',icol)

          coupler => coupler%next
        enddo
        ! sum of total mass and bc/ss
        if (option%nflowdof > 0) then
          string = ' SUM of All Water Mass and Flux'
          call OutputWriteToHeader(fid,string,'kg','',icol)
        endif

        ! Print the water mass [kg] and species mass [mol] in the specified regions (header)
        if (associated(output_option%mass_balance_region_list)) then
          cur_mbr => output_option%mass_balance_region_list
          do
            if (.not.associated(cur_mbr)) exit
            string = 'Region ' // trim(cur_mbr%region_name) // ' Water(L)'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            string = 'Region ' // trim(cur_mbr%region_name) // ' Water(G)'
            call OutputWriteToHeader(fid,string,'kg','',icol)
            string = 'Region ' // trim(cur_mbr%region_name) // ' Water(S)'
            call OutputWriteToHeader(fid,string,'kg','',icol)


            cur_mbr => cur_mbr%next
          enddo
        endif
      
        write(fid,'(a)') ''
      
      else
        open(unit=fid,file=filename,action="write",status="old",position="append")
      endif
    
    endif

100 format(100es16.8)
110 format(100es16.8)

    ! ----------------------------------
    ! write time
    if (option%myrank == option%io_rank) then
      write(fid,100,advance="no") option%time/output_option%tconv
    endif
  
    if (option%nflowdof > 0) then
      if (option%myrank == option%io_rank) &
        write(fid,100,advance="no") option%flow%dt/output_option%tconv
    endif
  
    ! ----------------------------------
    ! print out mass balance data

    if (option%nflowdof > 0) then
      select type(realization_base)
        class is(realization_subsurface_type)
          select case(option%iflowmode)
            case(MPFLOW_MODE)
              call MpFlowComputeMassBalance(realization_base,mass_balance)
          end select
        class default
          option%io_buffer = 'Unrecognized realization class in MassBalance().'
          call printErrMsg(option)
      end select

      sum_local = mass_balance(:,ispec)
      sum_allwater = 0.d0

      int_mpi = option%flow%nspecflow*option%flow%nfluid
      call MPI_Reduce(sum_local,sum_global, &
                    int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                    option%io_rank,option%mycomm,ierr)

      if (option%myrank == option%io_rank) then
        do iphase = 1, option%flow%nfluid
          write(fid,110,advance="no") sum_global(iphase)
            !
          sum_allwater = sum_allwater + sum_global(iphase)
        enddo
      endif
    endif
  
    coupler => patch%boundary_condition_list%first
    mpflow_auxvars_bc_or_ss => patch%aux%MpFlow%auxvars_bc
    bcs_done = PETSC_FALSE
    do
      if (.not.associated(coupler)) then
        if (bcs_done) then
          exit
        else
          bcs_done = PETSC_TRUE
          if (associated(patch%source_sink_list)) then
            coupler => patch%source_sink_list%first
            if (.not.associated(coupler)) exit
            mpflow_auxvars_bc_or_ss => patch%aux%MpFlow%auxvars_ss
          else
            exit
          endif
        endif
      endif

      offset = coupler%connection_set%offset
      if (option%nflowdof > 0 &
         .and. option%compute_mass_balance_new) then
        ! print out cumulative H2O flux
        sum_local = 0.d0
        do iconn = 1, coupler%connection_set%num_connections
          sum_local = sum_local + mpflow_auxvars_bc_or_ss(offset+iconn)%mass_balance(:,ispec)
        enddo
        int_mpi = option%flow%nfluid
        call MPI_Reduce(sum_local,sum_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
        if (option%myrank == option%io_rank) then
          ! change sign for positive in / negative out
          write(fid,110,advance="no") -sum_global
          !
          do iphase = 1, option%flow%nfluid
              sum_allwater = sum_allwater+sum_global(iphase)
          enddo
        endif

        ! print out H2O flux (mass)
        sum_local = 0.d0
        do iconn = 1, coupler%connection_set%num_connections
          sum_local = sum_local + mpflow_auxvars_bc_or_ss(offset+iconn)%mass_balance_delta(:,ispec)
        enddo
        ! mass_balance_delta units = delta kmol h2o; must convert to delta kg h2o
        sum_local = sum_local*FMWH2O

        int_mpi = option%flow%nfluid
        call MPI_Reduce(sum_local,sum_global, &
                          int_mpi,MPI_DOUBLE_PRECISION,MPI_SUM, &
                          option%io_rank,option%mycomm,ierr)
                              
        if (option%myrank == option%io_rank) then
          ! change sign for positive in / negative out
          write(fid,110,advance="no") -sum_global*output_option%tconv   ! ? need checking here
        endif

      endif

      coupler => coupler%next
    enddo
    ! sum of all mass-balance items
    ! print out sum of global water mass and bc/ss fluxes
    ! (NOTE: THIS amount MUST BE constant, otherwise error in mass-conservation)
    if (option%myrank == option%io_rank) then
      write(fid,110,advance="no") sum_allwater
    endif
  
    ! Print the total water and component mass in the specified regions (data)
    if (associated(output_option%mass_balance_region_list)) then
      cur_mbr => output_option%mass_balance_region_list
      do
        if (.not.associated(cur_mbr)) exit
        call PatchGetWaterMassInRegion(cur_mbr%region_cell_ids, &
                                     cur_mbr%num_cells,patch,option, &
                                     global_water_mass)
        write(fid,110,advance="no") global_water_mass
        cur_mbr => cur_mbr%next
      enddo
    endif
  
    if (option%myrank == option%io_rank) then
      write(fid,'(a)') ''
      close(fid)
    endif
  
  end do  ! ispec =1, nspecflow

  mass_balance_first = PETSC_FALSE

end subroutine OutputMassBalance
! ************************************************************************** !
! ************************************************************************** !

end module Output_MassBalance_module
