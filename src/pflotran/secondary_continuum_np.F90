! Added by Amphos 21
! For technical details contact albert.nardi@amphos21.com
! ========================================================

module Secondary_Continuum_NP_module
  
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Secondary_Continuum_Aux_module

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none

  private

  ! secondary continuum cell type
  PetscInt, parameter, public :: SLAB = 0
  PetscInt, parameter, public :: NESTED_CUBE = 1
  PetscInt, parameter, public :: NESTED_SPHERE = 2

  PetscReal, parameter :: perturbation_tolerance = 1.d-5

  public :: SecondaryContinuumType, &
            SecondaryContinuumSetProperties, &
            SecondaryRTAuxVarInit, &
            SecondaryRTResJacMulti_NP, &
            SecondaryRTAuxVarComputeMulti, &
            THCSecHeatAuxVarCompute, &
            THSecHeatAuxVarCompute, &
            MphaseSecHeatAuxVarCompute, &
            SecondaryRTUpdateIterate_NP, &
            SecondaryRTUpdateEquilState, &
            SecondaryRTUpdateKineticState, &
            SecondaryRTTimeCut, &
            SecondaryRTGetVariable, &
            SecondaryRTSetVariable, &
            ComputeElectricPotentialTotalComponent_NP

contains

! ************************************************************************** !

subroutine SecondaryContinuumType(sec_continuum,nmat,aream, &
            volm,dm1,dm2,aperture,epsilon,log_spacing,outer_spacing, &
            interfacial_area,option)
  ! 
  ! The area, volume, grid sizes for secondary continuum
  ! are calculated based on the input dimensions and geometry
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/11/12
  ! 

  use Option_module

  implicit none
  
  type(sec_continuum_type) :: sec_continuum

  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: igeom, nmat, m
  PetscReal :: aream(nmat), volm(nmat), dm1(nmat), dm2(nmat)
  PetscReal :: dy, r0, r1, aream0, am0, vm0, interfacial_area
  PetscReal :: num_density, aperture, epsilon, fracture_spacing
  PetscReal :: outer_spacing, matrix_block_size
  PetscReal :: grid_spacing(nmat)
  PetscBool :: log_spacing
  PetscReal :: sum

  PetscInt, save :: icall

  data icall/0/

  igeom = sec_continuum%itype
  option%nsec_cells = nmat
    
  select case (igeom)      
    case(SLAB)
    
      dy = sec_continuum%slab%length/nmat
      aream0 = sec_continuum%slab%area
      do m = 1, nmat
        volm(m) = dy*aream0
      enddo
      am0 = aream0
      vm0 = nmat*dy*aream0
      interfacial_area = am0/vm0
     
      do m = 1, nmat
        aream(m) = aream0
        dm1(m) = 0.5d0*dy
        dm2(m) = 0.5d0*dy
      enddo

      if (icall == 0 .and. OptionPrintToFile(option)) then
        icall = 1
        string = 'DCDM Multiple Continuum Model'
        write(option%fid_out,'(/,2x,a,/)') trim(string)
        string = 'Slab'
        write(option%fid_out,'(2x,a,/)') trim(string)
        num_density = (1.d0-epsilon)/vm0
        write(option%fid_out,'(2x,"number density: ",11x,1pe12.4," m^(-3)")') num_density
        write(option%fid_out,'(2x,"matrix block size: ",8x,1pe12.4," m")') sec_continuum%slab%length
        write(option%fid_out,'(2x,"epsilon: ",18x,1pe12.4)') epsilon
        write(option%fid_out,'(2x,"specific interfacial area: ",1pe12.4," m^(-1)")') interfacial_area
        do m = 1, nmat
          if (m == 1) write(option%fid_out,'(/,2x,"node matrix volume fraction")') 
          write(option%fid_out,'(2x,i3,3x,1pe12.4)') m,volm(m)/vm0 !*(1.d0 - epsilon)
        enddo
!       aperture = r0*(1.d0/(1.d0-epsilon)**(1.d0/3.d0)-1.d0)
!       write(option%fid_out,'(2x,"aperture: ",17x,1pe12.4," m")') aperture
      endif
      
      ! Store the distances
      sec_continuum%distance(1) = dm1(1)
      do m = 2, nmat
        sec_continuum%distance(m) = sec_continuum%distance(m-1) + &
                                      dm2(m-1) + dm1(m)
      enddo
          
    case(NESTED_CUBE)

      if (sec_continuum%nested_cube%fracture_spacing > 0.d0) then

        fracture_spacing = sec_continuum%nested_cube%fracture_spacing
!        override epsilon if aperture defined
        if (aperture > 0.d0) then
          r0 = fracture_spacing - aperture
          epsilon = 1.d0 - (1.d0 + aperture/r0)**(-3.d0)
        else if (epsilon > 0.d0) then
          r0 = fracture_spacing*(1.d0-epsilon)**(1.d0/3.d0)
          aperture = r0*((1.d0-epsilon)**(-1.d0/3.d0)-1.d0)
        endif
                                            
      else if (sec_continuum%nested_cube%matrix_block_size > 0.d0) then

        r0 = sec_continuum%nested_cube%matrix_block_size

!        override epsilon if aperture defined
        if (aperture > 0.d0) then
          fracture_spacing = r0 + aperture
          epsilon = 1.d0 - (1.d0 + aperture/r0)**(-3.d0)
        else if (epsilon > 0.d0) then
          fracture_spacing = r0*(1.d0-epsilon)**(-1.d0/3.d0)
          aperture = fracture_spacing - r0
        endif
      endif
      
      if (log_spacing) then 
        
        matrix_block_size = r0
        call SecondaryContinuumCalcLogSpacing(matrix_block_size,outer_spacing, &
                                              nmat,grid_spacing,option)
        
        r0 = 2.d0*grid_spacing(1)
        dm1(1) = 0.5d0*grid_spacing(1)
        dm2(1) = 0.5d0*grid_spacing(1)
        volm(1) = r0**3.d0
        aream(1) = 6.d0*r0**2.d0
        do m = 2, nmat
          dm1(m) = 0.5d0*grid_spacing(m)
          dm2(m) = 0.5d0*grid_spacing(m)
          r1 = r0 + 2.d0*(dm1(m) + dm2(m))
          volm(m) = r1**3.d0 - r0**3.d0
          aream(m) = 6.d0*r1**2.d0
          r0 = r1
        enddo
        r0 = matrix_block_size
        am0 = 6.d0*r0**2.d0
        vm0 = r0**3.d0
        interfacial_area = am0/vm0

      else
        dy = r0/nmat/2.d0
     
        r0 = 2.d0*dy
        volm(1) = r0**3.d0
        do m = 2, nmat
          r1 = r0 + 2.d0*dy
          volm(m) = r1**3.d0 - r0**3.d0
          r0 = r1
        enddo
      
        r0 = 2.d0*dy
        aream(1) = 6.d0*r0**2.d0
        dm1(1) = 0.5d0*dy
        dm2(1) = 0.5d0*dy
        do m = 2, nmat
          dm1(m) = 0.5d0*dy
          dm2(m) = 0.5d0*dy
          r0 = r0 + 2.d0*dy
          aream(m) = 6.d0*r0**2.d0
        enddo
        r0 = real(2.d0*nmat)*dy
        am0 = 6.d0*r0**2.d0
        vm0 = r0**3.d0
        interfacial_area = am0/vm0
      endif

      if (icall == 0 .and. OptionPrintToFile(option)) then
        icall = 1
        string = 'DCDM Multiple Continuum Model'
        write(option%fid_out,'(/,2x,a,/)') trim(string)
        string = 'Nested Cubes'
        write(option%fid_out,'(2x,a,/)') trim(string)
        num_density = (1.d0-epsilon)/vm0
        write(option%fid_out,'(2x,"number density: ",11x,1pe12.4," m^(-3)")') num_density
        write(option%fid_out,'(2x,"matrix block size: ",8x,1pe12.4," m")') r0
        write(option%fid_out,'(2x,"epsilon: ",18x,1pe12.4)') epsilon
        write(option%fid_out,'(2x,"specific interfacial area: ",1pe12.4," m^(-1)")') interfacial_area
        write(option%fid_out,'(2x,"fracture aperture: ",8x,1pe12.4," m")') aperture
        write(option%fid_out,'(2x,"fracture spacing: ",9x,1pe12.4," m")') fracture_spacing
        write(option%fid_out,'(/,2x,"node  vol. frac.      dm1         dm2         aream       dy          y")')
        r0 = 0.d0
        do m = 1, nmat
          if (m == 1) then
            r0 = r0 + dm1(m)
          else
            r0 = r0 + dm2(m-1)+dm1(m)
          endif
          write(option%fid_out,'(2x,i3,3x,1p6e12.4)') m,volm(m)/vm0,dm1(m),dm2(m),aream(m), &
          dm1(m)+dm2(m),r0
        enddo
      endif

      ! Store the distances
      sec_continuum%distance(1) = dm1(1)
      do m = 2, nmat
        sec_continuum%distance(m) = sec_continuum%distance(m-1) + &
                                      dm2(m-1) + dm1(m)
      enddo     

    case(NESTED_SPHERE)
    
      dy = sec_continuum%nested_sphere%radius/nmat
      r0 = dy

      volm(1) = 4.d0/3.d0*pi*r0**3.d0
      do m = 2, nmat
        r1 = r0 + dy
        volm(m) = 4.d0/3.d0*pi*(r1**3.d0 - r0**3.d0)
        r0 = r1
      enddo
      
      r0 = dy
      aream(1) = 4.d0*pi*r0**2.d0
      dm1(1) = 0.5d0*dy
      dm2(1) = 0.5d0*dy
      do m = 2, nmat
        r0 = r0 + dy
        dm1(m) = 0.5d0*dy
        dm2(m) = 0.5d0*dy
        aream(m) = 4.d0*pi*r0**2.d0
      enddo
      r0 = 0.5d0*real(2.d0*nmat)*dy
      am0 = 4.d0*pi*r0**2.d0
      vm0 = am0*r0/3.d0
      interfacial_area = am0/vm0

      if (icall == 0 .and. OptionPrintToFile(option)) then
        icall = 1
        string = 'DCDM Multiple Continuum Model'
        write(option%fid_out,'(/,2x,a,/)') trim(string)
        string = 'Nested Spheres'
        write(option%fid_out,'(2x,a,/)') trim(string)
        num_density = (1.d0-epsilon)/vm0
        write(option%fid_out,'(2x,"number density: ",11x,1pe12.4," m^(-3)")') num_density
        write(option%fid_out,'(2x,"sphere radius: ",8x,1pe12.4," m")') sec_continuum%nested_sphere%radius
        write(option%fid_out,'(2x,"epsilon: ",18x,1pe12.4)') epsilon
        write(option%fid_out,'(2x,"specific interfacial area: ",1pe12.4," m^(-1)")') interfacial_area
        do m = 1, nmat
          if (m == 1) write(option%fid_out,'(/,2x,"node matrix volume fraction")') 
          write(option%fid_out,'(2x,i3,3x,1pe12.4)') m,volm(m)/vm0*(1.d0 - epsilon)
        enddo

!       aperture = r0*(1.d0/(1.d0-epsilon)**(1.d0/3.d0)-1.d0)
!       write(option%fid_out,'(2x,"aperture: ",17x,1pe12.4," m")') aperture
      endif
      
      ! Store the distances
      sec_continuum%distance(1) = dm1(1)
      do m = 2, nmat
        sec_continuum%distance(m) = sec_continuum%distance(m-1) + &
                                      dm2(m-1) + dm1(m)
      enddo
                        
  end select
  
  
  sum = 0.d0
  do m = 1,nmat
    if (volm(m)/vm0 > 1.d0) then
      print *, 'Error: volume fraction for cell', m, 'is greater than 1.'
      stop
    else 
      sum = sum + volm(m)/vm0
    endif
  enddo
  
  if (icall /= 2 .and. OptionPrintToFile(Option)) then
    icall = 2
    write(option%fid_out,'(/,"sum of volume fractions:",1x,1pe12.4)') sum
  endif
  
  if (abs(sum - 1.d0) > 1.d-6) then
    option%io_buffer = 'Error: Sum of the volume fractions of the' // &
                       ' secondary cells is not equal to 1.'
    call PrintErrMsg(option)
  endif
  
end subroutine SecondaryContinuumType

! ************************************************************************** !

subroutine SecondaryContinuumSetProperties(sec_continuum, &
                                           sec_continuum_name, & 
                                           sec_continuum_length, &
                                           sec_continuum_matrix_block_size, &
                                           sec_continuum_fracture_spacing, &
                                           sec_continuum_radius, &
                                           sec_continuum_area, &
                                           option)
  ! 
  ! The type, dimensions of the secondary
  ! continuum are set
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 07/17/12
  ! 
                                    
  use Option_module
  use String_module
  
  implicit none
  
  type(sec_continuum_type) :: sec_continuum
  type(option_type) :: option
  PetscReal :: sec_continuum_matrix_block_size
  PetscReal :: sec_continuum_fracture_spacing
  PetscReal :: sec_continuum_length
  PetscReal :: sec_continuum_area
  PetscReal :: sec_continuum_radius
  character(len=MAXWORDLENGTH) :: sec_continuum_name

  call StringToUpper(sec_continuum_name)
  
  select case(trim(sec_continuum_name))
    case("SLAB")
      sec_continuum%itype = SLAB
      sec_continuum%slab%length = sec_continuum_length
      if (Equal(sec_continuum_area,0.d0)) then
        option%io_buffer = 'Keyword "AREA" not specified for SLAB type ' // &
                           'under SECONDARY_CONTINUUM'
        call PrintErrMsg(option)
      endif
      sec_continuum%slab%area = sec_continuum_area
    case("NESTED_CUBES")
      sec_continuum%itype = NESTED_CUBE
      sec_continuum%nested_cube%matrix_block_size = sec_continuum_matrix_block_size
      sec_continuum%nested_cube%fracture_spacing = sec_continuum_fracture_spacing
    case("NESTED_SPHERES")
      sec_continuum%itype = NESTED_SPHERE
      sec_continuum%nested_sphere%radius = sec_continuum_radius
    case default
      option%io_buffer = 'Keyword "' // trim(sec_continuum_name) // '" not ' // &
                         'recognized in SecondaryContinuumSetProperties()'
      call PrintErrMsg(option)
  end select
      
end subroutine SecondaryContinuumSetProperties  

! ************************************************************************** !

subroutine SecondaryContinuumCalcLogSpacing(matrix_size,outer_grid_size, &
                                            sec_num_cells,grid_spacing,option)
  ! 
  ! Given the matrix block size and the
  ! grid spacing of the outer most secondary continuum cell, a geometric
  ! series is assumed and the grid spacing of the rest of the cells is
  ! calculated
  ! 
  ! Equation:
  ! \frac{1 - \rho}{1 - \rho_M}*\rho*(M-1) = \frac{2\Delta\xi_m}{l_M}
  !
  ! where
  !   \Delta\xi_m: Grid spacing of the outer most continuum cell (INPUT)
  !   l_M        : Matrix block size (INPUT)
  !   M          : Number of secondary continuum cells (INPUT)
  !   \rho       : Logarithmic grid spacing factor (COMPUTED)
  !
  ! Author: Satish Karra, LANL
  ! Date: 07/17/12
  ! 
                                              
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal :: matrix_size, outer_grid_size
  PetscInt :: sec_num_cells
  PetscReal :: grid_spacing(sec_num_cells)
  PetscReal :: delta, delta_new, inner_grid_size
  PetscReal :: F, dF
  PetscReal, parameter :: tol = 1.d-12
  PetscInt, parameter :: maxit = 50
  PetscInt :: i 
  
  
  if (mod(sec_num_cells,2) /= 0) then
     option%io_buffer = 'NUM_CELLS under SECONDARY_CONTINUUM has to be' // &
                        ' even for logarithmic grid spacing'
      call PrintErrMsg(option)
  endif
  
  delta = 0.99d0
  
  do i = 1, maxit
    F = (1.d0 - delta)/(1.d0 - delta**sec_num_cells)*delta**(sec_num_cells - 1.d0) - &
        2.d0*outer_grid_size/matrix_size
    dF = (1.d0 + sec_num_cells*(delta - 1.d0) - delta**sec_num_cells)/ &
         (delta**sec_num_cells - 1.d0)**2.d0*delta**(sec_num_cells - 2.d0)
    delta_new = delta + F/dF
    if ((abs(F) < tol)) exit
    delta = delta_new
    if (delta < 0.d0) delta = 0.5d0
!   if (delta > 1.d0) delta = 0.9d0
  enddo
  
  if (i == maxit) then
     option%io_buffer = 'Log Grid spacing solution has not converged' // &
                        ' with given fracture values.'
     call PrintErrMsg(option)
  endif

  inner_grid_size = outer_grid_size/delta**(sec_num_cells - 1)
  
  do i = 1, sec_num_cells
    grid_spacing(i) = inner_grid_size*delta**(i-1)
  enddo

!  write(option%fid_out,'("  Logarithmic grid spacing: delta = ",1pe12.4)') delta
    
end subroutine SecondaryContinuumCalcLogSpacing

! ************************************************************************** !

subroutine SecondaryRTTimeCut(realization)
  ! 
  ! Resets secondary concentrations to previous time
  ! step when there is a time cut
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/29/13
  ! 

  use Realization_Subsurface_class
  use Grid_module
  use Reaction_Aux_module
  
  implicit none
  class(realization_subsurface_type) :: realization
  class(reaction_rt_type), pointer :: reaction
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  type(grid_type), pointer :: grid
  
  PetscInt :: local_id, ghosted_id
  PetscInt :: ngcells, ncomp
  PetscInt :: cell, comp

  reaction => realization%reaction
  rt_sec_transport_vars => realization%patch%aux%SC_RT%sec_transport_vars
  grid => realization%patch%grid  
  
  ncomp = reaction%naqcomp
  
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (realization%patch%imat(ghosted_id) <= 0) cycle
    do comp = 1, ncomp
      ngcells = rt_sec_transport_vars(ghosted_id)%ncells
      do cell = 1, ngcells
        rt_sec_transport_vars(ghosted_id)%updated_conc(comp,cell) = &
          rt_sec_transport_vars(ghosted_id)%sec_rt_auxvar(cell)%pri_molal(comp)
      enddo
    enddo
  enddo
 
end subroutine SecondaryRTTimeCut

! ************************************************************************** !

subroutine SecondaryRTAuxVarInit(multicontinuum,epsilon,rt_sec_transport_vars,reaction, &
                                 initial_condition,constraint,option)
  ! 
  ! Initializes all the secondary continuum reactive
  ! transport variables
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 02/05/13
  ! 
  
  use Coupler_module
  use Transport_Constraint_module
  use Condition_module
  use Global_Aux_module
  use Material_module
  use Option_module
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_module
  use Transport_Constraint_RT_module
  
  use EOS_Water_module
  
  implicit none 
  
  type(sec_transport_type) :: rt_sec_transport_vars
  type(multicontinuum_property_type) :: multicontinuum
  class(reaction_rt_type), pointer :: reaction
  type(coupler_type), pointer :: initial_condition
  type(option_type), pointer :: option
  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  class(material_auxvar_type), allocatable :: material_auxvar
  class(tran_constraint_rt_type), pointer :: constraint
  type(flow_condition_type), pointer :: initial_flow_condition
  

  PetscReal :: equil_conc(reaction%mineral%nmnrl)
  PetscInt :: i, cell
  PetscReal :: area_per_vol
  PetscReal :: dum1
  PetscReal :: epsilon
  PetscInt :: num_iterations
  PetscErrorCode :: ierr
  
  num_iterations = 0

  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  material_auxvar%porosity = option%flow%reference_porosity

  call SecondaryContinuumSetProperties( &
        rt_sec_transport_vars%sec_continuum, &
        multicontinuum%name, &
        multicontinuum%length, &
        multicontinuum%matrix_block_size, &
        multicontinuum%fracture_spacing, &
        multicontinuum%radius, &
        multicontinuum%area, &
        option)
        
  rt_sec_transport_vars%ncells = multicontinuum%ncells
  rt_sec_transport_vars%aperture = multicontinuum%aperture
  rt_sec_transport_vars%epsilon = epsilon 
  rt_sec_transport_vars%log_spacing = multicontinuum%log_spacing
  rt_sec_transport_vars%outer_spacing = multicontinuum%outer_spacing    
        
  allocate(rt_sec_transport_vars%area(rt_sec_transport_vars%ncells))
  allocate(rt_sec_transport_vars%vol(rt_sec_transport_vars%ncells))
  allocate(rt_sec_transport_vars%dm_minus(rt_sec_transport_vars%ncells))
  allocate(rt_sec_transport_vars%dm_plus(rt_sec_transport_vars%ncells))
  allocate(rt_sec_transport_vars%sec_continuum% &
             distance(rt_sec_transport_vars%ncells))
    
  call SecondaryContinuumType(rt_sec_transport_vars%sec_continuum, &
                              rt_sec_transport_vars%ncells, &
                              rt_sec_transport_vars%area, &
                              rt_sec_transport_vars%vol, &
                              rt_sec_transport_vars%dm_minus, &
                              rt_sec_transport_vars%dm_plus, &
                              rt_sec_transport_vars%aperture, &
                              rt_sec_transport_vars%epsilon, &
                              rt_sec_transport_vars%log_spacing, &
                              rt_sec_transport_vars%outer_spacing, &
                              area_per_vol,option)                                
  rt_sec_transport_vars%interfacial_area = area_per_vol* &
         (1.d0 - rt_sec_transport_vars%epsilon)*multicontinuum% &
         area_scaling
  
  ! Initializing the secondary RT auxvars
  allocate(rt_sec_transport_vars%sec_rt_auxvar(rt_sec_transport_vars%ncells))
  do cell = 1, rt_sec_transport_vars%ncells
    call RTAuxVarInit(rt_sec_transport_vars%sec_rt_auxvar(cell),reaction,option)
  enddo

  allocate(rt_sec_transport_vars%sec_jac(reaction%naqcomp,reaction%naqcomp))    
           
  ! Allocate diagonal terms
  allocate(rt_sec_transport_vars%cxm(reaction%naqcomp,reaction%naqcomp,&
           rt_sec_transport_vars%ncells)) 
  allocate(rt_sec_transport_vars%cxp(reaction%naqcomp,reaction%naqcomp,&
           rt_sec_transport_vars%ncells))  
  allocate(rt_sec_transport_vars%cdl(reaction%naqcomp,reaction%naqcomp,&
           rt_sec_transport_vars%ncells)) 
  allocate(rt_sec_transport_vars% &
           r(reaction%naqcomp*rt_sec_transport_vars%ncells))
  allocate(rt_sec_transport_vars% &
           updated_conc(reaction%naqcomp,rt_sec_transport_vars%ncells))
           
  
  initial_flow_condition => initial_condition%flow_condition
  do cell = 1, rt_sec_transport_vars%ncells
    global_auxvar => initial_condition%tran_condition% &
                       constraint_coupler_list%global_auxvar
    rt_auxvar => rt_sec_transport_vars%sec_rt_auxvar(cell)
    if (associated(initial_flow_condition)) then
      if (associated(initial_flow_condition%pressure)) then
        if (associated(initial_flow_condition%pressure%dataset)) then
          global_auxvar%pres = &
            initial_flow_condition%pressure%dataset%rarray(1)
        else
          global_auxvar%pres = option%flow%reference_pressure
        endif
      else 
        global_auxvar%pres = option%flow%reference_pressure
      endif
      if (associated(initial_flow_condition%temperature)) then
        if (associated(initial_flow_condition%temperature%dataset)) then
          global_auxvar%temp  = &
            initial_flow_condition%temperature%dataset%rarray(1)
        else
          global_auxvar%temp = option%flow%reference_temperature
        endif
      else
        global_auxvar%temp = option%flow%reference_temperature
      endif
        
      call EOSWaterDensity(global_auxvar%temp, &
                           global_auxvar%pres(1), &
                           global_auxvar%den_kg(1), &
                           dum1,ierr)
    else
      global_auxvar%pres = option%flow%reference_pressure
      global_auxvar%temp = option%flow%reference_temperature
      global_auxvar%den_kg(option%liquid_phase) = &
        option%flow%reference_density(option%liquid_phase)

    endif
    global_auxvar%sat = option%flow%reference_saturation

    if (option%transport%nphase > option%nphase) then
      ! gas phase not considered explicitly on flow side
      global_auxvar%den_kg(option%gas_phase) = &
        option%flow%reference_density(option%gas_phase)
      global_auxvar%sat(option%gas_phase) = &
        1.d0 - global_auxvar%sat(option%liquid_phase)
    endif

    !Use multicontinuum sorption
    reaction%mc_flag = 1
    call ReactionEquilibrateConstraint(rt_auxvar,global_auxvar, &
                          material_auxvar, &
                          reaction,constraint, &
                          num_iterations, &
                          PETSC_FALSE,option)   
    reaction%mc_flag = 0
    
    rt_sec_transport_vars%updated_conc(:,cell) =  rt_auxvar%pri_molal   
       
  enddo                                    
  
  call MaterialAuxVarStrip(material_auxvar)
  deallocate(material_auxvar)
  
  rt_sec_transport_vars%sec_jac_update = PETSC_FALSE
  rt_sec_transport_vars%sec_jac = 0.d0
  rt_sec_transport_vars%cxm = 0.d0
  rt_sec_transport_vars%cxp = 0.d0
  rt_sec_transport_vars%cdl = 0.d0
  rt_sec_transport_vars%r = 0.d0
      
end subroutine SecondaryRTAuxVarInit  

! ************************************************************************** !

subroutine SecondaryRTResJacMulti_NP(sec_transport_vars,auxvar, &
                                  global_auxvar,prim_vol, &
                                  reaction,rt_parameter,diffusion_coefficient, &
                                  porosity,tortuosity,option,res_transport)
  ! 
  ! RTSecondaryTransportMulti:  Calculates the source term contribution due to
  ! secondary continuum in the primary continuum residual for multicomponent
  ! system assuming only aqueous reaction
  ! 
  ! Author: Albert Nardi, Amphos 21
  ! Date: 8/31/2021
  ! 
                               
                            
  use Option_module 
  use Global_Aux_module
  use Block_Solve_module
  use Block_Tridiag_module
  use Utility_module
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  type(sec_transport_type) :: sec_transport_vars
  type(reactive_transport_auxvar_type) :: auxvar
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type), pointer :: reaction
  type(option_type) :: option
  PetscReal :: coeff_left(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_diag(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right(reaction%naqcomp,reaction%naqcomp, &
                           sec_transport_vars%ncells)
  PetscReal :: res(sec_transport_vars%ncells*reaction%naqcomp)

  PetscReal :: curConc_up, curConc_dn, curConc 
  PetscReal :: curCharge
  PetscReal :: cur_diff 
  PetscReal :: curStoich 
  PetscReal :: sum_transp_up, sum_transp_dn
  PetscReal :: sum_denom_up, sum_denom_dn
  PetscReal :: factor_up_em(reaction%naqcomp) 
  PetscReal :: factor_dn_em(reaction%naqcomp) 
  PetscReal :: potent_up(reaction%naqcomp)
  PetscReal :: potent_dn(reaction%naqcomp)

  PetscReal :: rhs(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: D_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: identity(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: b_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: sec_jac(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: inv_D_M(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: conc_upd(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: total_upd(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: total_prev(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: conc_current_M(reaction%naqcomp)
  PetscReal :: total_current_M(reaction%naqcomp)
  PetscReal :: res_transport(reaction%naqcomp)
  PetscReal :: total_primary_node(reaction%naqcomp)
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)
  PetscReal :: res_react(reaction%naqcomp)
  PetscReal :: jac_react(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: dtotal(reaction%naqcomp,reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: dtotal_prim(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: mc_pri_molal(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: mc_sec_molal(reaction%neqcplx,sec_transport_vars%ncells)
  PetscInt :: i, j, k, n, l
  PetscInt :: ngcells, ncomp
  PetscInt :: ni
  PetscInt :: nicomp, icomp
  PetscReal :: area_fm
  PetscReal :: diffusion_coefficient
  PetscReal :: porosity
  PetscReal :: tortuosity
  PetscReal :: dC_prim(reaction%naqcomp,reaction%naqcomp)
  
  PetscReal :: arrhenius_factor
  PetscReal :: pordt, pordiff
  PetscReal :: pordiff_prim, pordiff_sec
  PetscReal :: portort
  PetscReal :: prim_vol ! volume of primary grid cell
  PetscReal :: dCsec_dCprim(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: dPsisec_dCprim(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: jcomp, lcomp, kcomp, icplx, ncompeq  
  PetscReal :: sec_sec_molal_M(reaction%neqcplx)   ! secondary species molality of secondary continuum
  
  PetscInt :: pivot(reaction%naqcomp,sec_transport_vars%ncells)
  PetscInt :: indx(reaction%naqcomp)
  PetscInt :: d, ier
  PetscReal :: m
  
  ! Quantities for numerical jacobian
  PetscReal :: conc_prim(reaction%naqcomp)
  PetscReal :: conc_prim_pert(reaction%naqcomp)
  PetscReal :: sec_jac_num(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: conc_current_M_pert(reaction%naqcomp)
  PetscReal :: total_current_M_pert(reaction%naqcomp)
  PetscReal :: res_transport_pert(reaction%naqcomp)
  PetscReal :: total_primary_node_pert(reaction%naqcomp)
  PetscReal :: dtotal_prim_num(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: dPsisec_dCprim_num(reaction%naqcomp,reaction%naqcomp)
  PetscReal :: pert
  PetscReal :: coeff_diag_dm(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_left_dm(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right_dm(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_left_pert(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_diag_pert(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right_pert(reaction%naqcomp,reaction%naqcomp, &
                           sec_transport_vars%ncells)
  PetscReal :: coeff_left_copy(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_diag_copy(reaction%naqcomp,reaction%naqcomp, &
                          sec_transport_vars%ncells)
  PetscReal :: coeff_right_copy(reaction%naqcomp,reaction%naqcomp, &
                           sec_transport_vars%ncells)

  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: lnQK, tempreal 
  PetscReal :: J_up, J_dn

  PetscReal :: total_sorb_upd(reaction%naqcomp,sec_transport_vars%ncells) 
  PetscReal :: total_sorb_prev(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: dtotal_sorb_upd(reaction%naqcomp,reaction%naqcomp,sec_transport_vars%ncells)

  class(material_auxvar_type), allocatable :: material_auxvar
  
  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  vol = sec_transport_vars%vol          
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus
  area_fm = sec_transport_vars%interfacial_area
  ncomp = reaction%naqcomp

  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  material_auxvar%porosity = porosity
  
  do j = 1, ncomp
    do i = 1, ngcells
      total_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total(j,1)
      if (reaction%neqsorb > 0) then 
        total_sorb_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total_sorb_eq(j)
      endif
    enddo
  enddo
  conc_upd = sec_transport_vars%updated_conc

    
  ! Note that sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j) units are in mol/kg
  ! Need to convert to mol/L since the units of total. in the Thomas 
  ! algorithm are in mol/L 
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  res = 0.d0
  res_transport = 0.d0
  res_transport_pert = 0.d0
  dC_prim = 0.0 
  rhs = 0.d0
  D_M = 0.d0
  identity = 0.d0
  b_M = 0.d0
  inv_D_M = 0.d0
  total_current_M = 0.d0
  dPsisec_dCprim = 0.d0
  dCsec_dCprim = 0.d0
  
  total_primary_node = auxvar%total(:,1) ! in mol/L 
  dtotal_prim = auxvar%aqueous%dtotal(:,:,1)

  ! Compute totals and its derivatives total_upd dtotal
  call RTAuxVarInit(rt_auxvar,reaction,option)
  do i = 1, ngcells
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i),option)
    rt_auxvar%pri_molal = conc_upd(:,i)
    mc_pri_molal(:, i) = conc_upd(:,i)
    
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    if (reaction%neqsorb > 0) then
       call RTotalSorb(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                       reaction%isotherm%multicontinuum_isotherm_rxn,option)
    endif
    total_upd(:,i) = rt_auxvar%total(:,1)           ! phase 1 liquid
    dtotal(:,:,i) = rt_auxvar%aqueous%dtotal(:,:,1) ! phase 1 liquid
    mc_sec_molal(:, i) = rt_auxvar%sec_molal(:)
    if (reaction%neqsorb > 0) then 
      total_sorb_upd(:,i) = rt_auxvar%total_sorb_eq(:)
      dtotal_sorb_upd(:,:,i) = rt_auxvar%dtotal_sorb_eq(:,:)
    endif
  enddo 
                          
!================ Calculate the secondary residual =============================        

  pordt = porosity/option%tran_dt*1d3

  do i = 1, ngcells
  
    ! D_j \nabla C_j
    do icomp = 1, ncomp
      n = icomp + (i-1)*ncomp
      
      pordiff = porosity*diffusion_coefficient*tortuosity*global_auxvar%den_kg(1)
      pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(icomp)*&
              tortuosity*global_auxvar%den_kg(1)
      
      ! Accumulation        
      res(n) = pordt*(total_upd(icomp,i) - total_prev(icomp,i))*vol(i)    ! in mol/L*m3/s
      if (reaction%neqsorb > 0) then 
        res(n) = res(n) + vol(i)/option%tran_dt*(total_sorb_upd(icomp,i) - total_sorb_prev(icomp,i))
      endif      

      ! Flux terms
      if (i.gt.1.and.i.lt.ngcells) then
        res(n) = res(n) - pordiff_prim*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                        (mc_pri_molal(icomp,i+1) - &
                        mc_pri_molal(icomp,i))
        res(n) = res(n) + pordiff_prim*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                        (mc_pri_molal(icomp,i) - &
                        mc_pri_molal(icomp,i-1))

      ! Apply boundary conditions
      ! Inner boundary
      else if (i.eq.1) then
        res(n) = res(n) - pordiff_prim*area(1)/(dm_minus(2) + dm_plus(1))* &
                        (mc_pri_molal(icomp,i+1) - &
                        mc_pri_molal(icomp,i))
           
      ! Outer boundary -- closest to primary node
      else !if (i.eq.ngcells) then
        res(n) = res(n) - &
                 pordiff_prim*area(ngcells)/dm_plus(ngcells)* &
                 (auxvar%pri_molal(icomp) - mc_pri_molal(icomp,i))
        res(n) = res(n) + &
                 pordiff_prim*area(ngcells-1)/(dm_minus(ngcells) &
                 + dm_plus(ngcells-1))* &
                 (mc_pri_molal(icomp,i) - mc_pri_molal(icomp,i-1))
      endif  

    enddo
    

    ! summatory D_i \nabla C_i   
    do icplx = 1, reaction%neqcplx
    
      pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)*&
    	       tortuosity*global_auxvar%den_kg(1)
      nicomp = reaction%eqcplxspecid(0,icplx)
      
      do ni = 1, nicomp
        icomp = reaction%eqcplxspecid(ni,icplx)        
        n = icomp + (i-1)*ncomp

        ! Flux terms
        if (i.gt.1.and.i.lt.ngcells) then
          res(n) = res(n) - pordiff_sec*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                   (mc_sec_molal(icplx,i+1) - &
                   mc_sec_molal(icplx,i))* &
                   reaction%eqcplxstoich(ni,icplx)
          res(n) = res(n) + pordiff_sec*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                   (mc_sec_molal(icplx,i) - &
                   mc_sec_molal(icplx,i-1))* &
                   reaction%eqcplxstoich(ni,icplx)
        
        ! Apply boundary conditions
        ! Inner boundary
        else if (i.eq.1) then
          res(n) = res(n) - pordiff_sec*area(1)/(dm_minus(2) + dm_plus(1))* &
                   (mc_sec_molal(icplx,i+1)- &
                   mc_sec_molal(icplx,i))* &
                   reaction%eqcplxstoich(ni,icplx)
        
        ! Outer boundary -- closest to primary node
        else !if (i.eq.ngcells) then
          res(n) = res(n) - &
                   pordiff_sec*area(ngcells)/dm_plus(ngcells)* &
                   (auxvar%sec_molal(icplx) - &
                   mc_sec_molal(icplx,i))* & 
                   reaction%eqcplxstoich(ni,icplx) 
          res(n) = res(n) + pordiff_sec* &
                   area(ngcells-1)/(dm_minus(ngcells) + dm_plus(ngcells-1))* &
                   (mc_sec_molal(icplx,i) - &
                   mc_sec_molal(icplx,i-1))* &
                   reaction%eqcplxstoich(ni,icplx)
        endif
      enddo
    enddo    
    
    ! Electromigration term
    call ComputeElectricPotentialTotalComponent_NP(i,reaction, &
                                            rt_parameter,  &
                                            rt_auxvar, &
                                            sec_transport_vars, &
                                            icomp,auxvar, &
                                            potent_dn, &
                                            potent_up, &
                                            mc_pri_molal, &
                                            mc_sec_molal)

    sum_denom_up = 1d-40
    sum_denom_dn = 1d-40
    
    do jcomp = 1, rt_parameter%naqcomp
      curCharge = reaction%primary_spec_Z(jcomp)
      sum_denom_up = sum_denom_up + potent_up(jcomp)*curCharge
      sum_denom_dn = sum_denom_dn + potent_dn(jcomp)*curCharge
    enddo

    do jcomp = 1, rt_parameter%naqcomp
      factor_up_em(jcomp) = potent_up(jcomp) / sum_denom_up
      factor_dn_em(jcomp) = potent_dn(jcomp) / sum_denom_dn
    enddo

    portort = porosity*tortuosity*global_auxvar%den_kg(1)
    
    sum_transp_dn = 0d0
    sum_transp_up = 0d0
    
    ! Primary l index
    do lcomp = 1, rt_parameter%naqcomp
      curCharge = reaction%primary_spec_Z(lcomp)
      cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)
      
      ! Flux terms
      if (i.gt.1.and.i.lt.ngcells) then
        curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                     (mc_pri_molal(lcomp, i+1) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                     (mc_pri_molal(lcomp, i) - &
                     mc_pri_molal(lcomp, i-1))      
      
      ! Apply boundary conditions
      ! Inner boundary
      else if (i.eq.1) then
        curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                     (mc_pri_molal(lcomp, i+1) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = 0.d0
      
      ! Outer boundary    
      else !if (i.eq.ngcells) then
        curConc_up = area(ngcells)/dm_plus(ngcells)* &
                     (auxvar%pri_molal(lcomp) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                     + dm_plus(ngcells-1))* &
                     (mc_pri_molal(lcomp, i) - &
                     mc_pri_molal(lcomp, i-1))
      endif
      
      sum_transp_up = sum_transp_up - curCharge*cur_diff*(curConc_up)
      sum_transp_dn = sum_transp_dn + curCharge*cur_diff*(curConc_dn)

      ! Secondary term k
      do icplx = 1, reaction%neqcplx
        nicomp = reaction%eqcplxspecid(0,icplx)
        do ni = 1, nicomp
          kcomp = reaction%eqcplxspecid(ni,icplx)
          if (kcomp == lcomp) then
            curCharge = reaction%primary_spec_Z(kcomp)
            cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
            curStoich = reaction%eqcplxstoich(ni,icplx)
            
            ! Flux terms
            if (i.gt.1.and.i.lt.ngcells) then
              curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                          (mc_sec_molal(icplx, i+1) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                          (mc_sec_molal(icplx, i) - &
                          mc_sec_molal(icplx, i-1))
            
            ! Apply boundary conditions
            ! Inner boundary
            else if (i.eq.1) then
              curConc_up = area(1)/(dm_minus(2) + dm_plus(1))* &
                          (mc_sec_molal(icplx, i+1) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = 0.d0
          
            ! Outer boundary   
            else !if (i.eq.ngcells) then 
              curConc_up = area(ngcells)/dm_plus(ngcells)* &
                          (auxvar%sec_molal(icplx) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                          + dm_plus(ngcells-1))* &
                          (mc_sec_molal(icplx, i) - &
                          mc_sec_molal(icplx, i-1))
            endif
            
            sum_transp_up = sum_transp_up - curCharge*cur_diff*curStoich*(curConc_up)
            sum_transp_dn = sum_transp_dn + curCharge*cur_diff*curStoich*(curConc_dn)
          endif
      
        enddo
      enddo
  
    enddo

    do icomp = 1, rt_parameter%naqcomp
      n = icomp + (i-1)*ncomp
      res(n) = res(n) - factor_dn_em(icomp)*portort*sum_transp_dn
      res(n) = res(n) - factor_up_em(icomp)*portort*sum_transp_up
    enddo
    
  enddo              
  
  !res = res*1.d3 ! Convert mol/L*m3/s to mol/s    

!================ Calculate the secondary jacobian =============================        

  do i = 1, ngcells
  
    ln_conc = log(rt_auxvar%pri_molal)
    ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
    lnQK = 0d0
    
    ! Accumulation
    do j = 1, ncomp
      do k = 1, ncomp 
        coeff_diag(j,k,i) = coeff_diag(j,k,i) + pordt*vol(i)*dtotal(j,k,i)
        if (reaction%neqsorb > 0) then
          coeff_diag(j,k,i) = coeff_diag(j,k,i) + vol(i)/option%tran_dt*(dtotal_sorb_upd(j,k,i))
        endif
      enddo
    enddo


    ! d (D_j * C_j) / dC_j =  D_j
    !do j = 1, rt_parameter%naqcomp
    do j = 1, ncomp
      do k = 1, ncomp 
        pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(j)*tortuosity &
                  *global_auxvar%den_kg(1)

        ! Flux terms
        if (i.gt.1.and.i.lt.ngcells) then
          coeff_diag(j,k,i) = coeff_diag(j,k,i) + &
                              pordiff_prim*area(i)/(dm_minus(i+1) + dm_plus(i)) + &
                              pordiff_prim*area(i-1)/(dm_minus(i) + dm_plus(i-1))
          coeff_left(j,k,i) = coeff_left(j,k,i) - &
                              pordiff_prim*area(i-1)/(dm_minus(i) + dm_plus(i-1))
          coeff_right(j,k,i) = coeff_right(j,k,i) - &
                              pordiff_prim*area(i)/(dm_minus(i+1) + dm_plus(i))                          
        
        ! Apply boundary conditions
        ! Inner boundary
        else if (i.eq.1) then
          coeff_diag(j,k,1) = coeff_diag(j,k,1) + &
                              pordiff_prim*area(1)/(dm_minus(2) + dm_plus(1))
          coeff_right(j,k,i) = coeff_right(j,k,1) - &
                              pordiff_prim*area(1)/(dm_minus(2) + dm_plus(1))                               
          
        ! Outer boundary -- closest to primary node
        else !if (i.eq.ngcells) then
          coeff_diag(j,k,ngcells) = coeff_diag(j,k,ngcells) + &
                                    pordiff_prim*area(ngcells-1)/(dm_minus(ngcells) &
                                    + dm_plus(ngcells-1)) &
                                    + pordiff_prim*area(ngcells)/dm_plus(ngcells)
          coeff_left(j,k,ngcells) = coeff_left(j,k,ngcells) - &
                                    pordiff_prim*area(ngcells-1)/(dm_minus(ngcells) &
                                    + dm_plus(ngcells-1))
        endif
      enddo
    enddo
  
  
    ! d (sumatory ( frac * D_i * C_i )) / dC_j =
    ! sumatory ( frac * D_i * dC_i/dC_j )
    do icplx = 1, reaction%neqcplx

      lnQK = -reaction%eqcplx_logK(icplx)*LOG_TO_LN

      ! activity of water
      if (reaction%eqcplxh2oid(icplx) > 0) then
        lnQK = lnQK + reaction%eqcplxh2ostoich(icplx)*rt_auxvar%ln_act_h2o
      endif

      nicomp = reaction%eqcplxspecid(0,icplx)
      do ni = 1, nicomp
        icomp = reaction%eqcplxspecid(ni,icplx)
        lnQK = lnQK + reaction%eqcplxstoich(ni,icplx)*ln_act(icomp)
      enddo

      pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)* &
                    tortuosity*global_auxvar%den_kg(1)

      nicomp = reaction%eqcplxspecid(0,icplx)
      do j = 1, nicomp
        jcomp = reaction%eqcplxspecid(j,icplx)        
        tempreal = reaction%eqcplxstoich(j,icplx)*exp(lnQK-ln_conc(jcomp))/ &
                                                 rt_auxvar%sec_act_coef(icplx)

        do k = 1, nicomp
          icomp = reaction%eqcplxspecid(k,icplx)
          
          ! Flux terms
          if (i.gt.1.and.i.lt.ngcells) then
            
            J_up = pordiff_sec*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                   reaction%eqcplxstoich(k,icplx)*tempreal
            J_dn = pordiff_sec*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                   reaction%eqcplxstoich(k,icplx)*tempreal
            
            coeff_diag(icomp,jcomp,i) = coeff_diag(icomp,jcomp,i) + &
                                        J_up + J_dn
            coeff_left(icomp,jcomp,i) = coeff_left(icomp,jcomp,i) - J_dn
            coeff_right(icomp,jcomp,i) = coeff_right(icomp,jcomp,i) - J_up

          ! Apply boundary conditions
          ! Inner boundary
          else if (i.eq.1) then
            
            J_up = pordiff_sec*area(1)/(dm_minus(2) + dm_plus(1))* &
                   reaction%eqcplxstoich(k,icplx)*tempreal
            
            coeff_diag(j,k,1) = coeff_diag(j,k,1) + J_up 
            coeff_right(j,k,1) = coeff_right(j,k,1) - J_up
        
          ! Outer boundary -- closest to primary node
          else !if (i.eq.ngcells) then
            
            J_up = pordiff_sec*area(ngcells)/dm_plus(ngcells)* &
                   reaction%eqcplxstoich(k,icplx)*tempreal
            J_dn = pordiff_sec*area(ngcells-1)/(dm_minus(ngcells) + dm_plus(ngcells-1))* &
                   reaction%eqcplxstoich(k,icplx)*tempreal
            
            coeff_diag(j,k,ngcells) = coeff_diag(j,k,ngcells) + &
                                      J_up + J_dn
            coeff_left(j,k,ngcells) = coeff_left(j,k,ngcells) - &
                                      J_dn
          endif

        enddo
      enddo
    enddo
  enddo

  
!====================== Add reaction contributions =============================        
  
  ! Reaction 
  do i = 1, ngcells
    res_react = 0.d0
    jac_react = 0.d0
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i), &
                      option)
    rt_auxvar%pri_molal = conc_upd(:,i) ! in mol/kg
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    material_auxvar%volume = vol(i)
    call RReaction(res_react,jac_react,PETSC_TRUE, &
                   rt_auxvar,global_auxvar,material_auxvar,reaction,option)                     
    do j = 1, ncomp
      res(j+(i-1)*ncomp) = res(j+(i-1)*ncomp) + res_react(j) 
    enddo
    coeff_diag(:,:,i) = coeff_diag(:,:,i) + jac_react  ! in kg water/s
  enddo  
  call MaterialAuxVarStrip(material_auxvar)
  deallocate(material_auxvar)
         
!============================== Forward solve ==================================        
                        
  rhs = -res   
  
  if (reaction%use_log_formulation) then
  ! scale the jacobian by concentrations
    i = 1
    do k = 1, ncomp
      coeff_diag(:,k,i) = coeff_diag(:,k,i)*conc_upd(k,i) ! m3/s*kg/L
      coeff_right(:,k,i) = coeff_right(:,k,i)*conc_upd(k,i+1)
    enddo
    do i = 2, ngcells-1
      do k = 1, ncomp
        coeff_diag(:,k,i) = coeff_diag(:,k,i)*conc_upd(k,i) ! m3/s*kg/L
        coeff_left(:,k,i) = coeff_left(:,k,i)*conc_upd(k,i-1)
        coeff_right(:,k,i) = coeff_right(:,k,i)*conc_upd(k,i+1)
      enddo
    enddo
    i = ngcells
      do k = 1, ncomp
        coeff_diag(:,k,i) = coeff_diag(:,k,i)*conc_upd(k,i) ! m3/s*kg/L
        coeff_left(:,k,i) = coeff_left(:,k,i)*conc_upd(k,i-1)
      enddo
  endif 
  
  ! First do an LU decomposition for calculating D_M matrix
  coeff_diag_dm = coeff_diag
  coeff_left_dm = coeff_left
  coeff_right_dm = coeff_right
  
  select case (option%secondary_continuum_solver)
    case(1) 
      do i = 2, ngcells
        coeff_left_dm(:,:,i-1) = coeff_left_dm(:,:,i)
      enddo
      coeff_left_dm(:,:,ngcells) = 0.d0
      call bl3dfac(ngcells,ncomp,coeff_right_dm,coeff_diag_dm,coeff_left_dm,pivot)  
    case(2)
      call decbt(ncomp,ngcells,ncomp,coeff_diag_dm,coeff_right_dm,coeff_left_dm,pivot,ier)
      if (ier /= 0) then
        print *,'error in matrix decbt: ier = ',ier
        stop
      endif
    case(3)
      ! Thomas algorithm for tridiagonal system
      ! Forward elimination
      if (ncomp /= 1) then
        option%io_buffer = 'THOMAS algorithm can be used only with single '// &
                           'component chemistry'
        call PrintErrMsg(option)
      endif
      do i = 2, ngcells
        m = coeff_left_dm(ncomp,ncomp,i)/coeff_diag_dm(ncomp,ncomp,i-1)
        coeff_diag_dm(ncomp,ncomp,i) = coeff_diag_dm(ncomp,ncomp,i) - &
                                    m*coeff_right_dm(ncomp,ncomp,i-1)
      enddo        
    case default
      option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                         'HINDMARSH or KEARST. For single component'// &
                         'chemistry THOMAS can be used.'
      call PrintErrMsg(option)
  end select
  
  ! Set the values of D_M matrix and create identity matrix of size ncomp x ncomp  
  do i = 1, ncomp
    do j = 1, ncomp
      D_M(i,j) = coeff_diag_dm(i,j,ngcells)
      if (j == i) then
        identity(i,j) = 1.d0
      else
        identity(i,j) = 0.d0
      endif
    enddo
  enddo
  
  ! Find the inverse of D_M
  call LUDecomposition(D_M,ncomp,indx,d) 
  do j = 1, ncomp
    call LUBackSubstitution(D_M,ncomp,indx,identity(1,j))
  enddo  
  inv_D_M = identity      
  
  if (option%numerical_derivatives_multi_coupling) then  
    ! Store the coeffs for numerical jacobian
    coeff_diag_copy = coeff_diag
    coeff_left_copy = coeff_left
    coeff_right_copy = coeff_right
  endif

  select case (option%secondary_continuum_solver)
    case(1) 
      do i = 2, ngcells
        coeff_left(:,:,i-1) = coeff_left(:,:,i)
      enddo
      coeff_left(:,:,ngcells) = 0.d0
      call bl3dfac(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot)  
      call bl3dsolf(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot, &
                    ONE_INTEGER,rhs)
    case(2)
      call decbt(ncomp,ngcells,ncomp,coeff_diag,coeff_right,coeff_left, &
                 pivot,ier)
      if (ier /= 0) then
        print *,'error in matrix decbt: ier = ',ier
        stop
      endif
      call solbtf(ncomp,ngcells,ncomp,coeff_diag,coeff_right,coeff_left, &
                  pivot,rhs)
    case(3)
      ! Thomas algorithm for tridiagonal system
      ! Forward elimination
      if (ncomp /= 1) then
        option%io_buffer = 'THOMAS algorithm can be used only with single '// &
                           'component chemistry'
        call PrintErrMsg(option)
      endif
      do i = 2, ngcells
        m = coeff_left(ncomp,ncomp,i)/coeff_diag(ncomp,ncomp,i-1)
        coeff_diag(ncomp,ncomp,i) = coeff_diag(ncomp,ncomp,i) - &
                                    m*coeff_right(ncomp,ncomp,i-1)
        rhs(i) = rhs(i) - m*rhs(i-1)
      enddo        
      rhs(ngcells) = rhs(ngcells)/coeff_diag(ncomp,ncomp,ngcells)
    case default
      option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                         'HINDMARSH or KEARST. For single component'// &
                         'chemistry THOMAS can be used.'
      call PrintErrMsg(option)
  end select
    
  ! Update the secondary concentrations
  do i = 1, ncomp
    if (reaction%use_log_formulation) then
      ! convert log concentration to concentration
      rhs(i+(ngcells-1)*ncomp) = dsign(1.d0,rhs(i+(ngcells-1)*ncomp))* &
        min(dabs(rhs(i+(ngcells-1)*ncomp)),reaction%max_dlnC)
      conc_current_M(i) = conc_upd(i,ngcells)*exp(rhs(i+(ngcells-1)*ncomp))
    else
      conc_current_M(i) = conc_upd(i,ngcells) + rhs(i+(ngcells-1)*ncomp)
    endif
  enddo

  ! Update the secondary continuum totals at the outer matrix node
  call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(ngcells), &
                    option)
  rt_auxvar%pri_molal = conc_current_M ! in mol/kg
  call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
  total_current_M = rt_auxvar%total(:,1)
  if (reaction%neqcplx > 0) sec_sec_molal_M = rt_auxvar%sec_molal
  call RTAuxVarStrip(rt_auxvar)
  

  
  b_m = porosity*tortuosity/dm_plus(ngcells)*area(ngcells)*inv_D_M ! in m·s/mol
   
  dCsec_dCprim = b_m*dtotal_prim

  ! Calculate the dervative of outer matrix node total with respect to the 
  ! primary node concentration
  if (reaction%use_log_formulation) then ! log formulation
    do j = 1, ncomp
      do l = 1, ncomp
        dPsisec_dCprim(j,l) = dCsec_dCprim(j,l)*conc_current_M(j)
      enddo
    enddo
   
    if (reaction%neqcplx > 0) then
      do icplx = 1, reaction%neqcplx
        ncompeq = reaction%eqcplxspecid(0,icplx)
        do j = 1, ncompeq
          jcomp = reaction%eqcplxspecid(j,icplx)
          do l = 1, ncompeq
            lcomp = reaction%eqcplxspecid(l,icplx)
            do k = 1, ncompeq
              kcomp = reaction%eqcplxspecid(k,icplx)
              dPsisec_dCprim(jcomp,lcomp) = dPsisec_dCprim(jcomp,lcomp) + &
                                            reaction%eqcplxstoich(j,icplx)* &
                                            reaction%eqcplxstoich(k,icplx)* &
                                            dCsec_dCprim(kcomp,lcomp)* &
                                            sec_sec_molal_M(icplx)                                
            enddo
          enddo      
        enddo
      enddo
    endif
   
  else   ! linear case  

    dPsisec_dCprim = dCsec_dCprim 
    
    if (reaction%neqcplx > 0) then
      do icplx = 1, reaction%neqcplx
        ncompeq = reaction%eqcplxspecid(0,icplx)
        do j = 1, ncompeq
          jcomp = reaction%eqcplxspecid(j,icplx)
          do l = 1, ncompeq
            lcomp = reaction%eqcplxspecid(l,icplx)
            do k = 1, ncompeq
              kcomp = reaction%eqcplxspecid(k,icplx)

              dPsisec_dCprim(jcomp,lcomp) = dPsisec_dCprim(jcomp,lcomp) + &
                                            reaction%eqcplxstoich(j,icplx)* &
                                            reaction%eqcplxstoich(k,icplx)* &
                                            dCsec_dCprim(kcomp,lcomp)* &
                                            sec_sec_molal_M(icplx)/ &
                                            conc_current_M(kcomp)
            enddo
          enddo      
        enddo
      enddo
    endif
  
  endif
  
  dPsisec_dCprim = dPsisec_dCprim*global_auxvar%den_kg(1)*1.d-3 ! in kg/L 
  
  ! D_j \nabla C_j
  do j = 1, ncomp    
    pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(j)* &
                   tortuosity*global_auxvar%den_kg(1)
    res_transport(j) = pordiff_prim*area_fm/dm_plus(ngcells)* &
                       (conc_current_M(j) - &
                       auxvar%pri_molal(j))*prim_vol ! in mol/s

  enddo
  
  ! summatory D_i \nabla C_i 
  do icplx = 1, reaction%neqcplx
    pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)* &
                  tortuosity*global_auxvar%den_kg(1)
    nicomp = reaction%eqcplxspecid(0,icplx)
    do ni = 1, nicomp
      icomp = reaction%eqcplxspecid(ni,icplx)
      res_transport(icomp) = res_transport(icomp) + &
                       pordiff_sec*area_fm/dm_plus(ngcells)* &
                       (sec_sec_molal_M(icplx) - &
                       auxvar%sec_molal(icplx))*prim_vol* &
                       reaction%eqcplxstoich(ni,icplx) ! in mol/s
    enddo
  enddo
  
  ! Electromigration term
  call ComputeElectricPotentialTotalComponent_NP(ngcells,reaction, &
                                          rt_parameter,  &
                                          rt_auxvar, &
                                          sec_transport_vars, &
                                          icomp,auxvar, &
                                          potent_dn, &
                                          potent_up, &
                                          mc_pri_molal, &
                                          mc_sec_molal)

  sum_denom_up = 1d-40
    
  do jcomp = 1, rt_parameter%naqcomp
    curCharge = reaction%primary_spec_Z(jcomp)
    sum_denom_up = sum_denom_up + potent_up(jcomp)*curCharge
  enddo

  do jcomp = 1, rt_parameter%naqcomp
    factor_up_em(jcomp) = potent_up(jcomp) / sum_denom_up
  enddo

  portort = porosity*tortuosity*global_auxvar%den_kg(1)

  sum_transp_up = 0d0
    
  ! Primary l index
  do lcomp = 1, rt_parameter%naqcomp
    curCharge = reaction%primary_spec_Z(lcomp)
    cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)    
    curConc_up = area_fm/dm_plus(ngcells)*prim_vol* &
              (auxvar%pri_molal(lcomp)-conc_current_M(lcomp))             
    sum_transp_up = sum_transp_up - curCharge*cur_diff*curConc_up

    ! Secondary term k
    do icplx = 1, reaction%neqcplx
      nicomp = reaction%eqcplxspecid(0,icplx)
      do ni = 1, nicomp
        kcomp = reaction%eqcplxspecid(ni,icplx)
        if (kcomp == lcomp) then                    
          curCharge = reaction%primary_spec_Z(kcomp)
          cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
          curStoich = reaction%eqcplxstoich(ni,icplx)          
          curConc_up = area_fm/dm_plus(ngcells)*prim_vol* &
                    (auxvar%sec_molal(icplx) - sec_sec_molal_M(icplx))   
          sum_transp_up = sum_transp_up - curCharge*cur_diff*curStoich*curConc_up
        endif

      enddo
    enddo
  enddo
  
  do icomp = 1, rt_parameter%naqcomp
    res_transport(icomp) = res_transport(icomp) - portort*factor_up_em(icomp)*sum_transp_up ! in mol/s
  enddo

  ! Calculate the jacobian contribution due to coupling term
  sec_jac = area_fm*pordiff/dm_plus(ngcells)*(dPsisec_dCprim - dtotal_prim)* &
            prim_vol !*1.d3 ! in kg water/s

  ! Store the contribution to the primary jacobian term
  sec_transport_vars%sec_jac = sec_jac 
  sec_transport_vars%sec_jac_update = PETSC_TRUE
  
  ! Store the coefficients from LU decomposition of the block tridiagonal
  ! sytem. These will be called later to perform backsolve to the get the
  ! updated secondary continuum concentrations at the end of the timestep
  sec_transport_vars%cxm = coeff_left
  sec_transport_vars%cxp = coeff_right
  sec_transport_vars%cdl = coeff_diag
  
  ! Store the solution of the forward solve
  sec_transport_vars%r = rhs
  

!============== Numerical jacobian for coupling term ===========================


  if (option%numerical_derivatives_multi_coupling) then

    call RTAuxVarInit(rt_auxvar,reaction,option)
    conc_prim = auxvar%pri_molal
    conc_prim_pert = conc_prim
  
    do l = 1, ncomp
      
      conc_prim_pert = conc_prim
      pert = conc_prim(l)*perturbation_tolerance
      conc_prim_pert(l) = conc_prim_pert(l) + pert
  
      res = 0.d0
      rhs = 0.d0
    
      coeff_diag_pert = coeff_diag_copy
      coeff_left_pert = coeff_left_copy
      coeff_right_pert = coeff_right_copy

      call RTAuxVarCopy(rt_auxvar,auxvar,option)
      rt_auxvar%pri_molal = conc_prim_pert ! in mol/kg
      call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
      total_primary_node_pert = rt_auxvar%total(:,1)
                          
!================ Calculate the secondary residual =============================        
      do i = 1, ngcells
  
        ! D_j \nabla C_j
        do icomp = 1, ncomp
      
          pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(icomp)*&
        	             tortuosity*global_auxvar%den_kg(1)
        	             
          n = icomp + (i-1)*ncomp
          
          ! Accumulation        
          res(n) = pordt*(total_upd(icomp,i) - total_prev(icomp,i))*vol(i)    ! in mol/L*m3/s
          if (reaction%neqsorb > 0) then 
            res(n) = res(n) + vol(i)/option%tran_dt*(total_sorb_upd(icomp,i) - total_sorb_prev(icomp,i))
          endif      
    
          ! Flux terms
          if (i.gt.1.and.i.lt.ngcells) then
            res(n) = res(n) - pordiff_prim*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                            (sec_transport_vars%sec_rt_auxvar(i+1)%pri_molal(icomp) - &
                            sec_transport_vars%sec_rt_auxvar(i)%pri_molal(icomp))
            res(n) = res(n) + pordiff_prim*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                            (sec_transport_vars%sec_rt_auxvar(i)%pri_molal(icomp) - &
                            sec_transport_vars%sec_rt_auxvar(i-1)%pri_molal(icomp)) 
    
          ! Apply boundary conditions
          ! Inner boundary
          else if (i.eq.1) then
            res(n) = res(n) - pordiff_prim*area(1)/(dm_minus(2) + dm_plus(1))* &
                            (sec_transport_vars%sec_rt_auxvar(2)%pri_molal(icomp) - &
                            sec_transport_vars%sec_rt_auxvar(1)%pri_molal(icomp))
               
          ! Outer boundary -- closest to primary node
          else !if (i.eq.ngcells) then
            res(n) = res(n) - &
                     pordiff_prim*area(ngcells)/dm_plus(ngcells)* &
                     (auxvar%pri_molal(icomp) - &
                     sec_transport_vars%sec_rt_auxvar(ngcells)%pri_molal(icomp))
            res(n) = res(n) + &
                     pordiff_prim*area(ngcells-1)/(dm_minus(ngcells) &
                     + dm_plus(ngcells-1))* &
                     (sec_transport_vars%sec_rt_auxvar(ngcells)%pri_molal(icomp) - &
                     sec_transport_vars%sec_rt_auxvar(ngcells-1)%pri_molal(icomp))
          endif  
    
        enddo
    
    
        ! summatory D_i \nabla C_i   
        do icplx = 1, reaction%neqcplx
        
          pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)*&
        	       tortuosity*global_auxvar%den_kg(1)
          nicomp = reaction%eqcplxspecid(0,icplx)
          
          do ni = 1, nicomp
            icomp = reaction%eqcplxspecid(ni,icplx)
            
            n = icomp + (i-1)*ncomp

            ! Flux terms
            if (i.gt.1.and.i.lt.ngcells) then
              res(n) = res(n) - pordiff_sec*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                       (sec_transport_vars%sec_rt_auxvar(i+1)%sec_molal(icplx) - &
                       sec_transport_vars%sec_rt_auxvar(i)%sec_molal(icplx))* &
                       reaction%eqcplxstoich(ni,icplx)
              res(n) = res(n) + pordiff_sec*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                       (sec_transport_vars%sec_rt_auxvar(i)%sec_molal(icplx) - &
                       sec_transport_vars%sec_rt_auxvar(i-1)%sec_molal(icplx))* &
                       reaction%eqcplxstoich(ni,icplx)
        
            ! Apply boundary conditions
            ! Inner boundary
            else if (i.eq.1) then
              res(n) = res(n) - pordiff_sec*area(1)/(dm_minus(2) + dm_plus(1))* &
                       (sec_transport_vars%sec_rt_auxvar(2)%sec_molal(icplx) - &
                       sec_transport_vars%sec_rt_auxvar(1)%sec_molal(icplx))* &
                       reaction%eqcplxstoich(ni,icplx)
        
            ! Outer boundary -- closest to primary node
            else !if (i.eq.ngcells) then
              res(n) = res(n) - &
                       pordiff_sec*area(ngcells)/dm_plus(ngcells)* &
                       (auxvar%sec_molal(icplx) - &
                       sec_transport_vars%sec_rt_auxvar(ngcells)%sec_molal(icplx))* & 
                       reaction%eqcplxstoich(ni,icplx)                                 
              res(n) = res(n) + pordiff_sec* &
                       area(ngcells-1)/(dm_minus(ngcells) + dm_plus(ngcells-1))* &
                       (sec_transport_vars%sec_rt_auxvar(ngcells)%sec_molal(icplx) - &
                       sec_transport_vars%sec_rt_auxvar(ngcells-1)%sec_molal(icplx))* &
                       reaction%eqcplxstoich(ni,icplx)
            endif    
          enddo
        enddo    
    
        ! Electromigration term
        call ComputeElectricPotentialTotalComponent_NP(i,reaction, &
                                                rt_parameter,  &
                                                rt_auxvar, &
                                                sec_transport_vars, &
                                                icomp,auxvar, &
                                                potent_dn, &
                                                potent_up, &
                                                mc_pri_molal, &
                                                mc_sec_molal)

        sum_denom_up = 1d-40
        sum_denom_dn = 1d-40
    
        do jcomp = 1, rt_parameter%naqcomp
          curCharge = reaction%primary_spec_Z(jcomp)
          sum_denom_dn = sum_denom_dn + potent_dn(jcomp)*curCharge
          sum_denom_up = sum_denom_up + potent_up(jcomp)*curCharge
        enddo

        do jcomp = 1, rt_parameter%naqcomp
          factor_up_em(jcomp) = potent_up(jcomp) / sum_denom_up
          factor_dn_em(jcomp) = potent_dn(jcomp) / sum_denom_dn
        enddo
    
        portort = porosity*tortuosity* &
        global_auxvar%den_kg(1)        
    
        ! Primary l index
        do lcomp = 1, rt_parameter%naqcomp
          curCharge = reaction%primary_spec_Z(lcomp)
          cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)
          
          ! Flux terms
          if (i.gt.1.and.i.lt.ngcells) then
            curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                         (sec_transport_vars%sec_rt_auxvar(i+1)%pri_molal(lcomp) - &
                         sec_transport_vars%sec_rt_auxvar(i)%pri_molal(lcomp))
            curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                         (sec_transport_vars%sec_rt_auxvar(i)%pri_molal(lcomp) - &
                         sec_transport_vars%sec_rt_auxvar(i-1)%pri_molal(lcomp))      
          
          ! Apply boundary conditions
          ! Inner boundary
          else if (i.eq.1) then
            curConc_up = area(1)/(dm_minus(2) + dm_plus(1))* &
                         (sec_transport_vars%sec_rt_auxvar(2)%pri_molal(lcomp) - &
                         sec_transport_vars%sec_rt_auxvar(1)%pri_molal(lcomp))
            curConc_dn = 0.d0
          
          ! Outer boundary    
          else !if (i.eq.ngcells) then
            curConc_up = area(ngcells)/dm_plus(ngcells)* &
                         (auxvar%pri_molal(lcomp) - &
                         sec_transport_vars%sec_rt_auxvar(ngcells)%pri_molal(lcomp))
            curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                         + dm_plus(ngcells-1))* &
                         (sec_transport_vars%sec_rt_auxvar(ngcells)%pri_molal(lcomp) - &
                         sec_transport_vars%sec_rt_auxvar(ngcells-1)%pri_molal(lcomp))
          endif
      
          sum_transp_up = sum_transp_up - curCharge*cur_diff*(curConc_up)
          sum_transp_dn = sum_transp_dn + curCharge*cur_diff*(curConc_dn)

          ! Secondary term k
          do icplx = 1, reaction%neqcplx
            nicomp = reaction%eqcplxspecid(0,icplx)
            do ni = 1, nicomp
              kcomp = reaction%eqcplxspecid(ni,icplx)
              if (kcomp == lcomp) then  
                curCharge = reaction%primary_spec_Z(kcomp)
                cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
                curStoich = reaction%eqcplxstoich(ni,icplx)
            
                ! Flux terms
                if (i.gt.1.and.i.lt.ngcells) then
                  curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                              (sec_transport_vars%sec_rt_auxvar(i+1)%sec_molal(icplx) - &
                              sec_transport_vars%sec_rt_auxvar(i)%sec_molal(icplx))
                  curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                              (sec_transport_vars%sec_rt_auxvar(i)%sec_molal(icplx) - &
                              sec_transport_vars%sec_rt_auxvar(i-1)%sec_molal(icplx))
                
                ! Apply boundary conditions
                ! Inner boundary
                else if (i.eq.1) then
                  curConc_up = area(1)/(dm_minus(2) + dm_plus(1))* &
                              (sec_transport_vars%sec_rt_auxvar(2)%sec_molal(icplx) - &
                              sec_transport_vars%sec_rt_auxvar(1)%sec_molal(icplx))
                  curConc_dn = 0.d0
              
                ! Outer boundary   
                else !if (i.eq.ngcells) then 
                  curConc_up = area(ngcells)/dm_plus(ngcells)* &
                              (auxvar%sec_molal(icplx) - &
                              sec_transport_vars%sec_rt_auxvar(ngcells)%sec_molal(icplx))
                  curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                              + dm_plus(ngcells-1))* &
                              (sec_transport_vars%sec_rt_auxvar(ngcells)%sec_molal(icplx) - &
                              sec_transport_vars%sec_rt_auxvar(ngcells-1)%sec_molal(icplx))
                endif
            
                sum_transp_up = sum_transp_up - curCharge*cur_diff*curStoich*(curConc_up)
                sum_transp_dn = sum_transp_dn + curCharge*cur_diff*curStoich*(curConc_dn)
              endif
            enddo
          enddo  
        enddo

        do icomp = 1, rt_parameter%naqcomp
          n = icomp + (i-1)*ncomp
          res(n) = res(n) - portort*factor_up_em(icomp)*sum_transp_up
          res(n) = res(n) - portort*factor_dn_em(icomp)*sum_transp_dn
        enddo
   
      enddo
                                                                
!============================== Forward solve ==================================        
                        
    rhs = -res   
           
    select case (option%secondary_continuum_solver)
      case(1) 
        call bl3dfac(ngcells,ncomp,coeff_right_pert,coeff_diag_pert, &
                      coeff_left_pert,pivot)  
        call bl3dsolf(ngcells,ncomp,coeff_right_pert,coeff_diag_pert, &
                       coeff_left_pert,pivot,ONE_INTEGER,rhs)
      case(2)
        call decbt(ncomp,ngcells,ncomp,coeff_diag_pert,coeff_right_pert, &
                    coeff_left_pert,pivot,ier)
        if (ier /= 0) then
          print *,'error in matrix decbt: ier = ',ier
          stop
        endif
        call solbtf(ncomp,ngcells,ncomp,coeff_diag_pert,coeff_right_pert, &
                     coeff_left_pert,pivot,rhs)
      case(3)
        ! Thomas algorithm for tridiagonal system
        ! Forward elimination
        if (ncomp /= 1) then
          option%io_buffer = 'THOMAS algorithm can be used only with '// &
                             'single component chemistry'
          call PrintErrMsg(option)
        endif
        do i = 2, ngcells
          m = coeff_left_pert(ncomp,ncomp,i)/coeff_diag_pert(ncomp,ncomp,i-1)
          coeff_diag_pert(ncomp,ncomp,i) = coeff_diag_pert(ncomp,ncomp,i) - &
                                      m*coeff_right_pert(ncomp,ncomp,i-1)
          rhs(i) = rhs(i) - m*rhs(i-1)
        enddo        
        rhs(ngcells) = rhs(ngcells)/coeff_diag(ncomp,ncomp,ngcells)
      case default
        option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                           'HINDMARSH or KEARST. For single component'// &
                           'chemistry THOMAS can be used.'
        call PrintErrMsg(option)
      end select      
    
      ! Update the secondary concentrations
      do i = 1, ncomp
        if (reaction%use_log_formulation) then
          ! convert log concentration to concentration
          rhs(i+(ngcells-1)*ncomp) = dsign(1.d0,rhs(i+(ngcells-1)*ncomp))* &
            min(dabs(rhs(i+(ngcells-1)*ncomp)),reaction%max_dlnC)
          conc_current_M_pert(i) = conc_upd(i,ngcells)* &
                                     exp(rhs(i+(ngcells-1)*ncomp))
        else
          conc_current_M_pert(i) = conc_upd(i,ngcells) + &
                                     rhs(i+(ngcells-1)*ncomp)
        endif
      enddo

      ! Update the secondary continuum totals at the outer matrix node
      call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(ngcells), &
                        option)
      rt_auxvar%pri_molal = conc_current_M_pert ! in mol/kg
      call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
      total_current_M_pert = rt_auxvar%total(:,1)
             
      ! ! Calculate the coupling term
      ! D_j \nabla C_j
      do j = 1, ncomp
    
        pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(j)* &
                   tortuosity*global_auxvar%den_kg(1)

        res_transport_pert(j) = res_transport_pert(j) + &
                       pordiff_prim*area_fm/dm_plus(ngcells)* &
                       (conc_current_M_pert(j) - &
                       auxvar%pri_molal(j))*prim_vol ! in mol/s

      enddo
  
      ! summatory D_i \nabla C_i 
      do icplx = 1, reaction%neqcplx
    
        pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)* &
                  tortuosity*global_auxvar%den_kg(1)
  
        nicomp = reaction%eqcplxspecid(0,icplx)
        do ni = 1, nicomp
          icomp = reaction%eqcplxspecid(ni,icplx)
      
          res_transport_pert(icomp) = res_transport_pert(icomp) + &
                       pordiff_sec*area_fm/dm_plus(ngcells)* &
                       (rt_auxvar%sec_molal(icplx) - &
                       auxvar%sec_molal(icplx))*prim_vol* &
                       reaction%eqcplxstoich(ni,icplx) ! in mol/s                       
        enddo
      enddo
  
      ! Electromigration term
      call ComputeElectricPotentialTotalComponent_NP(ngcells,reaction, &
                                          rt_parameter,  &
                                          rt_auxvar, &
                                          sec_transport_vars, &
                                          icomp,auxvar, &
                                          potent_dn, &
                                          potent_up, &
                                          mc_pri_molal, &
                                          mc_sec_molal)

      sum_denom_up = 1d-40
    
      do jcomp = 1, rt_parameter%naqcomp
        curCharge = reaction%primary_spec_Z(jcomp)
        sum_denom_up = sum_denom_up + potent_up(jcomp)*curCharge
      enddo

      do jcomp = 1, rt_parameter%naqcomp
        factor_up_em(jcomp) = potent_up(jcomp) / sum_denom_up
      enddo

      portort = porosity*tortuosity*global_auxvar%den_kg(1)
    
      sum_transp_up = 0d0
    
      ! Primary l index
      do lcomp = 1, rt_parameter%naqcomp
        curCharge = reaction%primary_spec_Z(lcomp)
        cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)
    
        curConc_up = area_fm/dm_plus(ngcells)*prim_vol* &
              (auxvar%pri_molal(lcomp) - conc_current_M_pert(lcomp))
        sum_transp_up = sum_transp_up - curCharge*cur_diff*curConc_up

        ! Secondary term k
        do icplx = 1, reaction%neqcplx
          nicomp = reaction%eqcplxspecid(0,icplx)
          do ni = 1, nicomp
            kcomp = reaction%eqcplxspecid(ni,icplx)
            if (kcomp == lcomp) then  
              curCharge = reaction%primary_spec_Z(kcomp)
              cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
              curStoich = reaction%eqcplxstoich(ni,icplx)
          
              curConc_up = area_fm/dm_plus(ngcells)*prim_vol* &
                  (auxvar%sec_molal(icplx) - rt_auxvar%sec_molal(icplx))   
              sum_transp_up = sum_transp_up - curCharge*cur_diff*curStoich*curConc_up
            endif
          enddo
        enddo
      enddo
  
      do icomp = 1, rt_parameter%naqcomp
        res_transport_pert(icomp) = res_transport_pert(icomp) - &
                                  portort*factor_up_em(icomp)*sum_transp_up ! in mol/s
      enddo
  
      dtotal_prim_num(:,l) = (total_primary_node_pert(:) - &
                               total_primary_node(:))/pert
  
      dPsisec_dCprim_num(:,l) = (total_current_M_pert(:) - &
                                  total_current_M(:))/pert
  
      sec_jac_num(:,l) = (res_transport_pert(:) - res_transport(:))/pert
    
    enddo    

    call RTAuxVarStrip(rt_auxvar)
    sec_transport_vars%sec_jac = sec_jac_num 

  endif
  

end subroutine SecondaryRTResJacMulti_NP

! ************************************************************************** !

subroutine SecondaryRTUpdateIterate_NP(snes,P0,dP,P1,dX_changed, &
                                    X1_changed,realization,ierr)
  ! 
  ! Checks update after the update is done
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 02/22/13
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Grid_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
 
  implicit none
  
  SNES :: snes
  Vec :: P0
  Vec :: dP
  Vec :: P1
  class(realization_subsurface_type) :: realization
  ! ignore changed flag for now.
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  
  type(reactive_transport_param_type), pointer :: rt_parameter
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(sec_transport_type), pointer :: rt_sec_transport_vars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(reaction_rt_type), pointer :: reaction
  PetscInt :: local_id, ghosted_id
  PetscReal :: sec_diffusion_coefficient
  PetscReal :: sec_porosity
  PetscReal :: sec_tortuosity
  
  PetscErrorCode :: ierr
  PetscReal :: inf_norm_sec
  PetscReal :: max_inf_norm_sec
  
  option => realization%option
  grid => realization%patch%grid
  rt_auxvars => realization%patch%aux%RT%auxvars
  global_auxvars => realization%patch%aux%Global%auxvars
  reaction => realization%reaction
  rt_parameter => realization%patch%aux%RT%rt_parameter
  if (option%use_sc) then
    rt_sec_transport_vars => realization%patch%aux%SC_RT%sec_transport_vars
  endif  
  
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  max_inf_norm_sec = 0.d0
  
  if (option%use_sc) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (realization%patch%imat(ghosted_id) <= 0) cycle
      sec_diffusion_coefficient = realization%patch% &
                                  material_property_array(1)%ptr% &
                                  multicontinuum%diff_coeff(1)
      sec_porosity = realization%patch%material_property_array(1)%ptr% &
                    multicontinuum%porosity

      sec_tortuosity = realization%patch%material_property_array(1)%ptr% &
                    multicontinuum%tortuosity                    

      call SecondaryRTAuxVarComputeMulti(&
                                    rt_sec_transport_vars(ghosted_id), &
                                    reaction, &
                                    option)              
 
      call SecondaryRTCheckResidual_np(rt_sec_transport_vars(ghosted_id), &
                                    rt_auxvars(ghosted_id), &
                                    global_auxvars(ghosted_id), &
                                    reaction,rt_parameter,sec_diffusion_coefficient, &
                                    sec_porosity, sec_tortuosity, &
                                    option,inf_norm_sec)
                                      
      max_inf_norm_sec = max(max_inf_norm_sec,inf_norm_sec)                                                                   
    enddo 
    call MPI_Allreduce(max_inf_norm_sec,option%infnorm_res_sec,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
  endif
  
      
end subroutine SecondaryRTUpdateIterate_NP

! ************************************************************************** !

subroutine SecondaryRTUpdateEquilState(sec_transport_vars,global_auxvars, &
                                       reaction,sec_porosity,option) 
  ! 
  ! Updates the equilibrium secondary continuum
  ! variables
  ! at the end of time step
  ! 
  ! Author: Satish Karra, LANL; Glenn Hammond (modification)
  ! Date: 02/22/13; 06/27/13
  ! 
                                     

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Reaction_module
  use Global_Aux_module
  use Material_Aux_module
 
  implicit none
  
  class(material_auxvar_type), allocatable :: material_auxvar
  type(option_type), pointer :: option
  type(sec_transport_type) :: sec_transport_vars
  type(global_auxvar_type) :: global_auxvars
  class(reaction_rt_type), pointer :: reaction
  PetscReal :: sec_porosity
  PetscInt :: ngcells,ncomp
  PetscInt :: i,j
  
  ngcells = sec_transport_vars%ncells
  ncomp = reaction%naqcomp

  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  material_auxvar%porosity = sec_porosity

  
  do j = 1, ncomp
    do i = 1, ngcells
      sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j) = sec_transport_vars%&
        updated_conc(j,i)
    enddo
  enddo
    
  do i = 1, ngcells
    call RTotalAqueous(sec_transport_vars%sec_rt_auxvar(i),global_auxvars, &
                       reaction,option)
    if (reaction%neqsorb > 0) then
      call RTotalSorb(sec_transport_vars%sec_rt_auxvar(i),global_auxvars, &
                      material_auxvar,reaction, &
                      reaction%isotherm%multicontinuum_isotherm_rxn,option)
    endif
  enddo
 
end subroutine SecondaryRTUpdateEquilState

! ************************************************************************** !

subroutine SecondaryRTUpdateKineticState(sec_transport_vars,global_auxvars, &
                                         reaction,porosity,option) 
  ! 
  ! Updates the kinetic secondary continuum
  ! variables at the end of time step
  ! 
  ! Author: Satish Karra, LANL; Glenn Hammond (modification)
  ! Date: 02/22/13; 06/27/13
  ! 
                                     

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Reaction_module
  use Global_Aux_module
  use Material_Aux_module
 
  implicit none
  

  type(option_type), pointer :: option
  type(sec_transport_type) :: sec_transport_vars
  type(global_auxvar_type) :: global_auxvars
  class(reaction_rt_type), pointer :: reaction
  PetscReal :: porosity
  PetscInt :: ngcells
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: res_react(reaction%naqcomp)
  PetscReal :: jac_react(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: i,j
  class(material_auxvar_type), allocatable :: material_auxvar
  
  ngcells = sec_transport_vars%ncells
  vol = sec_transport_vars%vol     
                   
  res_react = 0.d0
  jac_react = 0.d0 ! These are not used anyway
  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  do i = 1, ngcells
    material_auxvar%porosity = porosity
    material_auxvar%volume = vol(i)
    call RReaction(res_react,jac_react,PETSC_FALSE, &
                   sec_transport_vars%sec_rt_auxvar(i), &
                   global_auxvars,material_auxvar,reaction,option)
  enddo
  call MaterialAuxVarStrip(material_auxvar)
  deallocate(material_auxvar)
  
  if (reaction%mineral%nkinmnrl > 0) then
    do i = 1, ngcells
      do j = 1, reaction%mineral%nkinmnrl
        sec_transport_vars%sec_rt_auxvar(i)%mnrl_volfrac(j) = &
          sec_transport_vars%sec_rt_auxvar(i)%mnrl_volfrac(j) + &
          sec_transport_vars%sec_rt_auxvar(i)%mnrl_rate(j)* &
          reaction%mineral%kinmnrl_molar_vol(j)* &
          option%tran_dt
          if (sec_transport_vars%sec_rt_auxvar(i)%mnrl_volfrac(j) < 0.d0) &
            sec_transport_vars%sec_rt_auxvar(i)%mnrl_volfrac(j) = 0.d0
      enddo
    enddo
  endif

  
end subroutine SecondaryRTUpdateKineticState

! ************************************************************************** !


subroutine SecondaryRTCheckResidual_np(sec_transport_vars,auxvar, &
                                    global_auxvar, &
                                    reaction,rt_parameter, &
                                    diffusion_coefficient, &
                                    porosity, tortuosity, option, &
                                    inf_norm_sec)
  ! 
  ! The residual of the secondary domain are checked
  ! to ensure convergence
  ! 
  ! Author: Albert Nardi, Amphos 21
  ! Date: 8/31/2021
  ! 
                                    
  use Option_module 
  use Global_Aux_module
  use Block_Solve_module
  use Block_Tridiag_module
  use Utility_module
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_module

  implicit none
  
  type(reactive_transport_param_type) :: rt_parameter
  type(sec_transport_type) :: sec_transport_vars
  type(reactive_transport_auxvar_type) :: auxvar
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type), pointer :: reaction
  type(option_type) :: option
  
  PetscReal :: res(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: curConc_up, curConc_dn 
  PetscReal :: sum_transp_up, sum_transp_dn
  PetscReal :: curCharge 
  PetscReal :: cur_diff 
  PetscReal :: curStoich
  PetscReal :: sum_denom_up, sum_denom_dn
  PetscReal :: potent_up(reaction%naqcomp), potent_dn(reaction%naqcomp)
  PetscReal :: factor_up_em(reaction%naqcomp), factor_dn_em(reaction%naqcomp)
  PetscReal :: conc_upd(reaction%naqcomp, sec_transport_vars%ncells) 
  PetscReal :: total_upd(reaction%naqcomp, sec_transport_vars%ncells)
  PetscReal :: total_prev(reaction%naqcomp, sec_transport_vars%ncells)
  PetscReal :: mc_pri_molal(reaction%naqcomp, sec_transport_vars%ncells)
  PetscReal :: mc_sec_molal(reaction%neqcplx, sec_transport_vars%ncells)
  PetscReal :: total_primary_node(reaction%naqcomp)
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: vol(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)
  PetscReal :: res_react(reaction%naqcomp)
  PetscReal :: jac_react(reaction%naqcomp,reaction%naqcomp)
  PetscInt :: i, j, k, n
  PetscInt :: ngcells, ncomp
  PetscInt :: ni
  PetscInt :: nicomp, icomp
  PetscReal :: area_fm
  PetscReal :: diffusion_coefficient
  PetscReal :: porosity
  PetscReal :: arrhenius_factor
  PetscReal :: pordt, pordiff
  PetscReal :: pordiff_prim, pordiff_sec
  PetscReal :: tortuosity
  PetscReal :: portort
  PetscInt :: jcomp, lcomp, kcomp, icplx, ncompeq
  PetscReal :: inf_norm_sec
  class(material_auxvar_type), allocatable :: material_auxvar

  PetscReal :: total_sorb_upd(reaction%naqcomp,sec_transport_vars%ncells) 
  PetscReal :: total_sorb_prev(reaction%naqcomp,sec_transport_vars%ncells)
  
  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  vol = sec_transport_vars%vol          
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus
  area_fm = sec_transport_vars%interfacial_area
  ncomp = reaction%naqcomp

  do j = 1, ncomp
    do i = 1, ngcells
      total_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total(j,1)
      if (reaction%neqsorb > 0) then
        total_sorb_prev(j,i) = sec_transport_vars%sec_rt_auxvar(i)%total_sorb_eq(j)
      endif
    enddo
  enddo
  conc_upd = sec_transport_vars%updated_conc
    
  ! Note that sec_transport_vars%sec_rt_auxvar(i)%pri_molal(j) units are in mol/kg
  ! Need to convert to mol/L since the units of total. in the Thomas 
  ! algorithm are in mol/L
  res = 0.d0
  
  total_primary_node = auxvar%total(:,1)                         ! in mol/L 
  pordt = porosity/option%tran_dt
  !pordiff = porosity*diffusion_coefficient
  !pordiff = porosity*diffusion_coefficient*tortuosity

  call RTAuxVarInit(rt_auxvar,reaction,option)
  do i = 1, ngcells
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i),option)
    rt_auxvar%pri_molal = conc_upd(:,i)
    ! sec_transport_vars%sec_rt_auxvar(i)%pri_molal =  conc_upd(:,i)
    mc_pri_molal(:, i) =  conc_upd(:,i)
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    if (reaction%neqsorb > 0) then
      ! call SecondaryRTotalSorb(rt_auxvar,global_auxvar,material_auxvar,reaction,option)
    endif
    total_upd(:,i) = rt_auxvar%total(:,1)
    ! sec_transport_vars%sec_rt_auxvar(i)%sec_molal = rt_auxvar%sec_molal
    mc_sec_molal(:, i) =  rt_auxvar%sec_molal
    if (reaction%neqsorb > 0) then 
      total_sorb_upd(:,i) = rt_auxvar%total_sorb_eq(:)
    endif
  enddo
                                    
!================ Calculate the secondary residual =============================        
  
  do i = 1, ngcells
  
    ! D_j \nabla C_j
    do icomp = 1, ncomp

      pordiff_prim = porosity*rt_parameter%pri_spec_diff_coef(icomp)*&
    	             tortuosity*global_auxvar%den_kg(1)/1000.d0
      n = icomp + (i-1)*ncomp
      
      ! Accumulation        
      res(n) = pordt*(total_upd(icomp,i) - total_prev(icomp,i))*vol(i)    ! in mol/L*m3/s
      if (reaction%neqsorb > 0) then 
        res(n) = res(n) + vol(i)/option%tran_dt*(total_sorb_upd(icomp,i) - total_sorb_prev(icomp,i))
      endif      

      ! Flux terms
      if (i.gt.1.and.i.lt.ngcells) then
        res(n) = res(n) - pordiff_prim*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                        (mc_pri_molal(icomp, i+1) - &
                        mc_pri_molal(icomp, i))
        res(n) = res(n) + pordiff_prim*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                        (mc_pri_molal(icomp, i) - &
                        mc_pri_molal(icomp, i-1)) 

      ! Apply boundary conditions
      ! Inner boundary
      else if (i.eq.1) then
        res(n) = res(n) - pordiff_prim*area(1)/(dm_minus(2) + dm_plus(1))* &
                        (mc_pri_molal(icomp, i+1) - &
                        mc_pri_molal(icomp, i))
           
      ! Outer boundary -- closest to primary node
      else !if (i.eq.ngcells) then
        res(n) = res(n) - &
                 pordiff_prim*area(ngcells)/dm_plus(ngcells)* &
                 (auxvar%pri_molal(icomp) - &
                 mc_pri_molal(icomp, i))
        res(n) = res(n) + &
                 pordiff_prim*area(ngcells-1)/(dm_minus(ngcells) &
                 + dm_plus(ngcells-1))* &
                 (mc_pri_molal(icomp, i) - &
                 mc_pri_molal(icomp, i-1))
      endif  

    enddo
    
    
    ! summatory D_i \nabla C_i   
    do icplx = 1, reaction%neqcplx
    
      pordiff_sec = porosity*rt_parameter%sec_spec_diff_coef(icplx)*&
    	       tortuosity*global_auxvar%den_kg(1)/1000
      nicomp = reaction%eqcplxspecid(0,icplx)
      
      do ni = 1, nicomp
        icomp = reaction%eqcplxspecid(ni,icplx)
        
        n = icomp + (i-1)*ncomp

        ! Flux terms
        if (i.gt.1.and.i.lt.ngcells) then
          res(n) = res(n) - pordiff_sec*area(i)/(dm_minus(i+1) + dm_plus(i))* &
                   (mc_sec_molal(icplx, i+1) - &
                   mc_sec_molal(icplx, i))* &
                   reaction%eqcplxstoich(ni,icplx)
          res(n) = res(n) + pordiff_sec*area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                   (mc_sec_molal(icplx, i) - &
                   mc_sec_molal(icplx, i-1))* &
                   reaction%eqcplxstoich(ni,icplx)
        
        ! Apply boundary conditions
        ! Inner boundary
        else if (i.eq.1) then
          res(n) = res(n) - pordiff_sec*area(1)/(dm_minus(2) + dm_plus(1))* &
                   (mc_sec_molal(icplx, i+1) - &
                   mc_sec_molal(icplx, i))* &
                   reaction%eqcplxstoich(ni,icplx)
        
        ! Outer boundary -- closest to primary node
        else !if (i.eq.ngcells) then
          res(n) = res(n) - &
                   pordiff_sec*area(ngcells)/dm_plus(ngcells)* &
                   (auxvar%sec_molal(icplx) - &
                   mc_sec_molal(icplx, i))* & 
                   reaction%eqcplxstoich(ni,icplx)
          res(n) = res(n) + pordiff_sec* &
                   area(ngcells-1)/(dm_minus(ngcells) + dm_plus(ngcells-1))* &
                   (mc_sec_molal(icplx, i) - &
                   mc_sec_molal(icplx, i-1))* &
                   reaction%eqcplxstoich(ni,icplx)
        endif

      enddo
    enddo
    
    
    ! Electromigration term
    call ComputeElectricPotentialTotalComponent_NP(i,reaction, &
                                            rt_parameter,  &
                                            rt_auxvar, &
                                            sec_transport_vars, &
                                            icomp,auxvar, &
                                            potent_dn, &
                                            potent_up, &
                                            mc_pri_molal, &
                                            mc_sec_molal)

    sum_denom_dn = 1d-40
    sum_denom_up = 1d-40
    
    do jcomp = 1, rt_parameter%naqcomp
      curCharge = reaction%primary_spec_Z(jcomp)
      sum_denom_up = sum_denom_up + potent_up(jcomp)*curCharge
      sum_denom_dn = sum_denom_dn + potent_dn(jcomp)*curCharge
    enddo

    do jcomp = 1, rt_parameter%naqcomp
      factor_up_em(jcomp) = potent_up(jcomp) / sum_denom_up
      factor_dn_em(jcomp) = potent_dn(jcomp) / sum_denom_dn
    enddo

    portort = porosity*tortuosity*global_auxvar%den_kg(1)/1000
    
    sum_transp_dn = 0d0
    sum_transp_up = 0d0
    
    ! Primary l index
    do lcomp = 1, rt_parameter%naqcomp
      curCharge = reaction%primary_spec_Z(lcomp)
      cur_diff = rt_parameter%pri_spec_diff_coef(lcomp)
      
      ! Flux terms
      if (i.gt.1.and.i.lt.ngcells) then
        curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                     (mc_pri_molal(lcomp, i+1) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                     (mc_pri_molal(lcomp, i) - &
                     mc_pri_molal(lcomp, i-1))
      
      ! Apply boundary conditions
      ! Inner boundary
      else if (i.eq.1) then
        curConc_up = area(1)/(dm_minus(2) + dm_plus(1))* &
                     (mc_pri_molal(lcomp, i+1) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = 0.d0
      
      ! Outer boundary    
      else !if (i.eq.ngcells) then
        curConc_up = area(ngcells)/dm_plus(ngcells)* &
                     (auxvar%pri_molal(lcomp) - &
                     mc_pri_molal(lcomp, i))
        curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                     + dm_plus(ngcells-1))* &
                     (mc_pri_molal(lcomp, i) - &
                     mc_pri_molal(lcomp, i-1))
      endif
      
      sum_transp_up = sum_transp_up - curCharge*cur_diff*(curConc_up)
      sum_transp_dn = sum_transp_dn + curCharge*cur_diff*(curConc_dn)

      ! Secondary term k
      do icplx = 1, reaction%neqcplx
        nicomp = reaction%eqcplxspecid(0,icplx)
        do ni = 1, nicomp
          kcomp = reaction%eqcplxspecid(ni,icplx)
          if (kcomp == lcomp) then  
            curCharge = reaction%primary_spec_Z(kcomp)
            cur_diff = rt_parameter%sec_spec_diff_coef(icplx)
            curStoich = reaction%eqcplxstoich(ni,icplx)
            
            ! Flux terms
            if (i.gt.1.and.i.lt.ngcells) then
              curConc_up = area(i)/(dm_minus(i+1) + dm_plus(i))* &
                          (mc_sec_molal(icplx, i+1) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = area(i-1)/(dm_minus(i) + dm_plus(i-1))* &
                          (mc_sec_molal(icplx, i) - &
                          mc_sec_molal(icplx, i-1))
            
            ! Apply boundary conditions
            ! Inner boundary
            else if (i.eq.1) then
              curConc_up = area(1)/(dm_minus(2) + dm_plus(1))* &
                          (mc_sec_molal(icplx, i+1) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = 0.d0
          
            ! Outer boundary   
            else !if (i.eq.ngcells) then 
              curConc_up = area(ngcells)/dm_plus(ngcells)* &
                          (auxvar%sec_molal(icplx) - &
                          mc_sec_molal(icplx, i))
              curConc_dn = area(ngcells-1)/(dm_minus(ngcells) &
                          + dm_plus(ngcells-1))* &
                          (mc_sec_molal(icplx, i) - &
                          mc_sec_molal(icplx, i-1))
            endif
            
            sum_transp_up = sum_transp_up - curCharge*cur_diff*curStoich*(curConc_up)
            sum_transp_dn = sum_transp_dn + curCharge*cur_diff*curStoich*(curConc_dn)
          endif
        enddo
      enddo      
  
    enddo

    do icomp = 1, rt_parameter%naqcomp
      n = icomp + (i-1)*ncomp
      res(n) = res(n) - portort*factor_dn_em(icomp)*sum_transp_dn
      res(n) = res(n) - portort*factor_up_em(icomp)*sum_transp_up
    enddo
   
  enddo
                        
  res = res*1.d3 ! Convert mol/L*m3/s to mol/s                                                                                 
                                    
!====================== Add reaction contributions =============================        
  
  ! Reaction 
  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  do i = 1, ngcells
    res_react = 0.d0
    jac_react = 0.d0
    call RTAuxVarCopy(rt_auxvar,sec_transport_vars%sec_rt_auxvar(i), &
                      option)
    rt_auxvar%pri_molal = conc_upd(:,i) ! in mol/kg
    call RTotalAqueous(rt_auxvar,global_auxvar,reaction,option)
    material_auxvar%porosity = porosity
    material_auxvar%volume = vol(i)
    call RReaction(res_react,jac_react,PETSC_FALSE, &
                   rt_auxvar,global_auxvar,material_auxvar,reaction,option)                     
    do j = 1, ncomp
      res(j+(i-1)*ncomp) = res(j+(i-1)*ncomp) + res_react(j) 
    enddo
  enddo  
  call MaterialAuxVarStrip(material_auxvar)         
  deallocate(material_auxvar)
  
 ! Need to decide how to scale the residual with volumes
  do i = 1, ngcells
    do j = 1, ncomp
      res(j+(i-1)*ncomp) = res(j+(i-1)*ncomp)/vol(i)
    enddo
  enddo
    
  inf_norm_sec = maxval(abs(res))  
  call RTAuxVarStrip(rt_auxvar)  
                                                                      
end subroutine SecondaryRTCheckResidual_np   

! ************************************************************************** !



subroutine SecondaryRTAuxVarComputeMulti(sec_transport_vars,reaction, &
                                         option)
  ! 
  ! Updates the secondary continuum
  ! concentrations at end of each time step for multicomponent system
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 2/1/13
  ! 
                               
                            
  use Option_module 
  use Reaction_Aux_module
  use Reaction_module
  use Reactive_Transport_Aux_module
  use Block_Solve_module
  use Block_Tridiag_module
  use Utility_module
  

  implicit none
  
  type(sec_transport_type) :: sec_transport_vars
  class(reaction_rt_type), pointer :: reaction
  type(option_type) :: option
  PetscReal :: coeff_left(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells)
  PetscReal :: coeff_diag(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells)
  PetscReal :: coeff_right(reaction%naqcomp,reaction%naqcomp, &
                 sec_transport_vars%ncells)
  PetscReal :: rhs(sec_transport_vars%ncells*reaction%naqcomp)
  PetscReal :: conc_upd(reaction%naqcomp,sec_transport_vars%ncells) 
  PetscInt :: i, j, n
  PetscInt :: ngcells, ncomp
  PetscInt :: pivot(reaction%naqcomp,sec_transport_vars%ncells)
  PetscInt :: indx(reaction%naqcomp)
  PetscInt :: d
    
  ngcells = sec_transport_vars%ncells
  ncomp = reaction%naqcomp
  ! Note that sec_transport_vars%sec_conc units are in mol/kg
  ! Need to convert to mol/L since the units of conc. in the Thomas 
  ! algorithm are in mol/L 
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0

  conc_upd = sec_transport_vars%updated_conc
           
  ! Use the stored coefficient matrices from LU decomposition of the
  ! block triagonal sytem
  coeff_left = sec_transport_vars%cxm
  coeff_right = sec_transport_vars%cxp
  coeff_diag = sec_transport_vars%cdl
  rhs = sec_transport_vars%r
    
  select case (option%secondary_continuum_solver)
    case(1) 
      call bl3dsolb(ngcells,ncomp,coeff_right,coeff_diag,coeff_left,pivot, &
                    ONE_INTEGER,rhs)
    case(2)
      call solbtb(ncomp,ngcells,ncomp,coeff_diag,coeff_right,coeff_left, &
                  pivot,rhs)
    case(3)
      do i = ngcells-1, 1, -1
        rhs(i) = (rhs(i) - coeff_right(ncomp,ncomp,i)*rhs(i+1))/ &
                             coeff_diag(ncomp,ncomp,i)
      enddo
    case default
      option%io_buffer = 'SECONDARY_CONTINUUM_SOLVER can be only ' // &
                         'HINDMARSH or KEARST. For single component'// &
                         'chemistry THOMAS can be used.'
      call PrintErrMsg(option)
  end select  
  
  do j = 1, ncomp
    do i = 1, ngcells
      n = j + (i - 1)*ncomp
      if (reaction%use_log_formulation) then 
        ! convert log concentration to concentration
        rhs(n) = dsign(1.d0,rhs(n))*min(dabs(rhs(n)),reaction%max_dlnC) 
        conc_upd(j,i) = exp(rhs(n))*conc_upd(j,i)
      else
        conc_upd(j,i) = rhs(n) + conc_upd(j,i)
      endif
      if (conc_upd(j,i) < 0.d0) conc_upd(j,i) = 1.d-8
    enddo
  enddo
  
  sec_transport_vars%updated_conc = conc_upd
    
end subroutine SecondaryRTAuxVarComputeMulti

! ************************************************************************** !

subroutine THCSecHeatAuxVarCompute(sec_heat_vars,global_auxvar, &
                                   therm_conductivity,dencpr, &
                                   option)
  ! 
  ! Computes secondary auxillary variables for each
  ! grid cell for heat transfer only
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/5/12
  ! 

  use Option_module 
  use Global_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: sec_temp(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
  temp_primary_node = global_auxvar%temp
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  sec_temp = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
  ! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0

                        
  rhs = sec_heat_vars%sec_temp  ! secondary continuum values from previous time step
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! Calculate temperature in the secondary continuum
  sec_temp(ngcells) = rhs(ngcells)/coeff_diag(ngcells)
  do i = ngcells-1, 1, -1
    sec_temp(i) = (rhs(i) - coeff_right(i)*sec_temp(i+1))/coeff_diag(i)
  enddo
  
  sec_heat_vars%sec_temp = sec_temp
            
end subroutine THCSecHeatAuxVarCompute

! ************************************************************************** !

subroutine THSecHeatAuxVarCompute(sec_heat_vars,global_auxvar, &
                                   therm_conductivity,dencpr, &
                                   option)
  ! 
  ! Computes secondary auxillary variables for each
  ! grid cell for heat transfer only
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/5/12
  ! 

  use Option_module 
  use Global_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: sec_temp(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
  temp_primary_node = global_auxvar%temp
  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  sec_temp = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
  ! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0

                        
  rhs = sec_heat_vars%sec_temp  ! secondary continuum values from previous time step
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! Calculate temperature in the secondary continuum
  sec_temp(ngcells) = rhs(ngcells)/coeff_diag(ngcells)
  do i = ngcells-1, 1, -1
    sec_temp(i) = (rhs(i) - coeff_right(i)*sec_temp(i+1))/coeff_diag(i)
  enddo
  
  sec_heat_vars%sec_temp = sec_temp
            
end subroutine THSecHeatAuxVarCompute

! ************************************************************************** !

subroutine MphaseSecHeatAuxVarCompute(sec_heat_vars,auxvar,global_auxvar, &
                                   therm_conductivity,dencpr, &
                                   option)
  ! 
  ! Computes secondary auxillary variables in each
  ! grid cell for heat transfer only
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/28/12
  ! 

  use Option_module 
  use Global_Aux_module
  use Mphase_Aux_module
  
  implicit none
  
  type(sec_heat_type) :: sec_heat_vars
  type(mphase_auxvar_elem_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: sec_temp(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  
  
  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
  temp_primary_node = auxvar%temp

  
  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0
  sec_temp = 0.d0
  
  alpha = option%flow_dt*therm_conductivity/dencpr

  
  ! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo
  
  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))
  
  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0
                        
  rhs = sec_heat_vars%sec_temp  ! secondary continuum values from previous time step
  rhs(ngcells) = rhs(ngcells) + & 
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node
                
  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! Calculate temperature in the secondary continuum
  sec_temp(ngcells) = rhs(ngcells)/coeff_diag(ngcells)
  do i = ngcells-1, 1, -1
    sec_temp(i) = (rhs(i) - coeff_right(i)*sec_temp(i+1))/coeff_diag(i)
  enddo

! print *,'temp_dcdm= ',(sec_temp(i),i=1,ngcells)
  
  sec_heat_vars%sec_temp = sec_temp


end subroutine MphaseSecHeatAuxVarCompute

! ************************************************************************** !

subroutine SecondaryRTGetVariable(realization, vec, ivar, isubvar, mc_layer)
 
  ! Extracts a secondary continuum variable for a layer of secondary
  ! continuum cells. Similar to RealizationGetVariable, but now the "layer"
  ! of the cells needs to be specified.

#include "petsc/finclude/petscvec.h"
  use Variables_module
  use Grid_module
  use Patch_module
  use Realization_Subsurface_class
  use petscvec

  implicit none
  
  class(realization_subsurface_type) :: realization
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: mc_layer

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  PetscReal, pointer :: vec_p(:)
  PetscInt :: local_id
  PetscErrorCode :: ierr

  patch => realization%patch
  grid => patch%grid

  call VecGetArrayF90(vec,vec_p,ierr);CHKERRQ(ierr)
  
  select case(ivar)
    case(SECONDARY_CONTINUUM_UPDATED_CONC)
      do local_id=1,grid%nlmax
        vec_p(local_id) = &
          patch%aux%SC_RT%sec_transport_vars(local_id)% &
          updated_conc(isubvar,mc_layer)
      enddo
    case(MINERAL_VOLUME_FRACTION)
      do local_id=1, grid%nlmax
        vec_p(local_id) = &
          patch%aux%SC_RT%sec_transport_vars(local_id)% &
          sec_rt_auxvar(mc_layer)%mnrl_volfrac(isubvar)
      enddo
    case(REACTION_AUXILIARY)
      do local_id=1, grid%nlmax
        vec_p(local_id) = &
          patch%aux%SC_RT%sec_transport_vars(local_id)% &
          sec_rt_auxvar(mc_layer)%auxiliary_data(isubvar)
      enddo
    case(PRIMARY_ACTIVITY_COEF)
      do local_id=1, grid%nlmax
        vec_p(local_id) = &
          patch%aux%SC_RT%sec_transport_vars(local_id)% &
          sec_rt_auxvar(mc_layer)%pri_act_coef(isubvar)
      enddo
    case(SECONDARY_ACTIVITY_COEF)
      do local_id=1, grid%nlmax
        vec_p(local_id) = &
          patch%aux%SC_RT%sec_transport_vars(local_id)% &
          sec_rt_auxvar(mc_layer)%sec_act_coef(isubvar)
      enddo
  end select
  
  call VecRestoreArrayF90(vec,vec_p,ierr);CHKERRQ(ierr)

end subroutine SecondaryRTGetVariable

! ************************************************************************** !

subroutine SecondaryRTSetVariable(realization, vec, vec_format, ivar, isubvar, mc_layer)
  
  ! Sets a secondary continuum variable to a layer of secondary
  ! continuum cells. Similar to RealizationSetVariable, but now the "layer"
  ! of the cells needs to be specified.

#include "petsc/finclude/petscvec.h"
  use Variables_module
  use Grid_module
  use Patch_module
  use Realization_Subsurface_class
  use Option_module
  use petscvec

  implicit none
  
  class(realization_subsurface_type) :: realization
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: mc_layer

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  PetscReal, pointer :: vec_p(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscErrorCode :: ierr
  patch => realization%patch
  grid => patch%grid

  if (vec_format == NATURAL .or. vec_format == LOCAL) then
    call PrintErrMsg(realization%option,&
                     'NATURAL and LOCAL vector formats not supported by &
SecondaryRTSetVariable')
  endif

  call VecGetArrayF90(vec,vec_p,ierr);CHKERRQ(ierr)
  
  select case(ivar)
    case(SECONDARY_CONTINUUM_UPDATED_CONC)
      do local_id=1, grid%nlmax
        patch%aux%SC_RT%sec_transport_vars(local_id)% &
          updated_conc(isubvar,mc_layer) = vec_p(local_id)
      enddo
    case(MINERAL_VOLUME_FRACTION)
      do local_id=1, grid%nlmax
        patch%aux%SC_RT%sec_transport_vars(local_id)% &
          sec_rt_auxvar(mc_layer)%mnrl_volfrac(isubvar) = vec_p(local_id)
      enddo
    case(REACTION_AUXILIARY)
      do local_id=1, grid%nlmax
        patch%aux%SC_RT%sec_transport_vars(local_id)% &
          sec_rt_auxvar(mc_layer)%auxiliary_data(isubvar) = vec_p(local_id)
      enddo
    case(PRIMARY_ACTIVITY_COEF)
      do local_id=1, grid%nlmax
        patch%aux%SC_RT%sec_transport_vars(local_id)% &
          sec_rt_auxvar(mc_layer)%pri_act_coef(isubvar) = vec_p(local_id)
      enddo
    case(SECONDARY_ACTIVITY_COEF)
      do local_id=1, grid%nlmax
        patch%aux%SC_RT%sec_transport_vars(local_id)% &
          sec_rt_auxvar(mc_layer)%sec_act_coef(isubvar) = vec_p(local_id)
      enddo
  end select
  
  call VecRestoreArrayF90(vec,vec_p,ierr);CHKERRQ(ierr)

end subroutine SecondaryRTSetVariable

! ************************************************************************** !

subroutine ComputeElectricPotentialTotalComponent_NP(ic,reaction, &
                                                    rt_parameter, &
                                                    rt_auxvar, &
                                                    sec_transport_vars, &
                                                    icomp,auxvar, &
                                                    potent_dn, &
                                                    potent_up, &
                                                    mc_pri_molal, &
                                                    mc_sec_molal) 
!
! Computes the U electric potential vector term
!
! Author: Amphos21 - Barcelona Science
! Date: 29/06/20
!

  use Reaction_Aux_module 
  use Reactive_Transport_Aux_module 

  class(reaction_rt_type), pointer :: reaction 
  type(reactive_transport_param_type) :: rt_parameter
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(reactive_transport_auxvar_type) :: auxvar 
  type(sec_transport_type) :: sec_transport_vars
  PetscReal :: mc_pri_molal(reaction%naqcomp,sec_transport_vars%ncells)
  PetscReal :: mc_sec_molal(reaction%neqcplx,sec_transport_vars%ncells)
  
  PetscInt :: icomp, jcomp
  PetscReal :: potent_up(reaction%naqcomp), potent_dn(reaction%naqcomp)
  PetscReal :: curCoef, curCharge, curStoich
  PetscReal :: curConc, curConc_up, curConc_dn
  PetscInt :: icplx, i, j, ncomp
  PetscInt :: ic 
  PetscReal :: a,b,c
  PetscInt :: ngcells
  PetscReal :: area(sec_transport_vars%ncells)
  PetscReal :: dm_plus(sec_transport_vars%ncells)
  PetscReal :: dm_minus(sec_transport_vars%ncells)
  ngcells = sec_transport_vars%ncells
  area = sec_transport_vars%area
  dm_plus = sec_transport_vars%dm_plus
  dm_minus = sec_transport_vars%dm_minus

  ! Primary species contribution
  do icomp = 1, reaction%naqcomp
    curCoef = rt_parameter%pri_spec_diff_coef(icomp)
    curCharge = reaction%primary_spec_Z(icomp)
    
    ! Flux terms
    if (ic.gt.1.and.ic.lt.ngcells) then
      a = mc_pri_molal(icomp, ic-1)
      b = mc_pri_molal(icomp, ic)
      c = mc_pri_molal(icomp, ic+1)
      
      ! use of differential logarithmic avg
      if (log(b).eq.log(c)) then
        curConc_up = b
      else
        curConc_up = (c-b)/ (log(c) - log(b))
      endif
      
      if (log(a).eq.log(b)) then
        curConc_dn = b
      else
        curConc_dn = (b-a)/ (log(b) - log(a))
      endif
    
    ! Apply boundary conditions
    ! Inner boundary
    else if (ic.eq.1) then
      b = mc_pri_molal(icomp, ic)
      c = mc_pri_molal(icomp, ic+1)
      
      ! use of differential logarithmic avg
      if (log(b).eq.log(c)) then
        curConc_up = b
      else
        curConc_up = (c-b)/ (log(c) - log(b))
      endif
      
      curConc_dn = b
    
    ! Outer boundary -- closest to primary node
    else !if (i.eq.ngcells) then
      a = mc_pri_molal(icomp, ic-1)
      b = mc_pri_molal(icomp, ic)
      c = auxvar%pri_molal(icomp)
      
      ! use of differential logarithmic avg
      if (log(a).eq.log(b)) then
        curConc_dn = b
      else
        curConc_dn = (b-a)/ (log(b) - log(a))
      endif

      if (log(b).eq.log(c)) then
        curConc_up = b
      else
        curConc_up = (c-b)/ (log(c) - log(b))
      endif

    endif
    
    potent_up(icomp) = curCoef*curCharge*curConc_up
    potent_dn(icomp) = curCoef*curCharge*curConc_dn
  enddo


  ! Secondary species contribution
  do icplx = 1, reaction%neqcplx
    curCharge = reaction%eqcplx_Z(icplx)
    ncomp = reaction%eqcplxspecid(0,icplx)
    do i = 1, ncomp
      jcomp = reaction%eqcplxspecid(i,icplx)
      curCoef = rt_parameter%sec_spec_diff_coef(icplx)
      curStoich = reaction%eqcplxstoich(i,icplx)

      ! Flux terms
      if (ic.gt.1.and.ic.lt.ngcells) then
        a = mc_sec_molal(icplx, ic-1)
        b = mc_sec_molal(icplx, ic)
        c = mc_sec_molal(icplx, ic+1)
      
        ! use of differential logarithmic avg
        if (log(b).eq.log(c)) then
          curConc_up = b
        else
          curConc_up = (c-b)/ (log(c) - log(b))
        endif
      
        if (log(a).eq.log(b)) then
          curConc_dn = b
        else
          curConc_dn = (b-a)/ (log(b) - log(a))
        endif
      
      ! Apply boundary conditions
      ! Inner boundary
      else if (ic.eq.1) then
        b = mc_sec_molal(icplx, ic)
        c = mc_sec_molal(icplx, ic+1)
      
        ! use of differential logarithmic avg
        if (log(b).eq.log(c)) then
          curConc_up = b
        else
          curConc_up = (c-b)/ (log(c) - log(b))
        endif

        curConc_dn = b

      ! Outer boundary -- closest to primary node
      else !if (i.eq.ngcells) then
        a = mc_sec_molal(icplx, ic-1)
        b = mc_sec_molal(icplx, ic)
        c = auxvar%sec_molal(icplx)
      
        ! use of differential logarithmic avg
        if (log(b).eq.log(c)) then
          curConc_up = b
        else
          curConc_up = (c-b)/ (log(c) - log(b))
        endif
      
        if (log(a).eq.log(b)) then
          curConc_dn = b
        else
          curConc_dn = (b-a)/ (log(b) - log(a))
        endif

      endif 

      potent_up(jcomp) = potent_up(jcomp) + curStoich * curCoef * curCharge * curConc_up
      potent_dn(jcomp) = potent_dn(jcomp) + curStoich * curCoef * curCharge * curConc_dn
    enddo
  enddo

end subroutine ComputeElectricPotentialTotalComponent_NP

end module Secondary_Continuum_NP_module
            