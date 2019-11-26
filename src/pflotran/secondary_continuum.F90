! added by S. Karra 07/11/12

module Secondary_Continuum_module
  
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
            THCSecHeatAuxVarCompute, &
            THSecHeatAuxVarCompute

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


end module Secondary_Continuum_module
            
