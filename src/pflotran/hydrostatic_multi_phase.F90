module HydrostaticMultiPhase_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Hydrostatic_Common_module
  use Lookup_Table_module

  implicit none

  private

  !module vriables
  PetscBool :: wat_gas_equil
  PetscBool :: oil_wat_zone
  PetscBool :: oil_gas_zone
  PetscBool :: wat_gas_zone
  PetscBool :: wat_present
  PetscBool :: oil_present
  PetscBool :: gas_present
  PetscBool :: datum_in_wat
  PetscBool :: datum_in_oil
  PetscBool :: datum_in_gas

  PetscReal :: datum_z
  PetscReal :: press_at_datum
  PetscReal :: owc_z
  PetscReal :: pcow_owc
  PetscReal :: ogc_z
  PetscReal :: pcog_ogc
  PetscReal :: wgc_z
  PetscReal :: pcwg_wgc
  PetscReal :: z_min
  PetscReal :: z_max
  
  PetscReal :: res_temp_const
  PetscBool :: res_temp_is_const

  PetscReal :: pb_const
  PetscBool :: pb_is_const

  PetscReal :: g_z

  PetscReal :: pw_owc
  PetscReal :: po_owc
  PetscReal :: po_ogc
  PetscReal :: pg_ogc
  PetscReal :: pw_wgc
  PetscReal :: pg_wgc
  
  class(lookup_table_general_type), pointer :: rtempvz_table => null()
  class(lookup_table_general_type), pointer :: pbvz_table => null()
  PetscBool :: temp_grad_given
  PetscReal :: temp_grad(3)
  PetscReal :: temp_at_datum

  !working arrays
  PetscReal, allocatable :: temp_vec(:)
  PetscReal, allocatable :: xm_nacl_vec(:)
  PetscReal, allocatable :: pb_vec(:)
  PetscReal, allocatable :: rv_vec(:)
  PetscReal, allocatable :: wat_press_vec(:)
  PetscReal, allocatable :: oil_press_vec(:)
  PetscReal, allocatable :: gas_press_vec(:)  
  PetscReal, allocatable :: wat_den_kg_vec(:)
  PetscReal, allocatable :: oil_den_kg_vec(:)
  PetscReal, allocatable :: gas_den_kg_vec(:)

  !parameters
  PetscReal, parameter :: z_eps = 1.0d-6

  type :: one_dim_grid_type
    PetscReal :: delta_z
    PetscReal :: min_z
    PetscReal :: max_z
    PetscReal, pointer :: z(:)
    PetscInt :: idatum
  contains 
    procedure :: ZLookup
  end type one_dim_grid_type


  public :: TOIHydrostaticUpdateCoupler, &
            HydrostaticMPUpdateCoupler
  !          HydrostaticTest 
 
contains

  ! ************************************************************************** !
subroutine HydrostaticMPUpdateCoupler(coupler,option,grid, &
                     characteristic_curves_array,sat_func_id,imat)
  ! 
  ! Computes the hydrostatic initial/boundary condition for multiphase modes
  ! given: 
  ! (1) oil water contact (owc) elevation
  ! (2) Pcow at owc
  ! (3) oil gas contact (ogc) elevation
  ! (4) Pcog at ogc
  ! (5) Reference pressure at a given datum
  !
  ! Author: Paolo Orsini
  ! Date: 01/05/19

  use Option_module
  use Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Characteristic_Curves_module
  use Utility_module

  implicit none

  type(coupler_type), intent(inout) :: coupler
  type(option_type) :: option
  type(grid_type), intent(in) :: grid
  type(characteristic_curves_ptr_type), intent(in) :: characteristic_curves_array(:)
  PetscInt, pointer, intent(in) :: sat_func_id(:)
  PetscInt, pointer, intent(in) :: imat(:)

  type(flow_condition_type), pointer :: condition
  class(characteristic_curves_type), pointer :: cc_ptr
  class(one_dim_grid_type), pointer :: one_d_grid
  PetscInt :: i_z
  PetscInt :: grid_size, i_press_start, i_n
  PetscReal :: dist_dn, dist_up
  PetscInt :: iconn, local_id, ghosted_id
  PetscReal :: dz_conn, z_cell
  PetscInt :: tbl_size

  PetscInt :: i_owc
  PetscInt :: i_ogc
  PetscInt :: i_wgc
  PetscReal :: dowc_dn, dowc_up
  PetscReal :: dogc_dn, dogc_up
  PetscReal :: dwgc_dn, dwgc_up
  
  PetscReal :: pw_cell
  PetscReal :: po_cell
  PetscReal :: pg_cell
  PetscReal :: pow_cell
  PetscReal :: pog_cell
  PetscReal :: pwg_cell
  PetscReal :: sw_cell
  PetscReal :: so_cell
  PetscReal :: sg_cell
  
  PetscReal :: pcow_min
  PetscReal :: pcow_max
  PetscReal :: pcow
  PetscReal :: pcog_min
  PetscReal :: pcog_max
  PetscReal :: pcog
  PetscReal :: pcwg_min
  PetscReal :: pcwg_max  
  PetscReal :: sw_max
  PetscReal :: sw_min
  PetscReal :: sg_max
  PetscReal :: sg_min
  PetscReal :: so_min
  PetscReal :: dsat_dp
  PetscReal :: dp_dsat

  PetscReal :: temp_owc
  PetscReal :: xm_nacl_owc
  PetscReal :: pb_owc
  PetscReal :: rv_owc
  
  PetscReal :: temp_ogc
  PetscReal :: xm_nacl_ogc
  PetscReal :: pb_ogc
  PetscReal :: rv_ogc
  
  PetscReal :: temp_wgc
  PetscReal :: xm_nacl_wgc
  PetscReal :: pb_wgc
  PetscReal :: rv_wgc
  
  PetscReal :: press_start
    
  PetscReal :: gravity_magnitude
  PetscReal :: datum_tmp_3v(3)
  
  !PetscReal :: max_z, min_z


  condition => coupler%flow_condition

  call HydrostaticPMLoader(condition,option)
  ! call HydrostaticPMLoader(condition,option,datum_z,press_at_datum, &
  !                          owc_z,pcow_owc,ogc_z,pcog_ogc,wgc_z,pcwg_wgc)  
  
  !determine which are the phase and transition zone present
  z_max = max(grid%z_max_global,datum_z)+1.d0 ! add 1m buffer
  z_min = min(grid%z_min_global,datum_z)-1.d0

  if (temp_grad_given) call RTempTableFromGrad(option)

  if ( associated(rtempvz_table) ) then
    if ( (rtempvz_table%axis1%values(1) - 1.0) > z_min ) then
      option%io_buffer = 'Hydrostatic equilibration: Temperature table - &
                          &range does not cover deepest layers'
      call printErrMsg(option)
    end if
    tbl_size = size(rtempvz_table%axis1%values(:))
    if ( (rtempvz_table%axis1%values(tbl_size) + 1.0) < z_max  ) then
      option%io_buffer = 'Hydrostatic equilibration: Temperature table - &
                          &range does not cover most shallow layers'
      call printErrMsg(option)
    end if
  endif
  ! determine what are the transition zones present
  if ( owc_z >= z_min .and. owc_z <= z_max ) then
    oil_wat_zone = PETSC_TRUE
  end if

  if ( ogc_z >= z_min .and. ogc_z <= z_max ) then
    oil_gas_zone = PETSC_TRUE
  end if

  !check that water is heavier than oil and oil heavier than gas
  if ( owc_z > ogc_z .and. oil_wat_zone .and. oil_gas_zone ) then
    option%io_buffer = 'Hydrostatic Equilibration input error: owc > ogc'
    call printErrMsg(option)
  end if

  if ( dabs(ogc_z - owc_z) < 0.1d0 .and. oil_wat_zone .and. oil_gas_zone ) then
    option%io_buffer = 'Hydrostatic Equilibration error: owc and ogc &
                        &less than 10 cm apart, consider to switch to wgc'
    call printErrMsg(option)
  end if

  ! wat_gas_zone determined in HydrostaticPMLoader
  if (wat_gas_equil) then
    oil_wat_zone = PETSC_FALSE
    oil_gas_zone = PETSC_FALSE
    if ( wgc_z >= z_min .and. wgc_z <= z_max) then
      wat_gas_zone = PETSC_TRUE
    end if
  end if
  
  !works out what are the phases present in the model
  if ( oil_wat_zone .or. wat_gas_zone ) then
    wat_present = PETSC_TRUE
  end if
    
  if ( oil_wat_zone .or. oil_gas_zone ) then
    oil_present = PETSC_TRUE
  end if

  if ( oil_gas_zone .or. wat_gas_zone ) then
    gas_present = PETSC_TRUE
  end if

  !detect single phase cases
  if ( (.not.oil_wat_zone) .and. (.not.oil_wat_zone) .and. &
     (  .not.wat_gas_equil) ) then
    if ( owc_z < z_min .and. ogc_z > z_max ) then
      oil_present = PETSC_TRUE
    else if ( owc_z < z_min .and. ogc_z < z_min ) then
      gas_present = PETSC_TRUE
    else if ( owc_z > z_max .and. ogc_z > z_max ) then
      wat_present = PETSC_TRUE
    end if
  end if

  if ( wat_gas_equil .and. (.not.wat_gas_zone) ) then
    if ( wgc_z > z_max ) then
      wat_present = PETSC_TRUE
    else if (wgc_z < z_min ) then
      gas_present = PETSC_TRUE
    end if
  end if

  !determine where is the datum
  datum_in_wat = PETSC_FALSE
  datum_in_oil = PETSC_FALSE
  datum_in_gas = PETSC_FALSE  
  
  if ( datum_z <= owc_z ) then
    datum_in_wat = PETSC_TRUE !datum is in the water region
  else if ( datum_z > owc_z .and. datum_z <= ogc_z ) then
    datum_in_oil = PETSC_TRUE
  else if ( datum_z > ogc_z ) then
    datum_in_gas = PETSC_TRUE
  end if

  !check gravity value
  gravity_magnitude = sqrt(DotProduct(option%gravity,option%gravity))
  g_z = option%gravity(Z_DIRECTION)
  if (dabs(gravity_magnitude-EARTH_GRAVITY) > 0.1d0) then
    option%io_buffer = 'Magnitude of gravity vector is not near 9.81.'
    call printErrMsg(option)
  endif

  !create 1d grid
  !replace datum_tmp_3v with datum_z when replacing TOIHydrostaticUpdateCoupler
  datum_tmp_3v(1:3) = 0.0d0
  datum_tmp_3v(3) = datum_z
  one_d_grid => CreateOneDimGrid(z_min,z_max,datum_tmp_3v)
  
  grid_size = size(one_d_grid%z(:))

  call AllocateGridArrays(grid_size)

  if (oil_wat_zone) then
    call GetZPointProps(owc_z,option,temp_owc,xm_nacl_owc,pb_owc,rv_owc)
    call one_d_grid%ZLookup(owc_z,i_owc,dowc_dn,dowc_up)
  end if
  if (oil_gas_zone) then
    call GetZPointProps(ogc_z,option,temp_ogc,xm_nacl_ogc,pb_ogc,rv_ogc)
    call one_d_grid%ZLookup(ogc_z,i_ogc,dogc_dn,dogc_up)
  end if
  if (wat_gas_zone) then
    call GetZPointProps(wgc_z,option,temp_wgc,xm_nacl_wgc,pb_wgc,rv_wgc)
    call one_d_grid%ZLookup(wgc_z,i_wgc,dwgc_dn,dwgc_up)
  end if  
  

  do i_z = 1,size(one_d_grid%z(:))
    call GetZPointProps(one_d_grid%z(i_z),option,temp_vec(i_z), &
                        xm_nacl_vec(i_z),pb_vec(i_z),rv_vec(i_z))
  end do
  
  ! compute hydrostatic pressure profiles for exisitng phases
  if (datum_in_wat) then
    call PhaseHydrostaticPressure(one_d_grid,option,LIQUID_PHASE, &
                press_at_datum,one_d_grid%idatum,temp_vec, &
                xm_nacl_vec,pb_vec,rv_vec, &
                wat_press_vec,wat_den_kg_vec)
    if ( oil_present ) then
      ! compute pw_owc and starts po computation from pw_owc
      pw_owc = wat_press_vec(i_owc) + wat_den_kg_vec(i_owc) * g_z * dowc_dn
      po_owc = pw_owc + pcow_owc
      i_press_start = i_owc + 1
      press_start = po_owc + PhaseDensity(option,OIL_PHASE,po_owc,temp_owc, &
                             xm_nacl_owc,pb_owc,rv_owc) * g_z * dowc_up
      !
      call PhaseHydrostaticPressure(one_d_grid,option,OIL_PHASE, &
                  press_start,i_press_start,temp_vec, &
                  xm_nacl_vec,pb_vec,rv_vec, &
                  oil_press_vec,oil_den_kg_vec)
    end if
    if ( gas_present .and. oil_present ) then
      !compute po_ogc and pg_ogc, then gas pressure profile from pg_ogc
      po_ogc = oil_press_vec(i_ogc) + oil_den_kg_vec(i_ogc) * g_z * dogc_dn
      pg_ogc = po_ogc + pcog_ogc
      i_press_start = i_ogc + 1
      press_start = pg_ogc + PhaseDensity(option,GAS_PHASE,pg_ogc,temp_ogc, &
                             xm_nacl_ogc,pb_ogc,rv_ogc) * g_z * dogc_up
    else if ( gas_present .and. (.not.oil_present) ) then
      !compute pw_wgc and pg_wgc, then gas pressure profile from pg_wgc
      pw_wgc = wat_press_vec(i_wgc) + wat_den_kg_vec(i_wgc) * g_z * dwgc_dn
      pg_wgc = pw_wgc + pcwg_wgc
      i_press_start = i_wgc + 1
      press_start = pg_wgc + PhaseDensity(option,GAS_PHASE,pg_wgc,temp_wgc, &
                             xm_nacl_wgc,pb_wgc,rv_wgc) * g_z * dwgc_up
    end if
    if (gas_present) then
      call PhaseHydrostaticPressure(one_d_grid,option,GAS_PHASE, &
                  press_start,i_press_start,temp_vec, &
                  xm_nacl_vec,pb_vec,rv_vec, &
                  gas_press_vec,gas_den_kg_vec)
    end if
  else if (datum_in_oil) then
    call PhaseHydrostaticPressure(one_d_grid,option,OIL_PHASE, &
                press_at_datum,one_d_grid%idatum,temp_vec, &
                xm_nacl_vec,pb_vec,rv_vec, &
                oil_press_vec,oil_den_kg_vec)
    if (wat_present) then
      po_owc = oil_press_vec(i_owc) + oil_den_kg_vec(i_owc) * g_z * dowc_dn
      pw_owc = po_owc - pcow_owc
      i_press_start = i_owc + 1
      press_start = pw_owc + PhaseDensity(option,LIQUID_PHASE,pw_owc,temp_owc, &
                             xm_nacl_owc,pb_owc,rv_owc) * g_z * dowc_up
      !
      call PhaseHydrostaticPressure(one_d_grid,option,LIQUID_PHASE, &
                  press_start,i_press_start,temp_vec, &
                  xm_nacl_vec,pb_vec,rv_vec, &
                  wat_press_vec,wat_den_kg_vec)
    end if  
    if (gas_present) then
      po_ogc = oil_press_vec(i_ogc) + oil_den_kg_vec(i_ogc) * g_z * dogc_dn
      pg_ogc = po_ogc + pcog_ogc
      i_press_start = i_ogc + 1
      press_start = pg_ogc + PhaseDensity(option,GAS_PHASE,pg_ogc,temp_ogc, &
                             xm_nacl_ogc,pb_ogc,rv_ogc) * g_z * dogc_up
      !
      call PhaseHydrostaticPressure(one_d_grid,option,GAS_PHASE, &
                  press_start,i_press_start,temp_vec, &
                  xm_nacl_vec,pb_vec,rv_vec, &
                  gas_press_vec,gas_den_kg_vec)
    end if
    !check what phase are above water
    !write z function lookup
  else if (datum_in_gas) then
    call PhaseHydrostaticPressure(one_d_grid,option,GAS_PHASE, &
                press_at_datum,one_d_grid%idatum,temp_vec, &
                xm_nacl_vec,pb_vec,rv_vec, &
                gas_press_vec,gas_den_kg_vec)
    if (oil_present) then
      !compute po_ogc then pressure profile
      pg_ogc = gas_press_vec(i_ogc) + gas_den_kg_vec(i_ogc) * g_z * dogc_dn
      po_ogc = pg_ogc - po_ogc
      i_press_start = i_ogc + 1
      press_start = po_ogc + PhaseDensity(option,OIL_PHASE,po_ogc,temp_ogc, &
                             xm_nacl_ogc,pb_ogc,rv_ogc) * g_z * dogc_up
      !
      call PhaseHydrostaticPressure(one_d_grid,option,OIL_PHASE, &
                  press_start,i_press_start,temp_vec, &
                  xm_nacl_vec,pb_vec,rv_vec, &
                  oil_press_vec,oil_den_kg_vec)
    end if
    if ( wat_present .and. oil_present ) then
      !compute po_owc and pw_owc, then water pressure profile from pw_owc
      po_owc = oil_press_vec(i_owc) + oil_den_kg_vec(i_owc) * g_z * dowc_dn
      pw_owc = po_owc - pcow_owc
      i_press_start = i_owc + 1
      press_start = pw_owc + PhaseDensity(option,LIQUID_PHASE,pw_owc,temp_owc, &
                             xm_nacl_owc,pb_owc,rv_owc) * g_z * dowc_up
    else if ( wat_present .and. (.not.oil_present) ) then
      !compute pg_wgc and pw_wgc, then water pressure profile
      pg_wgc = gas_press_vec(i_wgc) + gas_den_kg_vec(i_wgc) * g_z * dwgc_dn
      pw_wgc = pg_wgc - pcwg_wgc
      i_press_start = i_wgc + 1
      press_start = pw_wgc + PhaseDensity(option,LIQUID_PHASE,pw_wgc,temp_wgc, &
                             xm_nacl_wgc,pb_wgc,rv_wgc) * g_z * dwgc_up
    end if
    if (wat_present) then
      call PhaseHydrostaticPressure(one_d_grid,option,LIQUID_PHASE, &
                  press_start,i_press_start,temp_vec, &
                  xm_nacl_vec,pb_vec,rv_vec, &
                  wat_press_vec,wat_den_kg_vec)
    end if  
  end if
  
  ! compute saturation for exisitng transition zones
  do iconn=1,coupler%connection_set%num_connections 
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
 
    dz_conn = 0.0d0
    ! geh: note that this is a boundary connection, thus the entire distance is between
    ! the face and cell center
    if (associated(coupler%connection_set%dist)) then
      dz_conn = coupler%connection_set%dist(0,iconn) * &
                      coupler%connection_set%dist(3,iconn)
    endif

    ! note the negative (-) dz_conn is required due to the offset of the boundary face
    z_cell = grid%z(ghosted_id)-dz_conn
  
    !cell z lookup for
    call one_d_grid%ZLookup(z_cell,i_n,dist_dn,dist_up)

    !cell pressure computations
    if (wat_present) then
      pw_cell = wat_press_vec(i_n) + wat_den_kg_vec(i_n) * g_z * dist_dn
    end if   

    if (oil_present) then
      po_cell = oil_press_vec(i_n) + oil_den_kg_vec(i_n) * g_z * dist_dn
    end if

    if (gas_present) then
      pg_cell = gas_press_vec(i_n) + gas_den_kg_vec(i_n) * g_z * dist_dn
    end if
   
    cc_ptr => characteristic_curves_array(sat_func_id(ghosted_id))%ptr
    
    !assign default end point values for saturation away from transition zones
    sw_min = 0.0d0
    so_min = 0.0d0 !if there is no oil remains zero
    sg_min = 0.0d0 !if there is no gas remains zero
    sw_max = 1.0d0
    sg_max = 1.0d0
    sw_cell = 0.0d0
    so_cell = 0.0d0
    sg_cell = 0.0d0
    
    !minimum saturations computed if rel perm function available
    !if not rel perm => phase not present in the problem and s_alpha_min = 0
    !if (wat_present ) then
    if (associated(cc_ptr%wat_rel_perm_func_owg)) then
      sw_min = cc_ptr%wat_rel_perm_func_owg%GetConnateSaturation(option)
      sw_max = 1.0
    end if
    !end if

    !if (oil_present) then
    if ( associated(cc_ptr%ow_rel_perm_func_owg) ) then
      so_min =  cc_ptr%ow_rel_perm_func_owg%GetConnateSaturation()
    else if ( associated( cc_ptr%oil_rel_perm_func_owg ) ) then
      so_min =  cc_ptr%oil_rel_perm_func_owg%GetConnateSaturation()
    end if
    !end if
      
    !if (gas_present) then
    if ( associated(cc_ptr%gas_rel_perm_func_owg) ) then
      sg_min = cc_ptr%gas_rel_perm_func_owg%GetConnateSaturation(option)
    end if  
    !end if

    sw_max = 1.0 - so_min - sg_min
    sg_max = 1.0 - sw_min - so_min
   
    ! if ( wat_present .and. oil_present .and. gas_present ) then
    !   sw_max = 1.0 - so_min - sg_min
    !   sg_max = 1.0 - sw_min - so_min
    ! else if (wat_present .and. oil_present .and. ( .not.gas_present) ) then
    !   sw_max = 1.0 - so_min
    ! else if (oil_present .and. gas_present .and. ( .not.wat_present) ) then
    !   sg_max = 1.0 - so_min
    ! else if (wat_present .and. gas_present .and. ( .not.oil_present) ) then
    !   sw_max = 1.0 - sg_min
    !   sg_max = 1.0 - sw_min
    ! end if
    
    if (oil_wat_zone) then
      pow_cell = po_cell - pw_cell
      pcow_min = cc_ptr%oil_wat_sat_func%GetPcMin()
      pcow_max = cc_ptr%oil_wat_sat_func%GetPcMax()
      if ( pow_cell <= pcow_min ) then
        sw_cell = sw_max
        po_cell = pw_cell + pcow_min !follows Pw as water supports the pressure
      else if ( pow_cell > pcow_min .and. pow_cell < pcow_max )  then
        if ( .not.cc_ptr%oil_wat_sat_func%sat_func_of_pc_available ) then
           option%io_buffer = 'The Pcow function used for&
                              & hydrostatic equilibration is not valid.'
           call printErrMsg(option)
        end if
        call cc_ptr%oil_wat_sat_func%Saturation(pow_cell,sw_cell, &
                                                dsat_dp,option)
      else if ( pow_cell >= pcow_max ) then
        sw_cell = sw_min
      end if  
    end if

    if (oil_gas_zone) then
      pog_cell = pg_cell - po_cell
      pcog_min = cc_ptr%oil_gas_sat_func%GetPcMin()
      pcog_max = cc_ptr%oil_gas_sat_func%GetPcMax()
      if ( pog_cell <= pcog_min ) then
        sg_cell = sg_min
      else if( pog_cell > pcog_min .and. pog_cell < pcog_max ) then
        if ( .not.cc_ptr%oil_gas_sat_func%sat_func_of_pc_available ) then
          option%io_buffer = 'The Pcog function used for&
                             & hydrostatic equilibration is not valid.'
          call printErrMsg(option)
        end if  
        call cc_ptr%oil_gas_sat_func%Saturation(pog_cell,sg_cell, &
                                                dsat_dp,option)
      else if (pog_cell >= pcog_max ) then
        sg_cell = sg_max
        po_cell = pg_cell - pcog_max !follows the Pg as the gas supports pressure 
      end if    
    end if

    !currenty not used in TOWG as gas condensate is not supported
    !to be tested when gas condensate is supported
    if (wat_gas_zone) then
      pwg_cell = pg_cell - pw_cell
      pcwg_min = cc_ptr%gas_wat_sat_func%GetPcMin()
      pcwg_max = cc_ptr%gas_wat_sat_func%GetPcMax()
      if ( pwg_cell <= pcwg_min ) then
        sw_cell = sw_max
        pg_cell = pw_cell + pcwg_min !pressure supported by water
      else if ( pwg_cell > pcwg_min .and. pwg_cell < pcwg_max ) then
        !to finsh
        if ( .not.cc_ptr%gas_wat_sat_func%sat_func_of_pc_available ) then
          option%io_buffer = 'The Pcgw function used for&
                             & hydrostatic equilibration is not valid.'
          call printErrMsg(option)
        end if
        call cc_ptr%gas_wat_sat_func%Saturation(pwg_cell,sw_cell, &
                                                dsat_dp,option)
        sg_cell = 1.0 - sw_cell
      else if ( pwg_cell >= pcwg_max ) then
        sw_cell = sw_min
        pw_cell = pg_cell - pcwg_max !pressure supported by gas
      end if
    end if

    if ( wat_present .and. (.not.oil_present) .and. (.not.gas_present) ) then
      sw_cell = sw_max
      so_cell = so_min
      sg_cell = sg_min
      po_cell = pw_cell
      pg_cell = pw_cell
    else if (oil_present .and. (.not.wat_present).and.(.not.gas_present)) then
      so_cell =  1.0 - sw_min - sg_min
      sw_cell = sw_min
      sg_cell = sg_min
      pw_cell = po_cell
      pg_cell = po_cell
    else if (gas_present .and. (.not.oil_present).and.(.not.wat_present)) then
      sg_cell = sg_max
      so_cell = so_min
      sw_cell = sw_min
      po_cell = pg_cell
      pw_cell = pg_cell
    else if ( wat_present .and. oil_present .and. (.not.gas_present) ) then
      so_cell =  1.0 - sw_cell
      sg_cell = sg_min
      pg_cell = po_cell
    else if ( oil_present .and. gas_present .and. (.not.wat_present) ) then
      so_cell =  1.0 - sg_cell
      sw_cell =  sw_min
      pw_cell = po_cell
    else if ( wat_present .and. gas_present .and. (.not.oil_present) ) then
      sg_cell = 1.0 - sw_cell
      so_cell = so_min
      po_cell = pg_cell !should not be needed
    else if ( wat_present .and. oil_present .and. gas_present ) then
      so_cell =  1.0 - sw_cell - sg_cell
      if (so_cell < 0.0 ) then
        so_cell = 0.0d0
        pwg_cell = pg_cell - pw_cell
        pcwg_min = pcow_min + pcog_min
        pcwg_max = pcow_max + pcog_max
        if ( pwg_cell <= pcwg_min ) then
          sw_cell = 1.0 - sg_min
          !pg_cell = pw_cell + pcwg(sg) = pw_cell + Pcog(1-Sw) + Pcow(Sw)
          !assumes Po = Pg to assign a value to the oil pressure
          !po_cell = pw_cell + Pcow(Sw)
          call cc_ptr%oil_wat_sat_func%CapillaryPressure(sw_cell,pcow,&
                                                         dp_dsat,option)
          po_cell = pw_cell + pcow
        else if ( pwg_cell > pcwg_min .and. pwg_cell < pcwg_max ) then
          !call Non-linear solution of: Pcwg(Sw) = Pcog(1-Sw) + Pcow(Sw) = pwg_cell
          option%io_buffer = 'hydrostatic: Sw + Sg > 1 not supported yet'
          call printErrMsg(option)
        else if ( pwg_cell >= pcwg_max ) then
          sw_cell = sw_min
          sg_cell = 1.0 - sw_min
          !pw_cell = pg_cell - pcwg(sg) = pg_cell - Pcog(1-Sw) - Pcow(Sw)
          !assumes Po = Pw to assign a value to the oil pressure
          !po_cell = pg_cell - pcwg(sg) = po_cell - Pcog(1-Sw)
          call cc_ptr%oil_gas_sat_func%CapillaryPressure(sg_cell,pcog, &
                                                          dp_dsat,option)
          po_cell = pg_cell - pcog
        end if
      end if  
    end if

    call HydrostaticPMCellInit(option,z_cell,pw_cell,po_cell,pg_cell,sw_cell, &
                               so_cell,sg_cell,iconn,coupler)

  end do !end of connection loop

  call DeallocateGridArrays()

  if (.not.temp_grad_given) then
    nullify(rtempvz_table)
  else
    call LookupTableDestroy(rtempvz_table)
  end if
  nullify(pbvz_table)

end subroutine HydrostaticMPUpdateCoupler

! ************************************************************************** !

subroutine HydrostaticMPInit()
  ! 
  ! Initialise Multi Phase hydrostatic variables
  !
  ! Author: Paolo Orsini
  ! Date: 01/13/19

  implicit none

  wat_gas_equil = PETSC_FALSE
  oil_wat_zone = PETSC_FALSE
  oil_gas_zone = PETSC_FALSE
  wat_gas_zone = PETSC_FALSE
  wat_present = PETSC_FALSE
  oil_present = PETSC_FALSE
  gas_present = PETSC_FALSE  

  !assign default values to datum and pahse contacts
  datum_z = 0.0d0
  owc_z = 0.0d0
  pcow_owc = 0.0d0
  ogc_z = 0.0d0
  pcog_ogc = 0.0d0
  wgc_z = 0.0d0
  pcwg_wgc = 0.0d0
  z_min = 0.0d0
  z_max = 0.0d0

  press_at_datum = UNINITIALIZED_DOUBLE
  
  res_temp_const = UNINITIALIZED_DOUBLE
  res_temp_is_const = PETSC_FALSE

  pb_const = UNINITIALIZED_DOUBLE
  pb_is_const = PETSC_FALSE

  g_z = UNINITIALIZED_DOUBLE

  pw_owc = UNINITIALIZED_DOUBLE
  po_owc = UNINITIALIZED_DOUBLE
  po_ogc = UNINITIALIZED_DOUBLE
  pg_ogc = UNINITIALIZED_DOUBLE
  pw_wgc = UNINITIALIZED_DOUBLE

  pg_wgc = UNINITIALIZED_DOUBLE
  nullify(rtempvz_table)
  temp_grad_given = PETSC_FALSE
  temp_grad(1:3) = UNINITIALIZED_DOUBLE
  temp_at_datum = UNINITIALIZED_DOUBLE
  nullify(pbvz_table)

end subroutine HydrostaticMPInit

! ************************************************************************** !

subroutine HydrostaticPMLoader(condition,option)
!subroutine HydrostaticPMLoader(condition,option,datum_z,press_at_datum, &
!                               owc_z,pcow_owc,ogc_z,pcog_ogc,wgc_z,pcwg_wgc)
  ! 
  ! Extract input data fro equilibration from different modes
  !
  ! Author: Paolo Orsini
  ! Date: 01/07/19

  use Option_module
  use Condition_module

  implicit none

  type(flow_condition_type), intent(in) :: condition
  type(option_type) :: option
  !consider to declare all these input variables as private module variables
  !declare with the same names in the calling function to make this change easier
  !once the old TOIHydrostaticUpdateCoupler has been removed
  ! PetscReal, intent(out) :: datum_z
  ! PetscReal, intent(out) :: press_at_datum
  ! PetscReal, intent(out) :: owc_z
  ! PetscReal, intent(out) :: pcow_owc
  ! PetscReal, intent(out) :: ogc_z
  ! PetscReal, intent(out) :: pcog_ogc
  ! PetscReal, intent(out) :: wgc_z
  ! PetscReal, intent(out) :: pcwg_wgc
  !add temp_table and pb_table

  PetscBool :: press_at_datum_found
  PetscReal :: temperature
  !PetscReal :: temp_grad(3)
  !PetscInt :: tbl_size
  !PetscInt :: data_idx
  !character(len=MAXWORDLENGTH) :: units_word

  !assign default values
  press_at_datum_found = PETSC_FALSE
  temperature = UNINITIALIZED_DOUBLE
  temp_grad(1:3) = UNINITIALIZED_DOUBLE
  temp_grad_given = PETSC_FALSE
  ! datum_z = 0.0d0
  ! press_at_datum = 0.0d0
  ! owc_z = 0.0d0
  ! pcow_owc = 0.0d0
  ! ogc_z = 0.0d0
  ! pcog_ogc = 0.0d0
  ! wgc_z = 0.0d0
  ! pcwg_wgc = 0.0d0


  select case(option%iflowmode)
    case(TOWG_MODE)
      if ( associated(condition%towg%oil_pressure) ) then
        press_at_datum = condition%towg%oil_pressure%dataset%rarray(1)
        press_at_datum_found = PETSC_TRUE
      end if
      !check if datum defined twice then load value 

      !extract phase contact information
      if ( associated(condition%towg%owc_z) ) then
        owc_z = condition%towg%owc_z%dataset%rarray(1)
      end if
      if ( associated(condition%towg%pcow_owc) ) then
        pcow_owc = condition%towg%pcow_owc%dataset%rarray(1)
      end if
      if ( associated(condition%towg%ogc_z) ) then
        ogc_z = condition%towg%ogc_z%dataset%rarray(1)
      end if
      if (option%iflow_sub_mode == TOWG_TODD_LONGSTAFF) then
        !no solvent allowed in the initial solution - put ogc above domain
        ogc_z = 1.0d20
      end if
      if ( associated(condition%towg%pcog_ogc) ) then
        pcog_ogc = condition%towg%pcog_ogc%dataset%rarray(1)
      end if
      if ( associated(condition%towg%temperature) ) then
        !gradient or constant 
        temperature = condition%towg%temperature%dataset%rarray(1)
        if ( associated(condition%towg%temperature%gradient ) ) then
          temp_grad(1:3) = condition%towg%temperature%gradient%rarray(1:3)
          temp_grad_given = PETSC_TRUE
        end if
      end if

      select case (option%iflow_sub_mode)
        case(TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
          if ( associated(condition%towg%pbvz_table) ) then
            pbvz_table => condition%towg%pbvz_table
            pb_is_const = PETSC_FALSE
          else
            pb_const = press_at_datum
            pb_is_const = PETSC_TRUE
          end if
        case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF)
          !nothing to do pb not needed  
      end select
      
      if (condition%towg%is_wg_equilibration) then
        ogc_z = 1.0d20
        owc_z = -1.0d20
        wat_gas_equil = PETSC_TRUE
        !the 
        if ( associated(condition%towg%owc_z) ) then
          wgc_z = condition%towg%owc_z%dataset%rarray(1)
        end if
        if ( associated(condition%towg%pcow_owc) ) then
          pcwg_wgc = condition%towg%pcow_owc%dataset%rarray(1)
        end if            
      end if
      
    case(TOIL_IMS_MODE)
      ogc_z = 1.0d20 !ensure there is no gas
      if ( associated(condition%toil_ims%pressure) ) then
        press_at_datum = condition%toil_ims%pressure%dataset%rarray(1)
        press_at_datum_found = PETSC_TRUE
      end if
      if ( associated(condition%toil_ims%owc) .and. & 
           associated(condition%toil_ims%owc_z) ) then
        option%io_buffer = 'OWC defined twice in owc(1:3) and owc_z'
        call printErrMsg(option)      
      else if ( associated(condition%toil_ims%owc) ) then
         owc_z = condition%toil_ims%owc%dataset%rarray(Z_DIRECTION) 
      else if ( associated(condition%toil_ims%owc_z) ) then
        owc_z = condition%toil_ims%owc_z%dataset%rarray(1)
      end if  
      if ( associated(condition%toil_ims%pcow_owc) ) then
        pcow_owc = condition%toil_ims%pcow_owc%dataset%rarray(1)
      end if
      if ( associated(condition%toil_ims%temperature) ) then
        !gradient or constant 
        temperature = condition%toil_ims%temperature%dataset%rarray(1)
        if ( associated(condition%toil_ims%temperature%gradient ) ) then
          temp_grad(1:3) = condition%toil_ims%temperature%gradient%rarray(1:3)
          temp_grad_given = PETSC_TRUE
        end if
      end if
  end select

  if ( associated(condition%datum_z) .and. &
       associated(condition%datum) &
      ) then
    option%io_buffer = 'TOWG datum defined twice in datum and datum_z'
    call printErrMsg(option)
  else if ( associated(condition%datum_z) ) then
    datum_z = condition%datum_z%dataset%rarray(1)
  else if (associated(condition%datum)) then
    datum_z = condition%datum%rarray(Z_DIRECTION)
  end if

  select case(option%iflowmode)
    case(TOWG_MODE)
      select case (option%iflow_sub_mode)
        case(TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
          if ( .not.associated(pbvz_table) ) then
            if ( dabs(ogc_z - datum_z) > z_eps ) then
              option%io_buffer = 'Hydrostatic equilibration: PB table not defined, &
                                  &datum and OGC must be in the same location'
              call printErrMsg(option)
            end if  
          end if
      end select    
  end select

  if (.not.press_at_datum_found) then
    option%io_buffer = 'TOWG Equilibration condition input error: &
                        &a pressure value at datum must be input'
    call printErrMsg(option)
  end if
  
  if ( associated(condition%rtempvz_table) ) then
    rtempvz_table => condition%rtempvz_table
  end if
  
  !ensure that at least a constant temperature value is setup
  if ( .not.associated(rtempvz_table) .and. (.not.temp_grad_given) ) then
     if (Uninitialized(temperature)) then
       option%io_buffer = 'MP Equilibration input error: a temperature value &
                           &or table must be input to initialise the reservoir'
       call printErrMsg(option)
     end if
     res_temp_is_const = PETSC_TRUE
     res_temp_const = temperature
  end if
  
  !avoid that both a gradient and a table are defined for the temperature
  if ( associated(rtempvz_table) .and. temp_grad_given ) then
    option%io_buffer = 'MP Equilibration input error: define either  &
                        &a table or a gradient for the reservoir temperature'
    call printErrMsg(option)
  end if

  if ( associated(rtempvz_table) .and. Initialized(temperature) ) then
    option%io_buffer = 'MP Equilibration input warning: a constant temperature &
                    &and a table have been defined, the constant temperature &
                    &will be ignored'
    call printMsg(option)
  end if

  !if both temp_grad and rtempvz_table given, returns an error
  if (temp_grad_given) then
    if (Uninitialized(temperature)) then
      option%io_buffer = 'MP Equilibration input error: a temperature value &
                          &at datum must be input to initialise the rservoir &
                          &with a temperature gradient'
      call printErrMsg(option)
    end if
    temp_at_datum = temperature
    ! rtempvz_table => LookupTableCreateGeneral(ONE_INTEGER)
    ! call rtempvz_table%LookupTableVarsInit(TWO_INTEGER)
    ! data_idx = 1 !elevation/depth in the first column
    ! units_word = 'm'
    ! call rtempvz_table%CreateAddLookupTableVar(ONE_INTEGER, &
    !                                  units_word,units_word,data_idx,option)
    ! data_idx = 2 !temperature in the second column
    ! units_word = 'C'
    ! call rtempvz_table%CreateAddLookupTableVar(TWO_INTEGER, &
    !                                  units_word,units_word,data_idx,option)
    ! allocate(rtempvz_table%var_data(2,2))
    ! rtempvz_table%var_data(1,1) = z_min
    ! rtempvz_table%var_data(2,1) = (z_min - datum_z) * temp_grad(3) + &
    !                               temperature
    ! rtempvz_table%var_data(1,2) = z_max
    ! rtempvz_table%var_data(2,2) = (z_max - datum_z) * temp_grad(3) + &
    !                               temperature
    ! !
    ! call rtempvz_table%LookupTableVarConvFactors(option)
    ! call rtempvz_table%VarPointAndUnitConv(option)
    ! call rtempvz_table%SetupConstValExtrap(option)
    ! call rtempvz_table%LookupTableVarInitGradients(option)
    ! allocate(rtempvz_table%axis1)
    ! allocate(rtempvz_table%axis1%values(2))
    ! rtempvz_table%axis1%values(1:2) = rtempvz_table%var_data(1,1:2)
    ! rtempvz_table%dims(1) = 2
  end if

end subroutine HydrostaticPMLoader

! ************************************************************************** !

subroutine RTempTableFromGrad(option)
  
  use Option_module

  implicit none

  type(option_type) :: option

  PetscInt :: data_idx
  character(len=MAXWORDLENGTH) :: units_word

  rtempvz_table => LookupTableCreateGeneral(ONE_INTEGER)
  call rtempvz_table%LookupTableVarsInit(TWO_INTEGER)
  data_idx = 1 !elevation/depth in the first column
  units_word = 'm'
  call rtempvz_table%CreateAddLookupTableVar(ONE_INTEGER, &
                                   units_word,units_word,data_idx,option)
  data_idx = 2 !temperature in the second column
  units_word = 'C'
  call rtempvz_table%CreateAddLookupTableVar(TWO_INTEGER, &
                                   units_word,units_word,data_idx,option)
  allocate(rtempvz_table%var_data(2,2))
  rtempvz_table%var_data(1,1) = z_min
  rtempvz_table%var_data(2,1) = (z_min - datum_z) * temp_grad(3) + &
                                temp_at_datum
  rtempvz_table%var_data(1,2) = z_max
  rtempvz_table%var_data(2,2) = (z_max - datum_z) * temp_grad(3) + &
                                temp_at_datum
  !
  call rtempvz_table%LookupTableVarConvFactors(option)
  call rtempvz_table%VarPointAndUnitConv(option)
  call rtempvz_table%SetupConstValExtrap(option)
  call rtempvz_table%LookupTableVarInitGradients(option)
  allocate(rtempvz_table%axis1)
  allocate(rtempvz_table%axis1%values(2))
  rtempvz_table%axis1%values(1:2) = rtempvz_table%var_data(1,1:2)
  rtempvz_table%dims(1) = 2

end subroutine RTempTableFromGrad

! ************************************************************************** !

subroutine HydrostaticPMCellInit(option,z,pw,po,pg,sw,so,sg,iconn,coupler)
  
  use Option_module
  use Coupler_module
  use PM_TOWG_Aux_module
  use PM_TOilIms_Aux_module
    
  implicit none
  
  type(option_type) :: option
  PetscReal, intent(in) :: z
  PetscReal, intent(in) :: pw
  PetscReal, intent(in) :: po
  PetscReal, intent(in) :: pg
  PetscReal, intent(in) :: sw
  PetscReal, intent(in) :: so
  PetscReal, intent(in) :: sg
  PetscInt, intent(in) :: iconn
  type(coupler_type), intent(inout) :: coupler
  
  !PetscReal, intent(out) :: aux_real(:)
  !PetscInt, intent(out) :: aux_int(:)
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: temp,xm_nacl,pb,rv
  PetscInt :: state
  PetscReal, parameter :: eps_oil   = 1.0d-6
  PetscReal, parameter :: eps_gas   = 1.0d-6
  PetscReal, parameter :: epsp      = 1.0d3
  
  call GetZPointProps(z,option,temp,xm_nacl,pb,rv)
  
  select case(option%iflowmode)
    case(TOIL_IMS_MODE)
      coupler%flow_aux_real_var(ONE_INTEGER,iconn) = po
      coupler%flow_aux_mapping(TOIL_IMS_PRESSURE_INDEX) = ONE_INTEGER
      coupler%flow_aux_real_var(TWO_INTEGER,iconn) = so
      coupler%flow_aux_mapping(TOIL_IMS_OIL_SATURATION_INDEX) = TWO_INTEGER
      coupler%flow_aux_real_var(THREE_INTEGER,iconn) = temp
      coupler%flow_aux_mapping(TOIL_IMS_TEMPERATURE_INDEX) = THREE_INTEGER
    case(TOWG_MODE)
      coupler%flow_aux_real_var(ONE_INTEGER,iconn) = po
      coupler%flow_aux_mapping(TOWG_OIL_PRESSURE_INDEX) = ONE_INTEGER
      coupler%flow_aux_real_var(TWO_INTEGER,iconn) = so
      coupler%flow_aux_mapping(TOWG_OIL_SATURATION_INDEX) = TWO_INTEGER
      select case (option%iflow_sub_mode)
        case(TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
          !PO correction when the user enter Pb > Po for oil region
          if( (sg < eps_gas) .and. (pb >= po) .and. (so > eps_oil) ) then
            pb = po - epsp
            write(word,*) z
            option%io_buffer = 'MP Equilibration warning: oil saturated cell &
                  &found below OGC, resetting to undersaturated z = ' // word
            call printMsg(option)
          end if
          ! ---- to go on pm_towg_aux
          ! Put cells into saturated or undersaturated state (one or other in this case)
          state = TOWG_THREE_PHASE_STATE
          if( (sg < eps_gas) .and. (pb < po) .and. (so > eps_oil) ) then
             state=TOWG_LIQ_OIL_STATE
          end if
          !----
          if( state == TOWG_THREE_PHASE_STATE ) then
            coupler%flow_aux_real_var(THREE_INTEGER,iconn) = sg
            coupler%flow_aux_mapping(TOWG_GAS_SATURATION_INDEX) = THREE_INTEGER
          else
            coupler%flow_aux_real_var(THREE_INTEGER,iconn) = pb
            coupler%flow_aux_mapping(TOWG_BUBBLE_POINT_INDEX) = THREE_INTEGER
          endif
          coupler%flow_aux_int_var(TOWG_STATE_INDEX,iconn) = state
        case (TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF)
          coupler%flow_aux_real_var(THREE_INTEGER,iconn) = sg
          coupler%flow_aux_mapping(TOWG_GAS_SATURATION_INDEX) = THREE_INTEGER
          coupler%flow_aux_int_var(TOWG_STATE_INDEX,iconn) = &
                                                       TOWG_THREE_PHASE_STATE
      end select
      if ( option%iflow_sub_mode == TOWG_SOLVENT_TL ) then
        coupler%flow_aux_real_var(FOUR_INTEGER,iconn) = 0.0d0 !imposing zero sovent saturation
        coupler%flow_aux_mapping(TOWG_SOLV_SATURATION_INDEX) = FOUR_INTEGER
        coupler%flow_aux_real_var(FIVE_INTEGER,iconn) = temp
        coupler%flow_aux_mapping(TOWG_TEMPERATURE_INDEX) = FIVE_INTEGER
      else
        coupler%flow_aux_real_var(FOUR_INTEGER,iconn) = temp
        coupler%flow_aux_mapping(TOWG_TEMPERATURE_INDEX) = FOUR_INTEGER
      end if
  end select
  
end subroutine HydrostaticPMCellInit
  
! coupler%flow_aux_real_var(1,iconn) = po_cell      
! call characteristic_curves%oil_wat_sat_func%Saturation( &
!           pc_comp,sat_liq_comp,dsat_dpres,option)         
! coupler%flow_aux_real_var(2,iconn) = 1.0d0 - sat_liq_comp
! if (coupler%flow_aux_real_var(2,iconn) < 1.0d-6 ) &
!    coupler%flow_aux_real_var(2,iconn) = 1.0d-6 
! end if
! end if
! 
! ! assign temperature
! temperature = temperature_at_datum + &
!           temperature_grad(X_DIRECTION)*dist_x + & ! gradient in K/m
!           temperature_grad(Y_DIRECTION)*dist_y + &
!           temperature_grad(Z_DIRECTION)*dist_z 
! coupler%flow_aux_real_var(3,iconn) = temperature

! ************************************************************************** !
! 
! subroutine InterpolatePMInput(condition,option,one_d_grid)
! !subroutine InterpolatePMInput(condition,option,one_d_grid,press_at_datum, &
! !                             temp_vec,xm_nacl_vec,pb_vec,rs_vec)
!   ! 
!   ! Interpolate thermodynamics properties into 1d grid
!   !
!   ! Author: Paolo Orsini
!   ! Date: 01/09/19  
! 
!   use Option_module
!   use Condition_module
! 
!   implicit none
! 
!   type(flow_condition_type), intent(in) :: condition
!   type(option_type) :: option
!   class(one_dim_grid_type), intent(in) :: one_d_grid
!   ! PetscReal, intent(in) :: press_at_datum
!   ! PetscReal, intent(out) :: temp_vec(:)
!   ! PetscReal, intent(out) :: xm_nacl_vec(:)
!   ! PetscReal, intent(out) :: pb_vec(:)
!   ! PetscReal, intent(out) :: rs_vec(:)
! 
!   PetscReal :: xm_nacl, temperature
!   PetscReal :: temp_grad(3)
!   PetscInt :: grid_size
!   PetscBool :: temp_grad_given
! 
!   xm_nacl = 0.0d0
!   temp_grad = UNINITIALIZED_DOUBLE
!   temperature = UNINITIALIZED_DOUBLE
!   grid_size = size(one_d_grid%z(:))
!   temp_grad_given = PETSC_FALSE
! 
!   select case(option%iflowmode)
!     case(TOWG_MODE)
!       select case (option%iflow_sub_mode)
!         case(TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
!           !if pbvz_table then
!           ! do lookup for each z and copy values onto pb_vec
!           !else          
!           if (pb_is_const) then            
!             !pb_vec(1:grid_size) = press_at_datum
!             pb_vec(1:grid_size) = pb_const
!           end if
!           !end if
!           ! if ( associated(condition%towg%temperature) ) then
!           !   !gradient or constant 
!           !   temperature = condition%towg%temperature%dataset%rarray(1)
!           !   if ( associated(condition%towg%temperature%gradient ) ) then
!           !     temp_grad(1:3) = condition%towg%temperature%gradient%rarray(1:3)
!           !     temp_grad_given = PETSC_TRUE
!           !   end if  
!           ! end if
!         case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF)
!           !do nothing : pb is not required
!       end select
!     case(TOIL_IMS_MODE)
!       ! if ( associated(condition%toil_ims%temperature) ) then
!       !   !gradient or constant 
!       !   temperature = condition%toil_ims%temperature%dataset%rarray(1)
!       !   if ( associated(condition%toil_ims%temperature%gradient ) ) then
!       !     temp_grad(1:3) = condition%toil_ims%temperature%gradient%rarray(1:3)
!       !     temp_grad_given = PETSC_TRUE
!       !   end if  
!       ! end if
!   end select
! 
!   !check also if both temperature and rtempvz_table have been defined and return 
!   !an error. This needs to be implmenented in cocordination with condition read
!   !must allow to read read either temperature (with/wihtout gread) or rtempvz_table
!   ! if not temp, (and temp_grad) or rtempvz_table return an error
! 
!   !if (res_temp_is_const) then
!     temp_vec(1:grid_size) = res_temp_const
!   !else if (associated(rtempvz_table)) then
!   ! interpolation using rtempvz_table
!   !end if
! 
!     ! if (temp_grad_given) then
!     !   call CompVertTempProfile(one_d_grid,temp_grad,temperature, &
!     !                            temp_vec)
!     ! else
!     !   temp_vec(1:grid_size)=temperature
!     ! end if  
!   !end if 
! 
! 
!   !if (associated(condition%saltvz_table)) then
!   ! do lookup for each z and copy values onto xm_nacl_vec
!   !else
!     xm_nacl = option%m_nacl * FMWNACL
!     xm_nacl = xm_nacl /(1.d3 + xm_nacl)
!     xm_nacl_vec(1:grid_size) = xm_nacl
!   !
! 
! end subroutine InterpolatePMInput

! ************************************************************************** !
! 
! subroutine Populate1DGridDataArrays(one_d_grid,option)
!   ! do I need this?
!   use Option_module
! 
!   implicit none
! 
!   class(one_dim_grid_type), intent(in) :: one_d_grid
!   type(option_type) :: option
! 
!   PetscInt :: i_z
! 
!   do i_z = 1,size(one_d_grid%z(:))
!     call GetZPointProps(one_d_grid%z(i_z),option,temp_vec(i_z), &
!                         xm_nacl_vec(i_z),pb_vec(i_z),rv_vec(i_z))
!   end do
! 
! end subroutine Populate1DGridDataArrays
! 
! ************************************************************************** !

function CreateOneDimGrid(min_z,max_z,datum)
  ! 
  ! Computes 1D grid for interpolation needed in hydrostatic pressure
  ! computation
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 

  implicit none

  class(one_dim_grid_type), pointer :: one_d_grid
  PetscReal,intent(in) :: datum(3)
  PetscReal, intent(in) :: min_z, max_z

  class(one_dim_grid_type), pointer :: CreateOneDimGrid

  PetscInt :: num_z, i_z
  PetscReal :: dist_z

  allocate(one_d_grid)

  one_d_grid%min_z = min_z
  one_d_grid%max_z = max_z

  one_d_grid%delta_z = min((max_z-min_z)/500.d0,1.d0)
  ! if zero, assign 1.d0 to avoid divide by zero below. essentially the grid
  ! is flat.
  if (one_d_grid%delta_z < 1.d-40) one_d_grid%delta_z = 1.d0

  num_z = int((max_z-min_z)/one_d_grid%delta_z) + 1

  allocate(one_d_grid%z(num_z))

  one_d_grid%idatum = int((datum(Z_DIRECTION)-min_z)/(max_z-min_z) * &
                        dble(num_z))+1  

  one_d_grid%z(one_d_grid%idatum) = datum(Z_DIRECTION)  

  dist_z = 0.d0
  do i_z=one_d_grid%idatum+1,num_z
    dist_z = dist_z + one_d_grid%delta_z
    one_d_grid%z(i_z) = one_d_grid%z(one_d_grid%idatum) + dist_z
  enddo

  dist_z = 0.d0
  do i_z = one_d_grid%idatum-1,1,-1
    dist_z = dist_z + one_d_grid%delta_z
    one_d_grid%z(i_z) = one_d_grid%z(one_d_grid%idatum) - dist_z
  enddo

  CreateOneDimGrid => one_d_grid

end function CreateOneDimGrid

! ************************************************************************** !

subroutine ZLookup(this,z,i_node,dist_dn,dist_up)
  !
  ! Lookup z within the one dimentional grid and renturns:
  ! i_node : index of the lower closest nodes
  ! dist_z_i_node: distance from the closest lower node
  !
  ! Author: Paolo Orsini
  ! Date: 01/12/19
  !

  implicit none

  class(one_dim_grid_type) :: this
  PetscReal, intent(in) :: z
  PetscInt, intent(out) :: i_node
  PetscReal, intent(out) :: dist_dn
  PetscReal, intent(out) :: dist_up

  PetscReal :: dist_z

  dist_z = z - this%z(this%idatum)
  i_node = this%idatum+int(dist_z/this%delta_z)
  dist_dn = z - this%z(i_node)
  dist_up = this%z(i_node + 1) - z

end subroutine ZLookup

! ************************************************************************** !

subroutine DestroyOneDimGrid(one_d_grid)
  !
  ! destroys OneDimGrid
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/30/15
  ! 

  use Utility_module 

  implicit none

  class(one_dim_grid_type), pointer :: one_d_grid

  if (.not.associated(one_d_grid) ) return

  call DeallocateArray(one_d_grid%z)

  deallocate(one_d_grid)
  nullify(one_d_grid)

end subroutine DestroyOneDimGrid

! ************************************************************************** !

subroutine GetZPointProps(z,option,temp,xm_nacl,pb,rv)
  
  use Option_module
  
  implicit none
  
  PetscReal, intent(in) :: z
  type(option_type) :: option
  PetscReal, intent(out) :: temp
  PetscReal, intent(out) :: xm_nacl
  PetscReal, intent(out) :: pb
  PetscReal, intent(out) :: rv

  temp = UNINITIALIZED_DOUBLE
  xm_nacl = UNINITIALIZED_DOUBLE
  pb = UNINITIALIZED_DOUBLE
  rv = UNINITIALIZED_DOUBLE

  if (res_temp_is_const) then
    temp = res_temp_const
  else
    ! inteprolate using rtempvz_table
    rtempvz_table%axis1%saved_index = 1
    call rtempvz_table%SampleAndGradient(TWO_INTEGER,z)
    temp = rtempvz_table%var_array(TWO_INTEGER)%ptr%sample
  end if
  
  xm_nacl = option%m_nacl * FMWNACL
  xm_nacl = xm_nacl /(1.d3 + xm_nacl)

  select case(option%iflowmode)
    case(TOWG_MODE)
      select case (option%iflow_sub_mode)
        case(TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
          if (pb_is_const) then
            pb = pb_const
          else
            call pbvz_table%SampleAndGradient(TWO_INTEGER,z)
            pb = pbvz_table%var_array(TWO_INTEGER)%ptr%sample
          end if
        case (TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF)
          ! do nothng as pb is not required
      end select
    case(TOIL_IMS_MODE)
  end select

end subroutine GetZPointProps

! ************************************************************************** !

subroutine PhaseHydrostaticPressure(one_d_grid,option,iphase,press_start, &
                                    id_start,temp_v,xm_nacl_v,pb_v,rv_v, &
                                    press_v,den_kg_v)
  ! 
  ! Compute an hydrostatic pressure profile for a given phase and 1D
  ! 1D discretisation
  ! The "starting point" for the computation of the pressure profile
  ! is the elevation where a reference pressure is given, it can be the datum,
  ! or one of the phase contact interface (e.g. OWC, OGC, etc)
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15 - 01/10/2019
  !

  use Option_module

  implicit none

  class(one_dim_grid_type) :: one_d_grid
  !PetscReal, intent(in) :: gravity(:) ! this is option%gravity
  type(option_type) :: option
  PetscInt, intent(in) :: iphase
  PetscReal, intent(in) :: press_start
  PetscInt, intent(in) :: id_start
  PetscReal, intent(in) :: temp_v(:)
  PetscReal, intent(in) :: xm_nacl_v(:)
  PetscReal, intent(in) :: pb_v(:)
  PetscReal, intent(in) :: rv_v(:)
  PetscReal, intent(out) :: press_v(:)
  PetscReal, intent(out) :: den_kg_v(:)

  PetscReal :: gravity(3)
  PetscReal :: pressure, pressure0, rho, rho_one, rho_kg, rho_zero
  PetscInt :: ipressure, num_iteration

  gravity(1:3) = option%gravity(1:3)

  rho_kg = PhaseDensity(option,iphase,press_start,temp_v(id_start), &
                        xm_nacl_v(id_start),pb_v(id_start),rv_v(id_start))

  ! fill properties for reference pressure
  den_kg_v(id_start) = rho_kg
  press_v(id_start) = press_start

  pressure0 = press_start
  rho_zero = rho_kg
  do ipressure=id_start+1,size(one_d_grid%z(:))
    rho_kg = PhaseDensity(option,iphase,pressure0,temp_v(ipressure), &
                          xm_nacl_v(ipressure),pb_v(ipressure),rv_v(ipressure))
    num_iteration = 0
    do
      pressure = pressure0 + 0.5d0*(rho_kg+rho_zero) * &
                 gravity(Z_DIRECTION) * one_d_grid%delta_z
      ! if (iphase == OIL_PHASE) then
      !   select case(option%iflowmode)
      !     case(TOWG_MODE)
      !       select case (option%iflow_sub_mode)
      !         case(TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
      !           ! if ( one_d_grid%z(ipressure) >= ogc_z) then
      !           !   pressure = pb_v(ipressure)
      !           ! else
      !           !   if (pressure < pb_v(ipressure)) then
      !           !     pressure = pb_v(ipressure)
      !           !   end if
      !           ! end if
      !       end select
      !   end select
      ! end if !end if OIL_PHASE
      rho_one = PhaseDensity(option,iphase,pressure,temp_v(ipressure), &
                        xm_nacl_v(ipressure),pb_v(ipressure),rv_v(ipressure))
      !check convergence on density
      if (dabs(rho_kg-rho_one) < 1.d-10) exit
      rho_kg = rho_one
      num_iteration = num_iteration + 1
      if (num_iteration > 100) then
        print *,'Phase-Hydrostatic iteration failed to converge', &
                 num_iteration,rho_one,rho_kg
        !print *, condition%name, idatum
        !print *, pressure_array
        stop
      endif
    enddo
    rho_zero = rho_kg
    press_v(ipressure) = pressure
    den_kg_v(ipressure) = rho_kg
    pressure0 = pressure
  enddo

  ! compute pressures below one_d_grid%z(id_start), if any
  pressure0 = press_v(id_start)
  rho_zero = den_kg_v(id_start)
  do ipressure=id_start-1,1,-1
    rho_kg = PhaseDensity(option,iphase,pressure0,temp_v(ipressure), &
                          xm_nacl_v(ipressure),pb_v(ipressure),rv_v(ipressure))
    num_iteration = 0
    do                   ! notice the negative sign (-) here
      pressure = pressure0 - 0.5d0*(rho_kg+rho_zero) * &
                 gravity(Z_DIRECTION) * one_d_grid%delta_z
      ! if (iphase == OIL_PHASE) then
      !   select case(option%iflowmode)
      !     case(TOWG_MODE)
      !       select case (option%iflow_sub_mode)
      !         case(TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
      !           if ( one_d_grid%z(ipressure) >= ogc_z) then
      !             pressure = pb_v(ipressure)
      !           else
      !             if (pressure < pb_v(ipressure)) then
      !               pressure = pb_v(ipressure)
      !             end if
      !           end if
      !       end select
      !   end select
      ! end if !end if OIL_PHASE
      rho_one = PhaseDensity(option,iphase,pressure,temp_v(ipressure), &
                         xm_nacl_v(ipressure),pb_v(ipressure),rv_v(ipressure))
      !check convergence on density
      if (dabs(rho_kg-rho_one) < 1.d-10) exit
      rho_kg = rho_one
      num_iteration = num_iteration + 1
      if (num_iteration > 100) then
        print *,'Phase-Hydrostatic iteration failed to converge', &
                 num_iteration,rho_one,rho_kg
        !print *, condition%name, idatum
        !print *, pressure_array
        stop
      endif
    enddo
    rho_zero = rho_kg
    press_v(ipressure) = pressure
    den_kg_v(ipressure) = rho_kg
    pressure0 = pressure
  enddo

  !write(*,*) "den fine bottom", den_kg(1) 
  !write(*,*) "den fine top", den_kg(size(one_d_grid%z(:)))

end subroutine PhaseHydrostaticPressure

! ************************************************************************** !
function PhaseDensity(option,iphase,p,t,xm_nacl,pb,rv)
  !
  ! computes phase density given the specified phase
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15 - 01/10/2019
  ! 

  use Option_module
  use EOS_Oil_module
  use EOS_Water_module
  use EOS_Gas_module ! when gas is considered

  implicit none

  type(option_type) :: option
  PetscInt, intent(in) :: iphase
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: xm_nacl
  PetscReal, intent(in) :: pb
  PetscReal, intent(in) :: rv
  PetscReal :: PhaseDensity ! kg/m3 

  PetscReal :: dw_mol
  PetscReal :: aux(1)
  PetscReal :: rs_molar
  PetscReal :: xo, xg, cr, deno, crusp
  PetscReal :: pb_comp
  !PetscInt :: table_idx

  PetscErrorCode :: ierr

  !table_idx = 1
  xo = 0.0d0
  xg = 0.0d0
  cr = 0.0d0
  deno = 0.0d0
  crusp = 0.0d0

  select case(iphase)
    case(LIQUID_PHASE)
      if ( option%m_nacl < 1.0d-40) then !no salt in water
        call EOSWaterDensity(t,p,PhaseDensity,dw_mol,ierr)
      else
        aux(1) = xm_nacl
        call EOSWaterDensityExt(t,p,aux,PhaseDensity,dw_mol,ierr)
      end if
    case(GAS_PHASE)
      !note: curently no gas-condensate can be modeled
      call EOSGasDensity(t,p,PhaseDensity,ierr)
      !convert to mass density
      PhaseDensity = PhaseDensity * EOSGasGetFMW()
    case(OIL_PHASE)
      select case(option%iflowmode)
        case(TOIL_IMS_MODE)
          call EOSOilDensity(t,p,PhaseDensity,ierr)
          PhaseDensity = PhaseDensity * EOSOilGetFMW()
        case(TOWG_MODE)
          select case (option%iflow_sub_mode)
            case(TOWG_BLACK_OIL,TOWG_SOLVENT_TL)
              pb_comp = min(pb,p)
              !Look up the Rs value
              !table_idx =  1
              !call getBlackOilComposition(pb,t,table_idx,xo,xg)
              call EOSOilRS(t,pb_comp,rs_molar,ierr)
              xo=1.0d0   /(1.0d0+rs_molar)
              xg=rs_molar/(1.0d0+rs_molar)
              !Density and compressibility look-up at bubble point
              call EOSOilDensity (t,pb_comp,deno,ierr)
              call EOSOilCompressibility(t,pb_comp,cr,ierr)
              !Correct for undersaturation
              crusp=cr*(p-pb_comp)
              deno=deno*(1.0+crusp*(1.0+0.5*crusp))
              !correct for composition, i.e. add gas component contribution
              deno=deno/xo
              !convert to mass density
              PhaseDensity = deno * (xo*EOSOilGetFMW() + xg*EOSGasGetFMW() )
            case(TOWG_IMMISCIBLE,TOWG_TODD_LONGSTAFF)
              call EOSOilDensity(t,p,PhaseDensity,ierr)
              PhaseDensity = PhaseDensity * EOSOilGetFMW()
        end select
      end select
  end select

end function PhaseDensity

! ************************************************************************** !

subroutine CompVertTempProfile(one_d_grid,temp_grad,temp_at_datum, &
                               temp_profile)
  ! 
  ! Computes the temperature vertical profile on the 1D grid use for the 
  ! hydrostatic pressure calculation
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 

  implicit none

  class(one_dim_grid_type) :: one_d_grid
  PetscReal, intent(in) :: temp_grad(:)
  PetscReal, intent(in) :: temp_at_datum
  PetscReal, intent(out) :: temp_profile(:)
 
  PetscInt :: i_z

  temp_profile(one_d_grid%idatum) = temp_at_datum

  do i_z=one_d_grid%idatum+1,size(one_d_grid%z(:))
    temp_profile(i_z) = temp_profile(i_z-1) + &
                        temp_grad(Z_DIRECTION)*one_d_grid%delta_z   
  end do

  do i_z=one_d_grid%idatum-1,1,-1
    temp_profile(i_z) = temp_profile(i_z+1) - & ! note the (-) sign
                        temp_grad(Z_DIRECTION)*one_d_grid%delta_z   
  end do

end subroutine CompVertTempProfile

! ************************************************************************** !

subroutine AllocateGridArrays(grid_size)

  implicit none

  PetscInt, intent(in) :: grid_size

  allocate(temp_vec(grid_size))
  temp_vec = UNINITIALIZED_DOUBLE
  allocate(xm_nacl_vec(grid_size))
  xm_nacl_vec = 0.0d0
  allocate(pb_vec(grid_size))
  pb_vec = UNINITIALIZED_DOUBLE
  allocate(rv_vec(grid_size))
  rv_vec = UNINITIALIZED_DOUBLE
  allocate(wat_press_vec(grid_size))
  wat_press_vec = UNINITIALIZED_DOUBLE
  allocate(wat_den_kg_vec(grid_size))
  wat_den_kg_vec = UNINITIALIZED_DOUBLE
  allocate(oil_press_vec(grid_size))
  oil_press_vec = UNINITIALIZED_DOUBLE
  allocate(oil_den_kg_vec(grid_size))
  oil_den_kg_vec = UNINITIALIZED_DOUBLE
  allocate(gas_press_vec(grid_size))
  gas_press_vec = UNINITIALIZED_DOUBLE
  allocate(gas_den_kg_vec(grid_size))
  gas_den_kg_vec = UNINITIALIZED_DOUBLE

end subroutine AllocateGridArrays

! ************************************************************************** !

subroutine DeallocateGridArrays()

  implicit none

  deallocate(temp_vec)
  deallocate(xm_nacl_vec)
  deallocate(pb_vec)
  deallocate(rv_vec)
  deallocate(wat_press_vec)
  deallocate(wat_den_kg_vec)
  deallocate(oil_press_vec)
  deallocate(oil_den_kg_vec)
  deallocate(gas_press_vec)
  deallocate(gas_den_kg_vec)

end subroutine DeallocateGridArrays

! ************************************************************************** !

subroutine TOIHydrostaticUpdateCoupler(coupler,option,grid, &
              characteristic_curves_array,sat_func_id,imat)

  ! 
  ! Computes the hydrostatic initial/boundary condition for oil/water systems
  ! given the oil water contact (owc) elevation
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 
  ! Algorithm description (To add with a schematic drawing of OWC)
  !- Identify location of datum and OWC within 1D domain
  ! 
  ! IF DATUM FALLS IN WATER REGION:
  ! 1. COMPUTE WATER HYDROSTATIC PRESSURE from:
  ! - istart (idatum),
  ! - pressure_start (pressure_at_datum), this must be a water pressure
  ! - temperature_at_datum, or temp_array(npressure)
  !     
  ! 2. COMPUTE OIL HYDROSTATIC PRESSURE profile passing:
  ! - istart (OWC_loc_id)
  ! - pressure_start (oil pressure_at_OWC), computed from pw at OWC and pc,  
  !   which is an input or computed from a saturation function for So = Soir
  ! - temperature_at_datum, or temp_array(npressure)

  ! IF DATUM FALLS in OIL REGION
  ! 1. COMPUTE OIL HYDROSTATIC PRESSURE from:
  ! - istart (idatum),
  ! - pressure_start (pressure_at_datum), this must be an oil pressure
  ! - temperature_at_datum, or temp_array(npressure) 
  ! 2. COMPUTES WATER HYDROSTATIC PRESSURE profile passing:
  ! - istart (OWC_loc_id)
  ! - pressure_start (water pressure_at_OWC), computed from po at OWC and pc, 
  !   which is an input or computed from a saturation function for So = Soir
  ! - temperature_at_datum, or temp_array(npressure)

  !The oil and water pressure arrays on the 1D grid are used to interpolate
  ! the pressure and compute the saturation on the 3D grid cells
  ! - compute oil initial phase pressures
  ! - compute equilibrating capillary pressure: pc = po - pw 
  ! - compute saturation from inverse pc curve: sw = pc^(-1)(sw) 

  !PO note - this should be a member function of the coupler class!

  use EOS_Water_module

  use Option_module
  use Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Characteristic_Curves_module
  use Region_module
  !use Grid_Structured_module
  !use Utility_module, only : DotProduct
  use Utility_module
  use Dataset_Gridded_HDF5_class
  use Dataset_Ascii_class
 
  !use TOilIms_Aux_module
  use PM_TOilIms_Aux_module 
    ! to use constant paramters such as TOIL_IMS_PRESSURE_DOF
    ! could work something out to eliminate this dependency here 
  

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  type(characteristic_curves_ptr_type) :: characteristic_curves_array(:)
  PetscInt, pointer, intent(in) :: sat_func_id(:)
  PetscInt, pointer, intent(in) :: imat(:)

  PetscReal :: xm_nacl
  PetscInt :: local_id, ghosted_id, iconn, ghosted_id_min_dist
  PetscReal :: oil_press_grad(3), wat_press_grad(3), temperature_grad(3)  
  PetscReal :: datum(3), owc(3)
  PetscReal :: pressure_at_datum, temperature_at_datum, press_start
  PetscReal :: max_z, min_z
  PetscReal :: gravity_magnitude
  PetscReal :: pw_owc, po_owc, pc_owc, pw_cell, po_cell, temperature, temp_owc
  PetscInt  :: i_owc, ipressure, ipress_start ! id_loc_owc 
  PetscReal, pointer :: wat_pressure_array(:), wat_density_array(:)
  PetscReal, pointer :: oil_pressure_array(:), oil_density_array(:)  
  PetscReal, pointer :: temperature_array(:), xm_nacl_array(:)
  PetscReal, pointer :: dummy_array1(:)
  PetscReal :: dummy_val
  PetscReal :: dist_x, dist_y, dist_z, delta_z, dist_z_for_pressure
  PetscReal :: dist_owc_start
  PetscReal :: dist_z_owc, dx_conn, dy_conn, dz_conn
  PetscBool :: datum_in_water, pw_hydrostatic, po_hydrostatic
  PetscReal :: sat_liq_owc, pc_comp, sat_liq_comp, dsat_dpres
  PetscReal :: sat_ir(2)
  PetscReal :: dpc_dsatl

  class(one_dim_grid_type), pointer :: one_d_grid
  type(flow_condition_type), pointer :: condition
  class(characteristic_curves_type), pointer :: characteristic_curves
  !class(characteristic_curves_type) :: characteristic_curves

  if (coupler%connection_set%num_connections == 0 ) return

  pw_hydrostatic = PETSC_TRUE
  po_hydrostatic = PETSC_TRUE  

  condition => coupler%flow_condition
  
  nullify(wat_pressure_array)
  nullify(wat_density_array)
  nullify(oil_pressure_array)
  nullify(oil_density_array)
  nullify(temperature_array)
  nullify(xm_nacl_array)
  nullify(dummy_array1)
  

  ! fix indices for map to flow_aux_real_var
  coupler%flow_aux_mapping(TOIL_IMS_PRESSURE_INDEX) = 1
  coupler%flow_aux_mapping(TOIL_IMS_OIL_SATURATION_INDEX) = 2
  coupler%flow_aux_mapping(TOIL_IMS_TEMPERATURE_INDEX) = 3 

  xm_nacl = option%m_nacl * FMWNACL
  xm_nacl = xm_nacl /(1.d3 + xm_nacl)

  !initialise datum
  datum(:) = 0.0d0
  oil_press_grad(:) = 0.0d0
  wat_press_grad(:) = 0.0d0
  temperature_grad(:) = 0.0d0
  dummy_val = 0.0d0 
   
  if ( associated(condition%datum) ) then
    datum(1:3) = condition%datum%rarray(1:3) 
  end if   
  pressure_at_datum = &
      condition%toil_ims%pressure%dataset%rarray(1)
  ! gradient is in m/m; needs conversion to Pa/m
  if (associated(condition%toil_ims%pressure%gradient)) then
    oil_press_grad(1:3) = &
    condition%toil_ims%pressure%gradient%rarray(1:3)
  endif
  if (associated(condition%toil_ims%liq_press_grad)) then
    wat_press_grad(1:3) = &
      condition%toil_ims%liq_press_grad%dataset%rarray(1:3)
  endif

  temperature_at_datum = &
    condition%toil_ims%temperature%dataset%rarray(1)
  if (associated(condition%toil_ims%temperature%gradient)) then
     temperature_grad(1:3) = &
       condition%toil_ims%temperature%gradient%rarray(1:3)
  endif

  max_z = max(grid%z_max_global,datum(Z_DIRECTION))+1.d0 ! add 1m buffer
  min_z = min(grid%z_min_global,datum(Z_DIRECTION))-1.d0

  if (associated(condition%toil_ims%owc)) then
    owc(1:3) = condition%toil_ims%owc%dataset%rarray(1:3)
  else ! default condition assumes owc above domain and datum: water domain
    if ( datum(Z_DIRECTION) >= max_z ) then
      owc(Z_DIRECTION) = datum(Z_DIRECTION) + 1.d0
    else
      owc(Z_DIRECTION) = max_z + 1.d0 ! place owc 1 m above domain
    end if
    datum_in_water = PETSC_TRUE
  end if 
  
  ! finds out where datume is located
  if ( datum(Z_DIRECTION) >= owc(Z_DIRECTION) ) then
    datum_in_water = PETSC_FALSE !datum is in the oil region 
  else
    datum_in_water = PETSC_TRUE !datum is in the water region
  end if

  if (dabs(wat_press_grad(Z_DIRECTION)) >= 1.d-40) pw_hydrostatic = PETSC_FALSE
  if (dabs(oil_press_grad(Z_DIRECTION)) >= 1.d-40) po_hydrostatic = PETSC_FALSE

  ! flag off unnecessary hydrostatic computations - ensure owc in one_d_grid%z
  if (owc(Z_DIRECTION) > max_z ) then
    owc(Z_DIRECTION) = max_z
    po_hydrostatic = PETSC_FALSE
  end if
  if (owc(Z_DIRECTION) < min_z )  then
    owc(Z_DIRECTION) = min_z
    pw_hydrostatic = PETSC_FALSE
  end if

  !if one of the phases requires automatic hydrostatic pressure 
  ! the 1D domain and discretisation is needed
  !if ( (oil_press_grad(Z_DIRECTION) < 1.d-40) .or. &
  !     (wat_press_grad(Z_DIRECTION) < 1.d-40 ) &
  !   ) then
  if ( pw_hydrostatic .or. po_hydrostatic ) then

    gravity_magnitude = sqrt(DotProduct(option%gravity,option%gravity))
  
    if (dabs(gravity_magnitude-EARTH_GRAVITY) > 0.1d0) then
      option%io_buffer = 'Magnitude of gravity vector is not near 9.81.'
      call printErrMsg(option)
    endif
 
    ! ceat 1D domain and discretization needed for interpolations
    one_d_grid => CreateOneDimGrid(min_z,max_z,datum)

    allocate(temperature_array(size(one_d_grid%z(:))))
    temperature_array = 0.d0
    call CompVertTempProfile(one_d_grid,temperature_grad, &
                             temperature_at_datum,temperature_array)
    ! allocate pressure, density and temperature arrays
    !allocate(pressure_array(2,size(one_d_grid%z(:)))) 
    !allocate(density_array(2,size(one_d_grid%z(:))))
    if (pw_hydrostatic) then
      allocate(wat_pressure_array(size(one_d_grid%z(:))))
      allocate(wat_density_array(size(one_d_grid%z(:))))
    end if
    if (po_hydrostatic) then
      allocate(oil_pressure_array(size(one_d_grid%z(:))))
      allocate(oil_density_array(size(one_d_grid%z(:))))
    end if
  end if
  allocate(dummy_array1(size(one_d_grid%z(:))))
  dummy_array1 = 0.0d0
  allocate(xm_nacl_array(size(one_d_grid%z(:))))
  xm_nacl_array = xm_nacl
  
  dist_x = 0.d0
  dist_y = 0.d0
  dist_z = owc(Z_DIRECTION) - datum(Z_DIRECTION)
  po_owc = 0.d0
  pw_owc = 0.d0
  ! identify where the owc is located in the one_dim mesh
  if (pw_hydrostatic .or. po_hydrostatic) then 
    i_owc = one_d_grid%idatum+int(dist_z/one_d_grid%delta_z)
    dist_z_for_pressure = owc(Z_DIRECTION) - one_d_grid%z(i_owc)
    temp_owc = temperature_array(i_owc) + &
               dist_z_for_pressure * temperature_grad(Z_DIRECTION)
  end if

  ghosted_id_min_dist = &
      !GetCouplerCellOnPhaseConact(coupler,grid,imat,owc(Z_DIRECTION),option)
       GetCellOnPhaseConact(coupler%connection_set,grid, &
                            imat,owc(Z_DIRECTION),option)

  !write(*,*) "my_rank", option%myrank, "min_ghost", ghosted_id_min_dist, &
  !           "sat_fun_id", sat_func_id(ghosted_id_min_dist)
  !PO This approach must be review - pc at owc must be entered in the input
  characteristic_curves => &
      characteristic_curves_array(sat_func_id(ghosted_id_min_dist))%ptr
  !characteristic_curves = characteristic_curves_array(func_id)%ptr
  ! the OWC is assumed to be located where So = So_ir 
  sat_ir(:) = CharCurvesGetGetResidualSats(characteristic_curves,option) 
  sat_liq_owc = 1.0 - sat_ir(2)
      
  if ( (pw_hydrostatic .or. po_hydrostatic ) .and. &
       (.not.characteristic_curves%oil_wat_sat_func%sat_func_of_pc_available) &
     ) then
       option%io_buffer = 'The capillary pressure function used for&
                          & hydrostatic equilibration is not valid.'
      call printErrMsg(option)
  end if    

  call characteristic_curves%oil_wat_sat_func% &
           CapillaryPressure(sat_liq_owc,pc_owc,dpc_dsatl,option)

  ! compute pressure and density profiles for phases where hydrostatic pressure
  ! is imposed. And pressure (water or oil) at owc elevation
  if (datum_in_water) then 
    if (pw_hydrostatic) then
      ! call PhaseHydrostaticPressure(one_d_grid,option%gravity, &
      !           LIQUID_PHASE,pressure_at_datum, &
      !           one_d_grid%idatum,xm_nacl,temperature_array, &
      !           wat_pressure_array,wat_density_array)
      call PhaseHydrostaticPressure(one_d_grid,option,LIQUID_PHASE, &
                  pressure_at_datum,one_d_grid%idatum,temperature_array, &
                  xm_nacl_array,dummy_array1,dummy_array1, &
                  wat_pressure_array,wat_density_array)
                  
      pw_owc = PressInterp(i_owc,dist_x,dist_y,dist_z_for_pressure, &
                           option%gravity,wat_pressure_array, &
                           wat_density_array,wat_press_grad)
    else
      pw_owc = PressGrad(dist_x,dist_y,dist_z,pressure_at_datum,wat_press_grad) 
    end if
    ! test pc=0
    !po_owc = pw_owc
    po_owc = pw_owc + pc_owc

    if (po_hydrostatic) then
      ! compute oil press and denisty profiles
      !id_loc_owc = one_d_grid%ElevationIdLoc(owc(Z_DIRECTION))
      ipress_start = i_owc + 1
      !if ( ipress_start > size(one_d_grid%z(:)) ) & !need to avoid this computation if owc > z_max
      !  ipress_start = size(one_d_grid%z(:))
      dist_owc_start = one_d_grid%z(ipress_start) - owc(Z_DIRECTION)
      press_start = po_owc + dist_owc_start * &
                    ! PhaseDensity(OIL_PHASE,po_owc,temp_owc,xm_nacl) * &
                    PhaseDensity(option,OIL_PHASE,po_owc,temp_owc,xm_nacl, &
                            dummy_val,dummy_val) * option%gravity(Z_DIRECTION)
      ! call PhaseHydrostaticPressure(one_d_grid,option%gravity, &
      !           OIL_PHASE,press_start, &
      !           ipress_start,xm_nacl,temperature_array, &
      !           oil_pressure_array,oil_density_array)
      !
      call PhaseHydrostaticPressure(one_d_grid,option,OIL_PHASE, &
                  press_start,ipress_start,temperature_array, &
                  xm_nacl_array,dummy_array1,dummy_array1, &
                  oil_pressure_array,oil_density_array)
    end if
  else ! datum is in the oil region
    if (po_hydrostatic) then
      ! call PhaseHydrostaticPressure(one_d_grid,option%gravity, &
      !           OIL_PHASE,pressure_at_datum, &
      !           one_d_grid%idatum,xm_nacl,temperature_array, &
      !           oil_pressure_array,oil_density_array)
      !
      call PhaseHydrostaticPressure(one_d_grid,option,OIL_PHASE, &
                  pressure_at_datum,one_d_grid%idatum,temperature_array, &
                  xm_nacl_array,dummy_array1,dummy_array1, &
                  oil_pressure_array,oil_density_array)      
      po_owc = PressInterp(i_owc,dist_x,dist_y,dist_z_for_pressure, &
                           option%gravity,oil_pressure_array, &
                           oil_density_array,oil_press_grad)
    else
      po_owc = PressGrad(dist_x,dist_y,dist_z,pressure_at_datum,oil_press_grad) 
    end if
    ! testing pc = 0 
    !pw_owc = po_owc
    pw_owc = po_owc - pc_owc

    if (pw_hydrostatic) then
      ! conmpute water press and denisty profiles
      !id_loc_owc = one_d_grid%ElevationIdLoc(owc(Z_DIRECTION))
      ipress_start = i_owc + 1
      dist_owc_start = one_d_grid%z(ipress_start) - owc(Z_DIRECTION)
      press_start = pw_owc + dist_owc_start * &
                    !PhaseDensity(LIQUID_PHASE,pw_owc,temp_owc,xm_nacl) * &
                    PhaseDensity(option,LIQUID_PHASE,pw_owc,temp_owc,xm_nacl, &
                            dummy_val,dummy_val) * option%gravity(Z_DIRECTION)
      ! call PhaseHydrostaticPressure(one_d_grid,option%gravity, &
      !           LIQUID_PHASE,press_start, &
      !           ipress_start,xm_nacl,temperature_array, &
      !           wat_pressure_array,wat_density_array)
      !
      call PhaseHydrostaticPressure(one_d_grid,option,LIQUID_PHASE, &
                  press_start,ipress_start,temperature_array, &
                  xm_nacl_array,dummy_array1,dummy_array1, &
                  wat_pressure_array,wat_density_array)      
    end if
  end if
  
  !compute pressure values in the coupler region cells 

  dx_conn = 0.d0
  dy_conn = 0.d0
  dz_conn = 0.d0

  do iconn=1,coupler%connection_set%num_connections 
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
 
    ! geh: note that this is a boundary connection, thus the entire distance is between
    ! the face and cell center
    if (associated(coupler%connection_set%dist)) then
      dx_conn = coupler%connection_set%dist(0,iconn) * &
                      coupler%connection_set%dist(1,iconn)
      dy_conn = coupler%connection_set%dist(0,iconn) * &
                      coupler%connection_set%dist(2,iconn)
      dz_conn = coupler%connection_set%dist(0,iconn) * &
                      coupler%connection_set%dist(3,iconn)
    endif

    ! note the negative (-) d?_conn is required due to the offset of the boundary face
    dist_x = grid%x(ghosted_id)-dx_conn-datum(X_DIRECTION) !datume here can be datum or owc !!
    dist_y = grid%y(ghosted_id)-dy_conn-datum(Y_DIRECTION)
    dist_z = grid%z(ghosted_id)-dz_conn-datum(Z_DIRECTION)
    ! 
    if ( pw_hydrostatic .or. po_hydrostatic ) then
      !location in the press_arrays
      ipressure = one_d_grid%idatum+int(dist_z/one_d_grid%delta_z) 
      dist_z_for_pressure = grid%z(ghosted_id) - &
                            dz_conn - one_d_grid%z(ipressure) 
    end if
  
    if ( (.not.pw_hydrostatic) .or. (.not.po_hydrostatic) ) then 
      dist_z_owc = grid%z(ghosted_id)-dz_conn-owc(Z_DIRECTION) !z_ref is owc
    end if

    if (pw_hydrostatic) then
      pw_cell = PressInterp(ipressure,dist_x,dist_y,dist_z_for_pressure, &
                            option%gravity,wat_pressure_array, &
                            wat_density_array,wat_press_grad)
    else
      pw_cell = PressGrad(dist_x,dist_y,dist_z_owc,pw_owc,wat_press_grad)
    end if 

    if (po_hydrostatic) then
      po_cell = PressInterp(ipressure,dist_x,dist_y,dist_z_for_pressure, &
                            option%gravity,oil_pressure_array, &
                            oil_density_array,oil_press_grad)
    else
      po_cell = PressGrad(dist_x,dist_y,dist_z_owc,po_owc,oil_press_grad)
    end if 
    
    !COMPUTE HERE SATURATION from pressure difference
    if (owc(Z_DIRECTION) >= max_z) then !owc above domain (water only)
      coupler%flow_aux_real_var(1,iconn) = pw_cell
      coupler%flow_aux_real_var(2,iconn) = &
        coupler%flow_condition%toil_ims%saturation%dataset%rarray(1)
        !to avoid truncation erros
        if (coupler%flow_aux_real_var(2,iconn) < 1.0d-6 ) then
          coupler%flow_aux_real_var(2,iconn) = 1.0d-6
        end if
    else if (owc(Z_DIRECTION) <= min_z) then !owc below domain (oil only)
      !OIL PRESSURE
      coupler%flow_aux_real_var(1,iconn) = po_cell 
      coupler%flow_aux_real_var(2,iconn) = &
                coupler%flow_condition%toil_ims%saturation%dataset%rarray(1)
    else
      !use instgructions below when imposing pc=0 capillary pressure 
      !if ( grid%z(ghosted_id) > owc(Z_DIRECTION) ) then
      !  !OIL PRESSURE
      !  coupler%flow_aux_real_var(1,iconn) = po_cell 
      !  !OIL SATURATION
      !  coupler%flow_aux_real_var(2,iconn) = 1.0d0
      !else 
      !  coupler%flow_aux_real_var(1,iconn) = pw_cell
      !  coupler%flow_aux_real_var(2,iconn) = 1.0d-6 !to avoid truncation erros
      !end if
      !use insytruction below for pc /= 0 
      pc_comp = po_cell - pw_cell
      if ( pc_comp <= 0.d0 ) then ! water-only region 
        coupler%flow_aux_real_var(1,iconn) = pw_cell
        coupler%flow_aux_real_var(2,iconn) = 1.0d-6 !to avoid truncation erros
      !oil region - case of zero capillary pressure - can assign So < So_ir
      else if ( (pc_comp > 0.d0) .and. &
              !  (characteristic_curves%saturation_function%pcmax < 1.0d-40 ) &
               (characteristic_curves%oil_wat_sat_func%pcmax < 1.0d-40 ) &
              ) then  
        !OIL SATURATION from input
        coupler%flow_aux_real_var(1,iconn) = po_cell
        coupler%flow_aux_real_var(2,iconn) = &
                coupler%flow_condition%toil_ims%saturation%dataset%rarray(1)        
      !else if ( pc_comp >= characteristic_curves%saturation_function%pcmax ) &
      else if ( pc_comp >= characteristic_curves%oil_wat_sat_func%pcmax ) &
        then
        ! oil region: can consider here connate water if required, or Sw_ir
        ! sat_ir(1) currenlty used, this Sw_ir taken from owc characteristic curve
        ! should consider input from deck
        coupler%flow_aux_real_var(1,iconn) = po_cell
        coupler%flow_aux_real_var(2,iconn) = 1.0d0 - sat_ir(1)
      else
        ! water/oil transition zone
        coupler%flow_aux_real_var(1,iconn) = po_cell      
        call characteristic_curves%oil_wat_sat_func%Saturation( &
                  pc_comp,sat_liq_comp,dsat_dpres,option)         
        coupler%flow_aux_real_var(2,iconn) = 1.0d0 - sat_liq_comp
        if (coupler%flow_aux_real_var(2,iconn) < 1.0d-6 ) &
           coupler%flow_aux_real_var(2,iconn) = 1.0d-6 
      end if
    end if

    ! assign temperature
    temperature = temperature_at_datum + &
                  temperature_grad(X_DIRECTION)*dist_x + & ! gradient in K/m
                  temperature_grad(Y_DIRECTION)*dist_y + &
                  temperature_grad(Z_DIRECTION)*dist_z 
    coupler%flow_aux_real_var(3,iconn) = temperature

  enddo

  call DestroyOneDimGrid(one_d_grid)

  call DeallocateArray(temperature_array)
  call DeallocateArray(wat_pressure_array)
  call DeallocateArray(wat_density_array)
  call DeallocateArray(oil_pressure_array)
  call DeallocateArray(oil_density_array)
  

end subroutine TOIHydrostaticUpdateCoupler

! ************************************************************************** !


end module HydrostaticMultiPhase_module

