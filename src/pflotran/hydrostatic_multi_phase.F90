module HydrostaticMultiPhase_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Hydrostatic_Common_module

  implicit none

  private

  !module vriables
  !PetscReal :: pb_constant
  PetscBool :: oil_wat_zone
  PetscBool :: oil_gas_zone
  PetscBool :: wat_gas_zone
  PetscBool :: wat_present
  PetscBool :: oil_present
  PetscBool :: gas_present
  PetscBool :: datum_in_wat
  PetscBool :: datum_in_oil
  PetscBool :: datum_in_gas

  !parameters
  PetscReal, parameter :: z_eps = 1.0d-6

#if 0
  !LIQ = water to be consistent with the remainder fo the code
  PetscInt, parameter :: HYDRO_LIQ_PHASE = 1  
  PetscInt, parameter :: HYDRO_GAS_PHASE = 2
  PetscInt, parameter :: HYDRO_OIL_PHASE = 3 

  type :: one_dim_grid_type
    PetscReal :: delta_z
    PetscReal :: min_z
    PetscReal :: max_z
    PetscReal, pointer :: z(:)
    PetscInt :: idatum
  contains 
    procedure :: ElevationIdLoc
  end type one_dim_grid_type
#endif

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

  type(coupler_type), intent(in) :: coupler
  type(option_type) :: option
  type(grid_type), intent(in) :: grid
  type(characteristic_curves_ptr_type), intent(in) :: characteristic_curves_array(:)
  PetscInt, pointer, intent(in) :: sat_func_id(:)
  PetscInt, pointer, intent(in) :: imat(:)

  type(flow_condition_type), pointer :: condition
  class(one_dim_grid_type), pointer :: one_d_grid
  PetscInt :: grid_size
  PetscReal :: max_z, min_z
  PetscReal :: gravity_magnitude
  PetscReal :: datum_tmp_3v(3)
  PetscReal, pointer :: temp_vec(:)
  PetscReal, pointer :: xm_nacl_vec(:)
  PetscReal, pointer :: pb_vec(:)
  PetscReal, pointer :: rs_vec(:)
  PetscReal, pointer :: wat_press_vec(:)
  PetscReal, pointer :: oil_press_vec(:)
  PetscReal, pointer :: gas_press_vec(:)  
  PetscReal, pointer :: wat_den_kg_vec(:)
  PetscReal, pointer :: oil_den_kg_vec(:)
  PetscReal, pointer :: gas_den_kg_vec(:)

  !consider to declare all these input variables as private module variables
  !once the old TOIHydrostaticUpdateCoupler has been removed
  !declared with the same names in HydrostaticPMLoader to make this change easier
  PetscReal :: datum_z
  PetscReal :: press_at_datum
  PetscReal :: owc_z
  PetscReal :: pcow_owc
  PetscReal :: ogc_z
  PetscReal :: pcog_ogc
  PetscReal :: wgc_z
  PetscReal :: pcwg_wgc

  !module variables initialisation
  oil_wat_zone = PETSC_FALSE
  oil_gas_zone = PETSC_FALSE
  wat_gas_zone = PETSC_FALSE
  wat_present = PETSC_FALSE
  oil_present = PETSC_FALSE
  gas_present = PETSC_FALSE  

  condition => coupler%flow_condition

  call HydrostaticPMLoader(condition,option,datum_z,press_at_datum, &
                           owc_z,pcow_owc,ogc_z,pcog_ogc,wgc_z,pcwg_wgc)  
  
  !determine which are the phase and transition zone present
  max_z = max(grid%z_max_global,datum_z)+1.d0 ! add 1m buffer
  min_z = min(grid%z_min_global,datum_z)-1.d0
  
  ! determine what are the transition zones present
  if ( owc_z >= min_z .and. owc_z <= max_z ) then
    oil_wat_zone = PETSC_TRUE
  end if

  if ( ogc_z >= min_z .and. ogc_z <= max_z ) then
    oil_gas_zone = PETSC_TRUE
  end if

  ! wat_gas_zone determined in HydrostaticPMLoader
  if (wat_gas_zone) then
    oil_wat_zone = PETSC_FALSE
    oil_gas_zone = PETSC_FALSE    
  end if  
  
  !determin what are the phases present in the model
  if ( oil_wat_zone .or. wat_gas_zone ) then
    wat_present = PETSC_TRUE
  end if
    
  if ( oil_wat_zone .or. oil_gas_zone ) then
    oil_present = PETSC_TRUE
  end if

  if ( oil_gas_zone .or. wat_gas_zone ) then
    gas_present = PETSC_TRUE
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

  if (dabs(gravity_magnitude-EARTH_GRAVITY) > 0.1d0) then
    option%io_buffer = 'Magnitude of gravity vector is not near 9.81.'
    call printErrMsg(option)
  endif

  !create 1d grid
  !replace datum_tmp_3v with datum_z when replacing TOIHydrostaticUpdateCoupler
  datum_tmp_3v(1:3) = 0.0d0
  datum_tmp_3v(3) = datum_z
  one_d_grid => CreateOneDimGrid(min_z,max_z,datum_tmp_3v)
  
  grid_size = size(one_d_grid%z(:))

  !allocate and initialise one_d_grid data arrays
  ! nullify(temp_vec)
  ! allocate(temp_vec(grid_size))
  ! temp_vec = UNINITIALIZED_DOUBLE
  ! nullify(xm_nacl_vec)
  ! allocate(xm_nacl_vec(grid_size))
  ! xm_nacl_vec = 0.0d0
  ! nullify(pb_vec)
  ! allocate(pb_vec(grid_size))
  ! pb_vec = UNINITIALIZED_DOUBLE
  ! nullify(rs_vec)
  ! allocate(rs_vec(grid_size))
  ! rs_vec = UNINITIALIZED_DOUBLE
  ! nullify(wat_den_kg_vec)
  ! allocate(wat_den_kg_vec(grid_size))
  ! wat_den_kg_vec = UNINITIALIZED_DOUBLE
  ! nullify(oil_den_kg_vec)
  ! allocate(oil_den_kg_vec(grid_size))
  ! oil_den_kg_vec = UNINITIALIZED_DOUBLE
  ! nullify(gas_den_kg_vec)
  ! allocate(gas_den_kg_vec(grid_size))
  ! gas_den_kg_vec = UNINITIALIZED_DOUBLE

  !allocate and initialise one_d_grid data arrays
  call AllocateInitArray(temp_vec,grid_size,UNINITIALIZED_DOUBLE)
  call AllocateInitArray(xm_nacl_vec,grid_size,0.0d0)
  call AllocateInitArray(pb_vec,grid_size,UNINITIALIZED_DOUBLE)
  call AllocateInitArray(rs_vec,grid_size,UNINITIALIZED_DOUBLE)
  call AllocateInitArray(wat_press_vec,grid_size,UNINITIALIZED_DOUBLE)
  call AllocateInitArray(oil_press_vec,grid_size,UNINITIALIZED_DOUBLE)
  call AllocateInitArray(gas_press_vec,grid_size,UNINITIALIZED_DOUBLE)
  call AllocateInitArray(wat_den_kg_vec,grid_size,UNINITIALIZED_DOUBLE)
  call AllocateInitArray(oil_den_kg_vec,grid_size,UNINITIALIZED_DOUBLE)
  call AllocateInitArray(gas_den_kg_vec,grid_size,UNINITIALIZED_DOUBLE)


  ! Interpolate Pb and Temp on 1d grid
  call InterpolatePMInput(condition,option,one_d_grid,press_at_datum, &
                                  temp_vec,xm_nacl_vec,pb_vec,rs_vec)  
  
  ! compute hydrostatic pressure profiles for exisitng phases
  if (datum_in_wat) then
    ! call PhaseHydrostaticPressure(one_d_grid,option%gravity, &
    !           LIQUID_PHASE,pressure_at_datum, &
    !           one_d_grid%idatum,xm_nacl,temperature_array, &
    !           wat_pressure_array,wat_density_array)    
  else if (datum_in_oil) then
    
  else if (datum_in_gas) then
    
  end if
  
  ! compute saturation for exisitng transition zones
  
  ! if oil, water and gas present, check is (sg + sw > 1)
  
  !loop over connections and assign initial solution variables

  ! within connection loop
  ! call HydrostaticPMWriter(condition%flow_aux_real_var(:,iconn), &
  !                            condition%flow_aux_int_var(:,iconn) )

  call DeallocateArray(temp_vec)
  call DeallocateArray(xm_nacl_vec)
  call DeallocateArray(pb_vec)
  call DeallocateArray(rs_vec)
  call DeallocateArray(wat_press_vec)
  call DeallocateArray(oil_press_vec)
  call DeallocateArray(gas_press_vec)
  call DeallocateArray(wat_den_kg_vec)
  call DeallocateArray(oil_den_kg_vec)
  call DeallocateArray(gas_den_kg_vec)


end subroutine HydrostaticMPUpdateCoupler

! ************************************************************************** !

subroutine HydrostaticPMLoader(condition,option,datum_z,press_at_datum, &
                               owc_z,pcow_owc,ogc_z,pcog_ogc,wgc_z,pcwg_wgc)
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
  PetscReal, intent(out) :: datum_z
  PetscReal, intent(out) :: press_at_datum
  PetscReal, intent(out) :: owc_z
  PetscReal, intent(out) :: pcow_owc
  PetscReal, intent(out) :: ogc_z
  PetscReal, intent(out) :: pcog_ogc
  PetscReal, intent(out) :: wgc_z
  PetscReal, intent(out) :: pcwg_wgc
  !add temp_table and pb_table

  PetscBool :: press_at_datum_found

  !assign default values
  press_at_datum_found = PETSC_FALSE
  datum_z = 0.0d0
  press_at_datum = 0.0d0
  owc_z = 0.0d0
  pcow_owc = 0.0d0
  ogc_z = 0.0d0
  pcog_ogc = 0.0d0
  wgc_z = 0.0d0
  pcwg_wgc = 0.0d0


  select case(option%iflowmode)
    case(TOWG_MODE)
      if ( associated(condition%towg%oil_pressure) ) then
        press_at_datum = condition%towg%oil_pressure%dataset%rarray(1)
        press_at_datum_found = PETSC_TRUE
      end if  
      !check if datum defined twice then load value 
      !common to TOIL once defined condition%toil_ims%datum_z
      if ( associated(condition%towg%datum_z) .and. &
           associated(condition%datum) &
          ) then
        option%io_buffer = 'TOWG datum defined twice in datum and datum_z'
        call printErrMsg(option)
      end if
      if ( associated(condition%towg%datum_z) ) then
        datum_z = condition%towg%datum_z%dataset%rarray(1)
      else if (associated(condition%datum)) then
        datum_z = condition%datum%rarray(Z_DIRECTION)
      end if
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
      if ( associated(condition%towg%pcog_ogc) ) then
        pcog_ogc = condition%towg%pcog_ogc%dataset%rarray(1)
      end if

      !if ( .not.associated(condition%towg%pbvz_table) ) then 
        if ( dabs(ogc_z - datum_z) > z_eps ) then
          option%io_buffer = 'Hydrostatic equilibration: PB table not defined, &
                              &datum and OGC must be in the same location'
          call printErrMsg(option)
        end if  
      !end if
      
      if (condition%towg%is_wg_equilibration) then
        ogc_z = 1.0d20
        owc_z = -1.0d20
        wat_gas_zone = PETSC_TRUE
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
  end select

  if (.not.press_at_datum_found) then
    option%io_buffer = 'TOWG Equilibration condition input error: &
                        &a pressure value at datum must be input'
    call printErrMsg(option)
  end if

end subroutine HydrostaticPMLoader

! ************************************************************************** !

subroutine InterpolatePMInput(condition,option,one_d_grid,press_at_datum, &
                             temp_vec,xm_nacl_vec,pb_vec,rs_vec)
  ! 
  ! Interpolate thermodynamics properties into 1d grid
  !
  ! Author: Paolo Orsini
  ! Date: 01/09/19  

  use Option_module
  use Condition_module

  implicit none

  type(flow_condition_type), intent(in) :: condition
  type(option_type) :: option
  class(one_dim_grid_type), intent(in) :: one_d_grid
  PetscReal, intent(in) :: press_at_datum
  PetscReal, intent(out) :: temp_vec(:)
  PetscReal, intent(out) :: xm_nacl_vec(:)
  PetscReal, intent(out) :: pb_vec(:)
  PetscReal, intent(out) :: rs_vec(:)

  PetscReal :: xm_nacl, temperature
  PetscReal :: temp_grad(3)
  PetscInt :: grid_size
  PetscBool :: temp_grad_given

  xm_nacl = 0.0d0
  temp_grad = UNINITIALIZED_DOUBLE
  temperature = UNINITIALIZED_DOUBLE
  grid_size = size(one_d_grid%z(:))
  temp_grad_given = PETSC_FALSE

  select case(option%iflowmode)
    case(TOWG_MODE)
      !if pbvz_table then
      ! do lookup for each z and copy values onto pb_vec
      !else
      pb_vec(1:grid_size) = press_at_datum
      !end if
      if ( associated(condition%towg%temperature) ) then
        !gradient or constant 
        temperature = condition%towg%temperature%dataset%rarray(1)
        if ( associated(condition%towg%temperature%gradient ) ) then
          temp_grad(1:3) = condition%towg%temperature%gradient%rarray(1:3)
          temp_grad_given = PETSC_TRUE
        end if  
      end if
    case(TOIL_IMS_MODE)
      if ( associated(condition%toil_ims%temperature) ) then
        !gradient or constant 
        temperature = condition%toil_ims%temperature%dataset%rarray(1)
        if ( associated(condition%toil_ims%temperature%gradient ) ) then
          temp_grad(1:3) = condition%toil_ims%temperature%gradient%rarray(1:3)
          temp_grad_given = PETSC_TRUE
        end if  
      end if
  end select
  
  !check also if both temperature and rtempvz_table have been defined and return 
  !an error. This needs to be implmenented in cocordination with condition read
  !must allow to read read either temperature (with/wihtout gread) or rtempvz_table
  ! if not temp, (and temp_grad) or rtempvz_table return an error
  
  !if (associated(condition%rtempvz_table)) then
  !else if
    if (temp_grad_given) then
      call CompVertTempProfile(one_d_grid,temp_grad,temperature, &
                               temp_vec)
    else
      temp_vec(1:grid_size)=temperature
    end if  
  !end if 


  !if (associated(condition%saltvz_table)) then
  ! do lookup for each z and copy values onto xm_nacl_vec
  !else
    xm_nacl = option%m_nacl * FMWNACL
    xm_nacl = xm_nacl /(1.d3 + xm_nacl)
    xm_nacl_vec(1:grid_size) = xm_nacl
  !

end subroutine InterpolatePMInput

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

#if 0
! ************************************************************************** !

function GetCouplerCellOnPhaseConact(coupler,grid,imat,z_phase_contact,option)

  ! Returns the ghosted_id of the closest cell to phase contact  
  ! Where more cells have the same minimum distance from z_phase_contact,
  ! the first cell in the list will be selected as for the minloc function 
  !
  ! Author: Paolo Orsini
  ! Date: 12/31/15
  ! 

  use Grid_module
  use Coupler_module
  use Utility_module
  use Option_module

  implicit none

  type(coupler_type), intent(in) :: coupler
  type(grid_type), intent(in) :: grid
  PetscInt, pointer, intent(in) :: imat(:)
  PetscReal, intent(in) :: z_phase_contact
  type(option_type) :: option

  PetscInt :: GetCouplerCellOnPhaseConact

  PetscInt :: iconn, local_id, ghosted_id
  PetscInt :: iconn_min(1)
  PetscReal, pointer :: phase_contact_dist(:) 

  allocate(phase_contact_dist(coupler%connection_set%num_connections))
  phase_contact_dist = 1.d20
  
  do iconn=1,coupler%connection_set%num_connections
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    if (imat(ghosted_id) <= 0) cycle
    phase_contact_dist(iconn) = dabs(grid%z(ghosted_id) - z_phase_contact)
  end do

  iconn_min = minloc(phase_contact_dist(:))

  local_id = coupler%connection_set%id_dn(iconn_min(1))

  !write(*,*) "in get ghost -rank=", option%myrank, "icon_min=", iconn_min(1), &
  !           " ghost =", grid%nL2G(local_id), "size dist = ", size(phase_contact_dist(:))

  GetCouplerCellOnPhaseConact = grid%nL2G(local_id) 

  call DeallocateArray(phase_contact_dist)
 

end function GetCouplerCellOnPhaseConact

! ************************************************************************** !

subroutine PhaseHydrostaticPressure(one_d_grid,gravity,iphase,press_start, &
                                    id_start,xm_nacl,temp,press,den_kg)
  ! 
  ! Compute an hydrostatic pressure profile for a given phase and 1D 
  ! 1D discretisation
  ! The "starting point" for the computation of the pressure profile
  ! is the elevation where a reference pressure is given, it can be the datum,
  ! or one of the phase contact interface (e.g. OWC, OGC, etc)
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 

  implicit none

  class(one_dim_grid_type) :: one_d_grid
  PetscReal, intent(in) :: gravity(:) ! this is option%gravity
  PetscInt, intent(in) :: iphase
  PetscReal, intent(in) :: press_start
  PetscInt, intent(in) :: id_start
  PetscReal, intent(in) :: temp(:)
  PetscReal, intent(in) :: xm_nacl
  PetscReal, intent(out) :: press(:)
  PetscReal, intent(out) :: den_kg(:)

  PetscReal :: pressure, pressure0, rho, rho_one, rho_kg, rho_zero
  PetscInt :: ipressure, num_iteration
  
  if(iphase == HYDRO_GAS_PHASE ) then
    print *, "PhaseHydrostaticPressure does not support gas"
    stop
  end if

  rho_kg = PhaseDensity(iphase,press_start,temp(id_start),xm_nacl)

  ! fill properties for reference pressure
  den_kg(id_start) = rho_kg
  press(id_start) = press_start 

  pressure0 = press_start 
  rho_zero = rho_kg 
  do ipressure=id_start+1,size(one_d_grid%z(:))
    rho_kg = PhaseDensity(iphase,pressure0,temp(ipressure),xm_nacl)
    num_iteration = 0
    do 
      pressure = pressure0 + 0.5d0*(rho_kg+rho_zero) * &
                 gravity(Z_DIRECTION) * one_d_grid%delta_z
      rho_one = PhaseDensity(iphase,pressure,temp(ipressure),xm_nacl)
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
    press(ipressure) = pressure
    den_kg(ipressure) = rho_kg
    pressure0 = pressure
  enddo

  ! compute pressures below one_d_grid%z(id_start), if any
  pressure0 = press(id_start)
  rho_zero = den_kg(id_start)
  do ipressure=id_start-1,1,-1
    rho_kg = PhaseDensity(iphase,pressure0,temp(ipressure),xm_nacl)
    num_iteration = 0
    do                   ! notice the negative sign (-) here
      pressure = pressure0 - 0.5d0*(rho_kg+rho_zero) * &
                 gravity(Z_DIRECTION) * one_d_grid%delta_z
      rho_one = PhaseDensity(iphase,pressure,temp(ipressure),xm_nacl)
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
    press(ipressure) = pressure
    den_kg(ipressure) = rho_kg
    pressure0 = pressure
  enddo

end subroutine PhaseHydrostaticPressure

! ************************************************************************** !
function PhaseDensity(iphase,p,t,xm_nacl) 
  !
  ! computes phase density given the specified phase
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 

  use EOS_Oil_module
  use EOS_Water_module
  use EOS_Gas_module ! when gas is considered

  implicit none

  PetscInt, intent(in) :: iphase
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: xm_nacl
 
  PetscReal :: PhaseDensity ! kg/m3 
  PetscReal :: dw_mol
  PetscReal :: aux(1)

  PetscErrorCode :: ierr

  select case(iphase)
    case(HYDRO_LIQ_PHASE)
      aux(1) = xm_nacl
      call EOSWaterDensityExt(t,p,aux,PhaseDensity,dw_mol,ierr)
    case(HYDRO_GAS_PHASE)
      !call EOSGasDensityNoDerive(t,p,PhaseDensity,ierr)
      ! rho_kg = rho * GAS_FMW (to get gas FMW currenlty mode specific)
      ! gas_fmw should be defined in gas_eos 
    case(HYDRO_OIL_PHASE)
      call EOSOilDensity(t,p,PhaseDensity,ierr)
      PhaseDensity = PhaseDensity * EOSOilGetFMW() 
  end select

end function PhaseDensity 

! ************************************************************************** !

function PressInterp(ipressure,dist_x,dist_y,dist_z_for_pressure,gravity, &
                     pressure_array,density_array,pressure_gradient)

  ! Computes hydrostatic pressure over dz from a reference pressure and
  ! a density, combining possible horizontal gradients  
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15


  implicit none

  PetscInt, intent(in) :: ipressure
  PetscReal, intent(in) :: dist_x
  PetscReal, intent(in) :: dist_y
  PetscReal, intent(in) :: dist_z_for_pressure
  PetscReal, intent(in) :: gravity(:)
  PetscReal, intent(in) :: pressure_array(:)
  PetscReal, intent(in) :: density_array(:)
  PetscReal, intent(in) :: pressure_gradient(:)

  PetscReal :: PressInterp

  PressInterp = pressure_array(ipressure) + &
                density_array(ipressure) * gravity(Z_DIRECTION) * &
                dist_z_for_pressure + &
                pressure_gradient(X_DIRECTION) * dist_x + & ! gradient in Pa/m
                pressure_gradient(Y_DIRECTION) * dist_y

end function PressInterp

! ************************************************************************** !

!function HydrostaticPressOverDz(dz,gravity,press_ref,density_ref)
!
!  implicit none
!
!  PetscReal, intent(in) :: dz
!  PetscReal, intent(in) :: gravity(:)
!  PetscReal, intent(in) :: press_ref
!  PetscReal, intent(in) :: density_ref
!
!  PetscReal :: HydrostaticPressOverDz
!
!  HydrostaticPressOverDz = press_ref + density_ref * gravity(Z_DIRECTION) * dz
!
!end function HydrostaticPressOverDz

! ************************************************************************** !

function PressGrad(dist_x,dist_y,dist_z,press_ref,pressure_gradient)

  ! Computes pressure from a reference pressure, distances in the three 
  ! directions, and a 3d local pressure gradient  
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15


  implicit none

  PetscReal, intent(in) :: dist_x
  PetscReal, intent(in) :: dist_y
  PetscReal, intent(in) :: dist_z
  PetscReal, intent(in) :: press_ref
  PetscReal, intent(in) :: pressure_gradient(:)

  PetscReal :: PressGrad

  PressGrad = press_ref + &
              pressure_gradient(X_DIRECTION)*dist_x + & ! gradient in Pa/m
              pressure_gradient(Y_DIRECTION)*dist_y + &
              pressure_gradient(Z_DIRECTION)*dist_z 

end function PressGrad

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
function ElevationIdLoc(this,elevation)
  ! 
  ! Detect id location in one_d_grid given an absolute elevation
  ! Note: absolute elevation, not relative to the datum 
  !
  ! Author: Paolo Orsini
  ! Date: 12/21/15
  ! 

  implicit none

  class(one_dim_grid_type) :: this
  PetscReal, intent(in) :: elevation !this is the global elevation
  
  PetscInt :: ElevationIdLoc
   
  ElevationIdLoc = int( (elevation-this%min_z)/(this%max_z-this%min_z) * &
                        dble(size(this%z(:))) ) + 1 

  !  idatum = int((datum(Z_DIRECTION)-min_z)/(max_z-min_z) * &
  !               dble(num_pressures))+1

end function ElevationIdLoc

! ************************************************************************** !

subroutine DestroyOneDimGrid(one_d_grid)
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/30/15
  ! 
  ! destroys OneDimGrid

  use Utility_module 

  implicit none

  class(one_dim_grid_type), pointer :: one_d_grid

  if (.not.associated(one_d_grid) ) return

  call DeallocateArray(one_d_grid%z)

  deallocate(one_d_grid)
  nullify(one_d_grid)

end subroutine DestroyOneDimGrid

! ************************************************************************** !
#endif

end module HydrostaticMultiPhase_module

