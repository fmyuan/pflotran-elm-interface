module Test_Characteristic_Curves_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use pFUnit_mod
  use Characteristic_Curves_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_Common_module
  use Characteristic_Curves_VG_module
  use Characteristic_Curves_WIPP_module
  use Option_module

  implicit none


  public :: Test_Characteristic_Curves

  interface Test_Characteristic_Curves
     module procedure newTest_Characteristic_Curves
  end interface Test_Characteristic_Curves

! ************************************************************************** !
!  @TestCase
  type, extends(TestCase) :: Test_Characteristic_Curves
      type(option_type), pointer :: option
      class(characteristic_curves_type), pointer :: cc_const
      class(characteristic_curves_type), pointer :: cc_bcb
      class(characteristic_curves_type), pointer :: cc_bcm
      class(characteristic_curves_type), pointer :: cc_vgm
      class(characteristic_curves_type), pointer :: cc_vgb
      class(characteristic_curves_type), pointer :: cc_vgt
      class(characteristic_curves_type), pointer :: cc_lb
      class(characteristic_curves_type), pointer :: cc_lm
      class(characteristic_curves_type), pointer :: cc_mk3
      class(characteristic_curves_type), pointer :: cc_mk4
      class(characteristic_curves_type), pointer :: cc_krp1
      class(characteristic_curves_type), pointer :: cc_krp2
      class(characteristic_curves_type), pointer :: cc_krp3
      class(characteristic_curves_type), pointer :: cc_krp4
      class(characteristic_curves_type), pointer :: cc_krp5
      class(characteristic_curves_type), pointer :: cc_krp8
      class(characteristic_curves_type), pointer :: cc_krp9
      class(characteristic_curves_type), pointer :: cc_krp11
      class(characteristic_curves_type), pointer :: cc_krp12
      procedure(runMethod), pointer :: userMethod => null()
    contains
      procedure :: setUp     
      procedure :: tearDown
      procedure :: runMethod
  end type Test_Characteristic_Curves

contains

! ************************************************************************** !

  function newTest_Characteristic_Curves(name, userMethod) result(test)

    implicit none

    character(len=*), intent(in) :: name
    procedure(runMethod) :: userMethod

    type(Test_Characteristic_Curves) :: test

    call test%setName(name)
    test%userMethod => userMethod

  end function newTest_Characteristic_Curves

! ************************************************************************** !

  subroutine setUp(this)

    use PFLOTRAN_Constants_module, only : MAXSTRINGLENGTH

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this
    character(len=MAXSTRINGLENGTH) :: error_string
    PetscInt :: error

    this%option => OptionCreate()

  ! Setting up the characteristic curve for Brooks Corey Burdine
    this%cc_bcb => CharacteristicCurvesCreate()
    this%cc_bcb%saturation_function => SFBCCreate()
    this%cc_bcb%saturation_function%Sr = 0.2d0
    this%cc_bcb%saturation_function%pcmax = 0.999d8
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_BC_type)
        sf%lambda = 0.7d0
        sf%alpha = 9.869d-6
      class default
        print *, 'not BC type in setUp'
    end select
    this%cc_bcb%liq_rel_perm_function => RPFBurdineBCLiqCreate()
    this%cc_bcb%liq_rel_perm_function%Sr = 0.2d0
    select type(sf=>this%cc_bcb%liq_rel_perm_function)
      class is(rpf_Burdine_BC_liq_type)
        sf%lambda = 0.7d0
      class default
        print *, 'not Burdine BC Liq type'
    end select
    this%cc_bcb%gas_rel_perm_function => RPFBurdineBCGasCreate()
    this%cc_bcb%gas_rel_perm_function%Sr = 0.2d0
    select type(sf=>this%cc_bcb%gas_rel_perm_function)
      class is(rpf_Burdine_BC_gas_type)
        sf%lambda = 0.7d0
        sf%Srg = 1.d-5
      class default
        print *, 'not Burdine BC Gas type'
    end select
    error_string = 'pFUnit Brooks Corey saturation function'
    call this%cc_bcb%saturation_function%SetupPolynomials(this%option, &
                                                          error_string)

  ! Setting up the characteristic curve for Brooks Corey Mualem
    this%cc_bcm => CharacteristicCurvesCreate()
    ! do not set up saturation function as it is tested with cc_bcb
    this%cc_bcm%liq_rel_perm_function => RPFMualemBCLiqCreate()
    this%cc_bcm%liq_rel_perm_function%Sr = 0.2d0
    select type(sf=>this%cc_bcm%liq_rel_perm_function)
      class is(rpf_Mualem_BC_liq_type)
        sf%lambda = 0.7d0
      class default
        print *, 'not Mualem BC Liq type'
    end select
    this%cc_bcm%gas_rel_perm_function => RPFMualemBCGasCreate()
    this%cc_bcm%gas_rel_perm_function%Sr = 0.2d0
    select type(sf=>this%cc_bcm%gas_rel_perm_function)
      class is(rpf_Mualem_BC_gas_type)
        sf%lambda = 0.7d0
        sf%Srg = 1.d-5
      class default
        print *, 'not Mualem BC Gas type'
    end select

  ! Setting up the characteristic curve for van Genuchten Mualem
    this%cc_vgm => CharacteristicCurvesCreate()
    this%cc_vgm%saturation_function => SFVGCreate()
    this%cc_vgm%saturation_function%Sr = 0.143d0
    this%cc_vgm%saturation_function%pcmax = 0.999d8
    select type(sf=>this%cc_vgm%saturation_function)
      class is(sat_func_VG_type)
        error = sf%set_m(0.527d0)
        error = sf%set_alpha(5.1054d-5)
      class default
        print *, 'not VG type in setUp'
    end select
    this%cc_vgm%liq_rel_perm_function => RPFMualemVGLiqCreate()
    this%cc_vgm%liq_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgm%liq_rel_perm_function)
      class is(rpf_Mualem_VG_liq_type)
        error = sf%set_m(0.527d0)
      class default
        print *, 'not Mualem VG Liq type'
    end select
    this%cc_vgm%gas_rel_perm_function => RPFMualemVGGasCreate()
    this%cc_vgm%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgm%gas_rel_perm_function)
      class is(rpf_Mualem_VG_gas_type)
        error = sf%set_m(0.527d0)
        sf%Srg = 0.d0
      class default
        print *, 'not Mualem VG Gas type'
    end select

  ! Setting up the characteristic curve for van Genuchten Burdine
    this%cc_vgb => CharacteristicCurvesCreate()
    ! do not set up saturation function as it is tested with cc_vgm
    this%cc_vgb%liq_rel_perm_function => RPFBurdineVGLiqCreate()
    this%cc_vgb%liq_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgb%liq_rel_perm_function)
      class is(rpf_Burdine_VG_liq_type)
        error = sf%set_m(0.527d0)
      class default
        print *, 'not Burdine VG Liq type'
    end select
    this%cc_vgb%gas_rel_perm_function => RPFBurdineVGGasCreate()
    this%cc_vgb%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgb%gas_rel_perm_function)
      class is(rpf_Burdine_VG_gas_type)
        error = sf%set_m(0.527d0)
        sf%Srg = 0.01d0
      class default
        print *, 'not Burdine VG Gas type'
    end select

  ! Setting up the characteristic curve for van Genuchten TOUGH2 IRP7
    this%cc_vgt => CharacteristicCurvesCreate()
    ! do not set up saturation function as it is tested with cc_vgm
    this%cc_vgt%gas_rel_perm_function => RPFTOUGH2IRP7GasCreate()
    this%cc_vgt%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgt%gas_rel_perm_function)
      class is(rpf_TOUGH2_IRP7_gas_type)
        sf%Srg = 0.01d0
      class default
        print *, 'not TOUGH2 IRP7 Gas type'
    end select

  ! Setting up the characteristic curve for Linear Mualem
    this%cc_lm => CharacteristicCurvesCreate()
    this%cc_lm%saturation_function => SFLinearCreate()
    this%cc_lm%saturation_function%Sr = 0.143d0
    this%cc_lm%saturation_function%pcmax = 0.999d8
    select type(sf=>this%cc_lm%saturation_function)
      class is(sat_func_Linear_type)
        sf%alpha = 5.1054d-5
      class default
        print *, 'not Linear type in setUp'
    end select    
    this%cc_lm%liq_rel_perm_function => RPFMualemLinearLiqCreate()
    this%cc_lm%liq_rel_perm_function%Sr = 0.143d0
    this%cc_lm%gas_rel_perm_function => RPFMualemLinearGasCreate()
    this%cc_lm%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_lm%gas_rel_perm_function)
      class is(rpf_Mualem_Linear_gas_type)
        sf%Srg = 0.01d0
      class default
        print *, 'not Mualem Linear Gas type'
    end select

  ! Setting up the characteristic curve for Linear Burdine
    this%cc_lb => CharacteristicCurvesCreate()
    ! do not set up saturation function as it is tested with cc_lm
    this%cc_lb%liq_rel_perm_function => RPFBurdineLinearLiqCreate()
    this%cc_lb%liq_rel_perm_function%Sr = 0.143d0
    this%cc_lb%gas_rel_perm_function => RPFBurdineLinearGasCreate()
    this%cc_lb%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_lb%gas_rel_perm_function)
      class is(rpf_Burdine_Linear_gas_type)
        sf%Srg = 0.01d0
      class default
        print *, 'not Burdine Linear Gas type'
    end select

  ! Setting up the characteristic curve for 3-param modified Kosugi
    this%cc_mk3 => CharacteristicCurvesCreate()
    this%cc_mk3%saturation_function => SFmKCreate()
    this%cc_mk3%saturation_function%Sr = 1.53D-1
    this%cc_mk3%saturation_function%pcmax = 0.999d8
    select type(sf=>this%cc_mk3%saturation_function)
      class is(sat_func_mK_type)
        sf%nparam = 3
        sf%sigmaz = 3.37D-1
        sf%muz = -6.3D0
        sf%rmax = 2.52D-3
      class default
        print *, 'no mK3 saturation function type in setUp'
    end select
    this%cc_mk3%liq_rel_perm_function => RPFmKLiqCreate()
    this%cc_mk3%liq_rel_perm_function%Sr = 1.53D-1
    select type(sf=>this%cc_mk3%liq_rel_perm_function)
      class is(rpf_mK_liq_type)
        sf%sigmaz = 3.37D-1
      class default
        print *, 'no mK3 Liq rel perm type in setUp'
    end select
    this%cc_mk3%gas_rel_perm_function => RPFmKGasCreate()
    this%cc_mk3%gas_rel_perm_function%Sr = 1.53D-1
    select type(sf=>this%cc_mk3%gas_rel_perm_function)
      class is(rpf_mK_gas_type)
        sf%sigmaz = 3.37D-1
        sf%Srg = 1.0D-3
      class default
        print *, 'no mK3 Gas rel perm type in setUp'
    end select

  ! Setting up the characteristic curve for 4-param modified Kosugi
    this%cc_mk4 => CharacteristicCurvesCreate()
    this%cc_mk4%saturation_function => SFmKCreate()
    this%cc_mk4%saturation_function%Sr = 1.53D-1
    this%cc_mk4%saturation_function%pcmax = 0.999d8
    select type(sf=>this%cc_mk4%saturation_function)
      class is(sat_func_mK_type)
        sf%nparam = 4
        sf%sigmaz = 3.37D-1
        sf%muz = -6.3D0
        sf%rmax = 2.52D-3
        sf%r0 = 1.07D-4
      class default
        print *, 'no mK4 saturation function type in setUp'
    end select
    this%cc_mk4%liq_rel_perm_function => RPFmKLiqCreate()
    this%cc_mk4%liq_rel_perm_function%Sr = 1.53D-1
    select type(sf=>this%cc_mk4%liq_rel_perm_function)
      class is(rpf_mK_liq_type)
        sf%sigmaz = 3.37D-1
      class default
        print *, 'no mK4 Liq rel perm type in setUp'
    end select
    this%cc_mk4%gas_rel_perm_function => RPFmKGasCreate()
    this%cc_mk4%gas_rel_perm_function%Sr = 1.53D-1
    select type(sf=>this%cc_mk4%gas_rel_perm_function)
      class is(rpf_mK_gas_type)
        sf%sigmaz = 3.37D-1
        sf%Srg = 1.0D-3
      class default
        print *, 'no mK4 Gas rel perm type in setUp'
    end select

  ! Setting up the characteristic curve constant capillary pressure 
  ! saturation function
    this%cc_const => CharacteristicCurvesCreate()
    this%cc_const%saturation_function => SFConstantCreate()
    select type(sf=>this%cc_const%saturation_function)
      class is(sat_func_constant_type)
        sf%constant_capillary_pressure = 1.d5
        sf%constant_saturation = 0.5d0
      class default
        print *, 'no constant saturation function type in setUp'
    end select
    
  ! Setting up the characteristic curve for BRAGFLO_KRP1
    this%cc_krp1 => CharacteristicCurvesCreate()
    this%cc_krp1%saturation_function => SFKRP1Create()
    this%cc_krp1%saturation_function%Sr = 0.20d0
    this%cc_krp1%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp1%saturation_function)
      class is(sat_func_KRP1_type)
        sf%Srg = 0.25d0
        sf%m = 0.60d0
        sf%alpha = 5.0d-5
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 1
    end select
    this%cc_krp1%liq_rel_perm_function => RPFKRP1LiqCreate()
    this%cc_krp1%liq_rel_perm_function%Sr = 0.20d0
    select type(sf=>this%cc_krp1%liq_rel_perm_function)
      class is(rpf_KRP1_liq_type)
        sf%m = 0.60d0
        sf%Srg = 0.25d0
    end select
    this%cc_krp1%gas_rel_perm_function => RPFKRP1GasCreate()
    this%cc_krp1%gas_rel_perm_function%Sr = 0.25d0
    select type(sf=>this%cc_krp1%gas_rel_perm_function)
      class is(rpf_KRP1_gas_type)
        sf%m = 0.60d0
        sf%Srg = 0.25d0
    end select
    
  ! Setting up the characteristic curve for BRAGFLO_KRP2
    this%cc_krp2 => CharacteristicCurvesCreate()
    this%cc_krp2%saturation_function => SFKRP2Create()
    this%cc_krp2%saturation_function%Sr = 0.10d0
    this%cc_krp2%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp2%saturation_function)
      class is(sat_func_KRP2_type)
        sf%lambda = 0.25d0
        sf%alpha = 5.0d-5
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp2%liq_rel_perm_function => RPFKRP2LiqCreate()
    this%cc_krp2%liq_rel_perm_function%Sr = 0.10d0
    select type(sf=>this%cc_krp2%liq_rel_perm_function)
      class is(rpf_KRP2_liq_type)
        sf%lambda = 0.25d0
    end select
    this%cc_krp2%gas_rel_perm_function => RPFKRP2GasCreate()
    this%cc_krp2%gas_rel_perm_function%Sr = 0.10d0
    select type(sf=>this%cc_krp2%gas_rel_perm_function)
      class is(rpf_KRP2_gas_type)
        sf%lambda = 0.25d0
    end select
    
  ! Setting up the characteristic curve for BRAGFLO_KRP3
    this%cc_krp3 => CharacteristicCurvesCreate()
    this%cc_krp3%saturation_function => SFKRP3Create()
    this%cc_krp3%saturation_function%Sr = 0.15d0
    this%cc_krp3%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp3%saturation_function)
      class is(sat_func_KRP3_type)
        sf%Srg = 0.25d0
        sf%lambda = 0.35d0
        sf%alpha = 5.0d-3
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp3%liq_rel_perm_function => RPFKRP3LiqCreate()
    this%cc_krp3%liq_rel_perm_function%Sr = 0.15d0
    select type(sf=>this%cc_krp3%liq_rel_perm_function)
      class is(rpf_KRP3_liq_type)
        sf%lambda = 0.35d0
        sf%Srg = 0.25d0
    end select
    this%cc_krp3%gas_rel_perm_function => RPFKRP3GasCreate()
    this%cc_krp3%gas_rel_perm_function%Sr = 0.15d0
    select type(sf=>this%cc_krp3%gas_rel_perm_function)
      class is(rpf_KRP3_gas_type)
        sf%lambda = 0.35d0
        sf%Srg = 0.25d0
    end select
    
  ! Setting up the characteristic curve for BRAGFLO_KRP4
    this%cc_krp4 => CharacteristicCurvesCreate()
    this%cc_krp4%saturation_function => SFKRP4Create()
    this%cc_krp4%saturation_function%Sr = 0.15d0
    this%cc_krp4%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp4%saturation_function)
      class is(sat_func_KRP4_type)
        sf%Srg = 0.25d0
        sf%lambda = 0.35d0
        sf%alpha = 5.0d-3
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp4%liq_rel_perm_function => RPFKRP4LiqCreate()
    this%cc_krp4%liq_rel_perm_function%Sr = 0.15d0
    select type(sf=>this%cc_krp4%liq_rel_perm_function)
      class is(rpf_KRP4_liq_type)
        sf%lambda = 0.35d0
        sf%Srg = 0.25d0
    end select
    this%cc_krp4%gas_rel_perm_function => RPFKRP4GasCreate()
    this%cc_krp4%gas_rel_perm_function%Sr = 0.15d0
    select type(sf=>this%cc_krp4%gas_rel_perm_function)
      class is(rpf_KRP4_gas_type)
        sf%lambda = 0.35d0
        sf%Srg = 0.25d0
    end select
    
  ! Setting up the characteristic curve for BRAGFLO_KRP5
    this%cc_krp5 => CharacteristicCurvesCreate()
    this%cc_krp5%saturation_function => SFKRP5Create()
    this%cc_krp5%saturation_function%Sr = 0.05d0
    this%cc_krp5%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp5%saturation_function)
      class is(sat_func_KRP5_type)
        sf%alpha = 5.0d-4
        sf%Srg = 0.20d0
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp5%liq_rel_perm_function => RPFKRP5LiqCreate()
    this%cc_krp5%liq_rel_perm_function%Sr = 0.05d0
    select type(sf=>this%cc_krp5%liq_rel_perm_function)
      class is(rpf_KRP5_liq_type)
        sf%Srg = 0.20d0
    end select
    this%cc_krp5%gas_rel_perm_function => RPFKRP5GasCreate()
    this%cc_krp5%gas_rel_perm_function%Sr = 0.05d0
    select type(sf=>this%cc_krp5%gas_rel_perm_function)
      class is(rpf_KRP5_gas_type)
        sf%Srg = 0.20d0
    end select

  ! Setting up the characteristic curve for BRAGFLO_KRP8
    this%cc_krp8 => CharacteristicCurvesCreate()
    this%cc_krp8%saturation_function => SFKRP8Create()
    this%cc_krp8%saturation_function%Sr = 0.10d0
    this%cc_krp8%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp8%saturation_function)
      class is(sat_func_KRP8_type)
        sf%Srg = 0.25d0
        sf%m = 0.45d0
        sf%alpha = 5.0d-5
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp8%liq_rel_perm_function => RPFKRP8LiqCreate()
    this%cc_krp8%liq_rel_perm_function%Sr = 0.10d0
    select type(sf=>this%cc_krp8%liq_rel_perm_function)
      class is(rpf_KRP8_liq_type)
        sf%m = 0.45d0
    end select
    this%cc_krp8%gas_rel_perm_function => RPFKRP8GasCreate()
    this%cc_krp8%gas_rel_perm_function%Sr = 0.10d0
    select type(sf=>this%cc_krp8%gas_rel_perm_function)
      class is(rpf_KRP8_gas_type)
        sf%m = 0.45d0
    end select
    
  ! Setting up the characteristic curve for BRAGFLO_KRP9
    this%cc_krp9 => CharacteristicCurvesCreate()
    this%cc_krp9%saturation_function => SFKRP9Create()
    this%cc_krp9%saturation_function%Sr = 0.05d0
    this%cc_krp9%saturation_function%pcmax = 0.999d7
    this%cc_krp9%liq_rel_perm_function => RPFKRP9LiqCreate()
    this%cc_krp9%liq_rel_perm_function%Sr = 0.05d0
    this%cc_krp9%gas_rel_perm_function => RPFKRP9GasCreate()
    this%cc_krp9%gas_rel_perm_function%Sr = 0.05d0
    
  ! Setting up the characteristic curve for BRAGFLO_KRP11
    this%cc_krp11 => CharacteristicCurvesCreate()
    this%cc_krp11%saturation_function => SFKRP11Create()
    this%cc_krp11%saturation_function%Sr = 0.08d0
    this%cc_krp11%saturation_function%pcmax = 0.999d7
    this%cc_krp11%liq_rel_perm_function => RPFKRP11LiqCreate()
    this%cc_krp11%liq_rel_perm_function%Sr = 0.08d0
    select type(sf=>this%cc_krp11%liq_rel_perm_function)
      class is(rpf_KRP11_liq_type)
        sf%tolc = 0.10
        sf%Srg = 0.18d0
    end select
    this%cc_krp11%gas_rel_perm_function => RPFKRP11GasCreate()
    this%cc_krp11%gas_rel_perm_function%Sr = 0.08d0
    select type(sf=>this%cc_krp11%gas_rel_perm_function)
      class is(rpf_KRP11_gas_type)
        sf%tolc = 0.10
        sf%Srg = 0.18d0
    end select
    
  ! Setting up the characteristic curve for BRAGFLO_KRP12
    this%cc_krp12 => CharacteristicCurvesCreate()
    this%cc_krp12%saturation_function => SFKRP12Create()
    this%cc_krp12%saturation_function%Sr = 0.12d0
    this%cc_krp12%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp12%saturation_function)
      class is(sat_func_KRP12_type)
        sf%lambda = 0.2d0
        sf%alpha = 5.0d-4
        sf%s_min = 0.015d0
        sf%s_effmin = 0.005d0
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp12%liq_rel_perm_function => RPFKRP12LiqCreate()
    this%cc_krp12%liq_rel_perm_function%Sr = 0.12d0
    select type(sf=>this%cc_krp12%liq_rel_perm_function)
      class is(rpf_KRP12_liq_type)
        sf%lambda = 0.20d0
        sf%Srg = 0.20d0
    end select
    this%cc_krp12%gas_rel_perm_function => RPFKRP12GasCreate()
    this%cc_krp12%gas_rel_perm_function%Sr = 0.12d0
    select type(sf=>this%cc_krp12%gas_rel_perm_function)
      class is(rpf_KRP12_gas_type)
        sf%lambda = 0.20d0
        sf%Srg = 0.20d0
    end select

  end subroutine setUp

! ************************************************************************** !

  subroutine tearDown(this)

    implicit none
    class (Test_Characteristic_Curves), intent(inout) :: this

    call OptionDestroy(this%option)
    call CharacteristicCurvesDestroy(this%cc_bcb)

  end subroutine tearDown

! ************************************************************************** !

  subroutine runMethod(this)
    implicit none
    class (Test_Characteristic_Curves), intent(inout) :: this
    call this%userMethod()
  end subroutine runMethod

! ************************************************************************** !

!  @Test
  subroutine testSF_BC_SetupPolynomials(this)

    implicit none

    class(Test_characteristic_Curves), intent(inout) :: this
    
    PetscReal :: values(4)
    PetscReal, parameter :: tolerance = 1.d-8
    PetscInt :: i
    character(len=128) :: string

    ! pressure polynomial
    values = [-4.6122570934041036d0, 1.4087313882738163d-4, &
              -1.1098865178492886d-9, 2.6190024815166203d-15]
    do i = 1, 4
      write(string,*) i
      string = 'Brooks-Corey-Burdine pressure polynomial coefficient #' // &
               trim(adjustl(string))
#line 576 "test_characteristic_curves.pf"
  call assertEqual(values(i), this%cc_bcb%saturation_function%pres_poly%coefficients(i), dabs(values(i))*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 576) )
  if (anyExceptions()) return
#line 577 "test_characteristic_curves.pf"
    enddo

    ! saturation spline
    values = [-83508464.603000879d0, 173197055.36650354d0, &
              -89688590.763502657d0, 0.d0]
    do i = 1, 3
      write(string,*) i
      string = 'Brooks-Corey-Burdine saturation spline coefficient #' // &
               trim(adjustl(string))
#line 586 "test_characteristic_curves.pf"
  call assertEqual(values(i), this%cc_bcb%saturation_function%sat_poly%coefficients(i),  dabs(values(i))*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 586) )
  if (anyExceptions()) return
#line 587 "test_characteristic_curves.pf"
    enddo

  end subroutine testSF_BC_SetupPolynomials

! ************************************************************************** !

!  @Test
  subroutine testsf_Brooks_Corey(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: dsat_pres
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: dkr_p
    PetscReal :: relative_permeability
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! saturation = f(capillary_pressure) below the polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_BC_type)
        capillary_pressure = 0.94d0/sf%alpha
      class default
        print *, 'not bc type in testsf_Brooks_Corey'
    end select
    call this%cc_bcb%saturation_function% &
         Saturation(capillary_pressure, &
                    liquid_saturation, dsat_pres,this%option)
    call this%cc_bcb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'Brooks-Corey-Burdine saturation as a function of capillary &
             &pressure below polynomial fit'
    value = 1.d0
#line 627 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 627) )
  if (anyExceptions()) return
#line 628 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine relative permeability as a function of &
             &capillary pressure below polynomial fit'
    value = 1.d0
#line 631 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 631) )
  if (anyExceptions()) return
#line 632 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of saturation as a function &
             &of capillary pressure below polynomial fit'
    value = 0.d0
#line 635 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 635) )
  if (anyExceptions()) return
#line 636 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of relative permeability &
             &as a function of pressure below polynomial fit'
    value = 0.d0
#line 639 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 639) )
  if (anyExceptions()) return
#line 640 "test_characteristic_curves.pf"

    ! saturation = f(capillary_pressure) within the polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_BC_type)
        capillary_pressure = 0.96d0/sf%alpha
      class default
        print *, 'not bc type in testsf_Brooks_Corey'
    end select
    call this%cc_bcb%saturation_function% &
         Saturation(capillary_pressure, &
                    liquid_saturation, dsat_pres,this%option)
    call this%cc_bcb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'Brooks-Corey-Burdine saturation as a function of capillary &
             &pressure within polynomial fit'
    value = 0.99971176979312304d0
#line 658 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 658) )
  if (anyExceptions()) return
#line 659 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine relative permeability as a function of &
             &capillary pressure within polynomial fit'
    value = 0.99789158871529349d0
#line 662 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 662) )
  if (anyExceptions()) return
#line 663 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of saturation as a function &
             &of capillary pressure within polynomial fit'
    value = 5.6675690490728353d-7
#line 666 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 666) )
  if (anyExceptions()) return
#line 667 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of relative permeability &
             &as a function of pressure within polynomial fit'
    value = 4.1422137957785640d-006
#line 670 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 670) )
  if (anyExceptions()) return
#line 671 "test_characteristic_curves.pf"

    ! saturation = f(capillary_pressure) above the polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_BC_type)
        capillary_pressure = 1.06d0/sf%alpha
      class default
        print *, 'not bc type in testsf_Brooks_Corey'
    end select
    call this%cc_bcb%saturation_function% &
         Saturation(capillary_pressure, &
                    liquid_saturation, dsat_pres,this%option)
    call this%cc_bcb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'Brooks-Corey-Burdine saturation as a function of capillary &
             &pressure above polynomial fit'
    value = 0.96802592722174041d0
#line 689 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 689) )
  if (anyExceptions()) return
#line 690 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine relative permeability as a function of &
             &capillary pressure above polynomial fit'
    value = 0.78749164071142996d0
#line 693 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 693) )
  if (anyExceptions()) return
#line 694 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of saturation as a function &
             &of capillary pressure above polynomial fit'
    value = 5.0054278424773111d-6
#line 697 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 697) )
  if (anyExceptions()) return
#line 698 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of relative permeability &
             &as a function of pressure above polynomial fit'
    value = 3.0060561800889172d-005
#line 701 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 701) )
  if (anyExceptions()) return
#line 702 "test_characteristic_curves.pf"

  end subroutine testsf_Brooks_Corey

! ************************************************************************** !

!  @Test
  subroutine testcp_Brooks_Corey(this)

    implicit none

    class (Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
    PetscReal, parameter :: temperature = 25.d0

    ! capillary pressure = f(saturation) well within polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_BC_type)
        liquid_saturation = 1.00001d0**(-sf%lambda)
      class default
        print *, 'not bc type in capillary pressure Brooks Corey'
    end select
    call this%cc_bcb%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = 'Brooks-Corey capillary pressure as a function of &
             &liquid saturation barely within polynomial fit'
    value = 54.068777590990067d0
#line 735 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 735) )
  if (anyExceptions()) return
#line 736 "test_characteristic_curves.pf"
    string = 'Brooks-Corey derivative of capillary pressure as a function of &
             &liquid saturation barely within polynomial fit'
    value = -7.7231957793806121d6
#line 739 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 739) )
  if (anyExceptions()) return
#line 740 "test_characteristic_curves.pf"

    ! capillary pressure = f(saturation) slightly within polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_BC_type)
        liquid_saturation = 1.04d0**(-sf%lambda)
      class default
        print *, 'not bc type in capillary pressure Brooks Corey'
    end select
    call this%cc_bcb%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = 'Brooks-Corey capillary pressure as a function of &
             &liquid saturation well within polynomial fit'
    value = 106436.99642977261d0
#line 754 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 754) )
  if (anyExceptions()) return
#line 755 "test_characteristic_curves.pf"
    string = 'Brooks-Corey derivative of capillary pressure as a function of &
             &liquid saturation well within polynomial fit'
    value = -1.9672548074671443d5
#line 758 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 758) )
  if (anyExceptions()) return
#line 759 "test_characteristic_curves.pf"

    ! capillary pressure = f(saturation) above polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_BC_type)
        liquid_saturation = 1.06d0**(-sf%lambda)
      class default
        print *, 'not bc type in capillary pressure Brooks Corey'
    end select
    call this%cc_bcb%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = 'Brooks-Corey capillary pressure as a function of &
             &liquid saturation above within polynomial fit'
    value = 109024.42772683989d0
#line 773 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 773) )
  if (anyExceptions()) return
#line 774 "test_characteristic_curves.pf"
    string = 'Brooks-Corey derivative of capillary pressure as a function of &
             &liquid saturation above within polynomial fit'
    value = -2.0492439614025093d5
#line 777 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 777) )
  if (anyExceptions()) return
#line 778 "test_characteristic_curves.pf"

  end subroutine testcp_Brooks_Corey

! ************************************************************************** !

!  @Test
  subroutine testrpf_BC_Burdine(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: relative_permeability
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_bcb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'Brooks-Corey-Burdine liquid relative permeability as a &
             &function of liquid saturation'
    value = 3.1991918327000197d-3
#line 805 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 805) )
  if (anyExceptions()) return
#line 806 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 6.2460411971762310d-2
#line 809 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 809) )
  if (anyExceptions()) return
#line 810 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_bcb%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Brooks-Corey-Burdine gas relative permeability as a &
             &function of liquid saturation'
    value = 0.38173220142506209d0
#line 818 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 818) )
  if (anyExceptions()) return
#line 819 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.6412199910843073d0
#line 822 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 822) )
  if (anyExceptions()) return
#line 823 "test_characteristic_curves.pf"

  end subroutine testrpf_BC_Burdine

! ************************************************************************** !

!  @Test
  subroutine testrpf_BC_Mualem(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: relative_permeability
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_bcm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'Brooks-Corey-Mualem liquid relative permeability as a &
             &function of liquid saturation'
    value = 5.2242583862629442d-3
#line 850 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 850) )
  if (anyExceptions()) return
#line 851 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 9.3290328326124022d-2
#line 854 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 854) )
  if (anyExceptions()) return
#line 855 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_bcm%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Brooks-Corey-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.65126653365343257d0
#line 863 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 863) )
  if (anyExceptions()) return
#line 864 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Mualem derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.7243442005604193d0
#line 867 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 867) )
  if (anyExceptions()) return
#line 868 "test_characteristic_curves.pf"

  end subroutine testrpf_BC_Mualem

! ************************************************************************** !

!  @Test
  subroutine testsf_van_Genuchten(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: dkr_p
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! saturation = f(capillary_pressure) at low capillary pressure    
    capillary_pressure = 10.d0
    call this%cc_vgm%saturation_function%Saturation(capillary_pressure, &
                                         liquid_saturation, &
                                         dsat_pres,this%option)
    call this%cc_vgm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'van Genuchten-Mualem saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.99999995045230206d0
#line 902 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 902) )
  if (anyExceptions()) return
#line 903 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.99957025105913566d0
#line 906 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 906) )
  if (anyExceptions()) return
#line 907 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of saturation as a function &
             &of capillary pressure at low capillary pressure'
    value = 1.0475199529417896d-8
#line 910 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 910) )
  if (anyExceptions()) return
#line 911 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of relative permeability &
             &as a function of capillary pressure at low capillary pressure'
    value = 4.7878857031474202d-005
#line 914 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 914) )
  if (anyExceptions()) return
#line 915 "test_characteristic_curves.pf"

    ! saturation = f(capillary_pressure) at high capillary pressure
    select type(sf=>this%cc_vgm%saturation_function)
      class is(sat_func_VG_type)
        capillary_pressure = 10.d0/sf%get_alpha()
      class default
        print *, 'not vg type in testsf_van_Genuchten'
    end select
    call this%cc_vgm%saturation_function%Saturation(capillary_pressure, &
                                         liquid_saturation, &
                                         dsat_pres,this%option)
    call this%cc_vgm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'van Genuchten-Mualem saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.20862404282784081d0
#line 933 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 933) )
  if (anyExceptions()) return
#line 934 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 4.4900562293186444d-6
#line 937 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 937) )
  if (anyExceptions()) return
#line 938 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of saturation as a function &
             &of capillary pressure at high capillary pressure'
    value = 3.7043838142841442d-7
#line 941 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 941) )
  if (anyExceptions()) return
#line 942 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of relative permeability &
             &as a function of capillary pressure at high capillary pressure'
    value = 1.0903614584398361d-010
#line 945 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 945) )
  if (anyExceptions()) return
#line 946 "test_characteristic_curves.pf"

  end subroutine testsf_van_Genuchten

! ************************************************************************** !

!  @Test
  subroutine testcp_van_Genuchten(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
    PetscReal, parameter :: temperature = 25.d0

    ! capillary pressure = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_vgm%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = 'van Genuchten-Mualem capillary pressure as a function of &
             &saturation'
    value = 38910.985405751228d0
#line 974 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 974) )
  if (anyExceptions()) return
#line 975 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of capillary pressure as a &
              &function of saturation'
    value = -1.2074621994359563d5
#line 978 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 978) )
  if (anyExceptions()) return
#line 979 "test_characteristic_curves.pf"

  end subroutine testcp_van_Genuchten

! ************************************************************************** !

!  @Test
  subroutine testrpf_van_Genuchten_Mualem(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_vgm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'van Genuchten-Mualem liquid relative permeability as a &
             &function of liquid saturation'
    value = 7.1160141309814171d-3
#line 1006 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1006) )
  if (anyExceptions()) return
#line 1007 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 8.9580035202641822d-2
#line 1010 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1010) )
  if (anyExceptions()) return
#line 1011 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_vgm%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'van Genuchten-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.61184154078016839d0
#line 1019 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1019) )
  if (anyExceptions()) return
#line 1020 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.4149310375033495d0
#line 1023 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1023) )
  if (anyExceptions()) return
#line 1024 "test_characteristic_curves.pf"

  end subroutine testrpf_van_Genuchten_Mualem

! ************************************************************************** !

!  @Test
  subroutine testrpf_van_Genuchten_Burdine(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_vgb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'van Genuchten-Burdine liquid relative permeability as a &
             &function of liquid saturation'
    value = 1.8220963608953099d-2
#line 1051 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1051) )
  if (anyExceptions()) return
#line 1052 "test_characteristic_curves.pf"
    string = 'van Genuchten-Burdine derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 0.20400586616752553d0
#line 1055 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1055) )
  if (anyExceptions()) return
#line 1056 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_vgb%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'van Genuchten-Burdine gas relative permeability as a &
             &function of liquid saturation'
    value = 0.29870096712333277d0
#line 1064 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1064) )
  if (anyExceptions()) return
#line 1065 "test_characteristic_curves.pf"
    string = 'van Genuchten-Burdine derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.4207000510364920d0
#line 1068 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1068) )
  if (anyExceptions()) return
#line 1069 "test_characteristic_curves.pf"

  end subroutine testrpf_van_Genuchten_Burdine

! ************************************************************************** !

!  @Test
  subroutine testrpf_TOUGH2_IRP7_gas(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0

    ! gas relative permeability = f(saturation)
    call this%cc_vgt%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'TOUGH2 IRP7 gas relative permeability as a &
             &function of liquid saturation'
    value = 0.27522069402853439d0
#line 1098 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1098) )
  if (anyExceptions()) return
#line 1099 "test_characteristic_curves.pf"
    string = 'TOUGH2 IRP7 derivative of gas relative permeability as a &
             &function of liquid saturation'
    value = -1.4564360410147879d0
#line 1102 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1102) )
  if (anyExceptions()) return
#line 1103 "test_characteristic_curves.pf"
  end subroutine testrpf_TOUGH2_IRP7_gas

! ************************************************************************** !

!  @Test
  subroutine testsf_Linear(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: dkr_p
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! saturation = f(capillary_pressure) at low capillary pressure    
    capillary_pressure = 10.d0
    call this%cc_lm%saturation_function%Saturation(capillary_pressure, &
                                         liquid_saturation, &
                                         dsat_pres,this%option)
    call this%cc_lm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'Linear-Mualem saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0001679766584224d0
#line 1136 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1136) )
  if (anyExceptions()) return
#line 1137 "test_characteristic_curves.pf"
    string = 'Linear-Mualem relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0000000000000000d0
#line 1140 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1140) )
  if (anyExceptions()) return
#line 1141 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of saturation as a function &
             &of capillary pressure at low capillary pressure'
    value = 8.5802608854958100d-9
#line 1144 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1144) )
  if (anyExceptions()) return
#line 1145 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of relative permeability as a &
             &function of capillary pressure at low capillary pressure'
    value = 0.d0
#line 1148 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1148) )
  if (anyExceptions()) return
#line 1149 "test_characteristic_curves.pf"

    ! saturation = f(capillary_pressure) at high capillary pressure
    select type(sf=>this%cc_lm%saturation_function)
      class is(sat_func_Linear_type)
        capillary_pressure = 10.d0/sf%alpha
      class default
        print *, 'not linear type in testsf_Linear'
    end select
    call this%cc_lm%saturation_function%Saturation(capillary_pressure, &
                                         liquid_saturation, &
                                         dsat_pres,this%option)
    call this%cc_lm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'van Genuchten-Mualem saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.99848743785071759d0
#line 1167 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1167) )
  if (anyExceptions()) return
#line 1168 "test_characteristic_curves.pf"
    string = 'Linear-Mualem relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.21040639641042236d0
#line 1171 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1171) )
  if (anyExceptions()) return
#line 1172 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of saturation as a function &
             &of capillary pressure at high capillary pressure'
    value = 8.5802608854958100d-9
#line 1175 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1175) )
  if (anyExceptions()) return
#line 1176 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of relative permeability as a &
             &function of capillary pressure at high capillary pressure'
    value = 3.7741597839123353d-007
#line 1179 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1179) )
  if (anyExceptions()) return
#line 1180 "test_characteristic_curves.pf"

  end subroutine testsf_Linear

! ************************************************************************** !

!  @Test
  subroutine testcp_Linear(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
    PetscReal, parameter :: temperature = 25.d0

    ! capillary pressure = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_lm%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = 'Linear capillary pressure as a function of &
             &saturation'
    value = 58292873.507671818d0
#line 1208 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1208) )
  if (anyExceptions()) return
#line 1209 "test_characteristic_curves.pf"
    string = 'Linear derivative of capillary pressure as a function of &
             &saturation'
    value = -1.1654657280764198d8
#line 1212 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1212) )
  if (anyExceptions()) return
#line 1213 "test_characteristic_curves.pf"

  end subroutine testcp_Linear

! ************************************************************************** !

!  @Test
  subroutine testrpf_Linear_Mualem(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_lm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'Linear-Mualem liquid relative permeability as a &
             &function of liquid saturation'
    value = 9.8205932575543323d-4
#line 1240 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1240) )
  if (anyExceptions()) return
#line 1241 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 8.6657420154635095d-3
#line 1244 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1244) )
  if (anyExceptions()) return
#line 1245 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_lm%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'Linear-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.70258636467899827d0
#line 1253 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1253) )
  if (anyExceptions()) return
#line 1254 "test_characteristic_curves.pf"

  end subroutine testrpf_Linear_Mualem

! ************************************************************************** !

!  @Test
  subroutine testrpf_Linear_Burdine(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_lb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'Linear-Mualem liquid relative permeability as a &
             &function of liquid saturation'
    value = 0.41656942823803966d0
#line 1281 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1281) )
  if (anyExceptions()) return
#line 1282 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 1.1668611435239207d0
#line 1285 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1285) )
  if (anyExceptions()) return
#line 1286 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_lb%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'Linear-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.57851239669421495d0
#line 1294 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1294) )
  if (anyExceptions()) return
#line 1295 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = -1.1806375442739079d0
#line 1298 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1298) )
  if (anyExceptions()) return
#line 1299 "test_characteristic_curves.pf"

  end subroutine testrpf_Linear_Burdine

! ************************************************************************** !

!  @Test
  subroutine testsf_modified_kosugi_3param(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: dkr_p
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! saturation = f(capillary_pressure)
    capillary_pressure = 1.0D+4
    call this%cc_mk3%saturation_function%Saturation(capillary_pressure, &
                                         liquid_saturation, &
                                         dsat_pres,this%option)
    ! relative permeability = f(saturation)
    call this%cc_mk3%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = '3-parameter modified Kosugi saturation as a function of '//&
         &'capillary pressure'
    value = 0.97481347586415668d0
#line 1334 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1334) )
  if (anyExceptions()) return
#line 1335 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi relative permeability as a '//&
         &'function of capillary pressure'
    value = 0.92508223462292838d0
#line 1338 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1338) )
  if (anyExceptions()) return
#line 1339 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of saturation as '//&
         &'a function of capillary pressure'
    value = 4.7607337955577978D-5
#line 1342 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1342) )
  if (anyExceptions()) return
#line 1343 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of relative perm as'// &
         &' a function of capillary pressure'
    value = 1.2551629198938072D-4
#line 1346 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1346) )
  if (anyExceptions()) return
#line 1347 "test_characteristic_curves.pf"

  end subroutine testsf_modified_kosugi_3param

! ************************************************************************** !

!  @Test
  subroutine testcp_modified_kosugi_3param(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
    PetscReal, parameter :: temperature = 25.d0

    ! capillary pressure = f(saturation)
    liquid_saturation = 2.5d-1
    call this%cc_mk3%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = '3-parameter modified Kosugi capillary pressure as a '//&
         &'function of saturation'
    value = 17707.010438883957d0
#line 1375 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1375) )
  if (anyExceptions()) return
#line 1376 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of capillary ' // &
         'pressure as a function of saturation'
    value = -24500.406539923177d0
#line 1379 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1379) )
  if (anyExceptions()) return
#line 1380 "test_characteristic_curves.pf"

  end subroutine testcp_modified_kosugi_3param

! ************************************************************************** !

!  @Test
  subroutine testrpf_modified_kosugi_3param(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_mk3%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = '3-parameter modified Kosugi liquid relative permeability'//&
         &' as a function of liquid saturation'
    value = 0.18300229972367688d0
#line 1407 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1407) )
  if (anyExceptions()) return
#line 1408 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of liquid relative'//&
         &' permeability as a function of liquid saturation'
    value = 0.92477440290277491d0
#line 1411 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1411) )
  if (anyExceptions()) return
#line 1412 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_mk3%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = '3-parameter modified Kosugi gas relative permeability as'//&
         &' a function of liquid saturation'
    value = 0.35040465818058608d0
#line 1420 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1420) )
  if (anyExceptions()) return
#line 1421 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of gas relative'//&
         &' permeability as a function of liquid saturation'
    value = -1.2770282015052450d0
#line 1424 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1424) )
  if (anyExceptions()) return
#line 1425 "test_characteristic_curves.pf"

  end subroutine testrpf_modified_kosugi_3param

! ************************************************************************** !

!  @Test
  subroutine testcp_modified_kosugi_4param(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
    PetscReal, parameter :: temperature = 25.d0

    ! capillary pressure = f(saturation)
    liquid_saturation = 2.5d-1
    call this%cc_mk4%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = '4-parameter modified Kosugi capillary pressure as a '//&
         &'function of saturation'
    value = 15671.942098006928d0
#line 1453 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1453) )
  if (anyExceptions()) return
#line 1454 "test_characteristic_curves.pf"
    string = 'r-parameter modified Kosugi derivative of capillary ' // &
         'pressure as a function of saturation'
    value = -19192.362637570313d0
#line 1457 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1457) )
  if (anyExceptions()) return
#line 1458 "test_characteristic_curves.pf"

  end subroutine testcp_modified_kosugi_4param

! ************************************************************************** !

!  @Test
  subroutine testrpf_modified_kosugi_4param(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_mk4%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = '4-parameter modified Kosugi liquid relative permeability'//&
         &' as a function of liquid saturation'
    value = 0.18300229972367688d0  ! same as 3-param
#line 1485 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1485) )
  if (anyExceptions()) return
#line 1486 "test_characteristic_curves.pf"
    string = '4-parameter modified Kosugi derivative of liquid relative'//&
         &' permeability as a function of liquid saturation'
    value = 0.92477440290277491d0  ! same as 3-param
#line 1489 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1489) )
  if (anyExceptions()) return
#line 1490 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_mk4%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = '4-parameter modified Kosugi gas relative permeability as'//&
         &' a function of liquid saturation'
    value = 0.35040465818058608d0  ! same as 3-param
#line 1498 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1498) )
  if (anyExceptions()) return
#line 1499 "test_characteristic_curves.pf"
    string = '4-parameter modified Kosugi derivative of gas relative'//&
         &' permeability as a function of liquid saturation'
    value = -1.2770282015052450d0   ! same as 3-param
#line 1502 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1502) )
  if (anyExceptions()) return
#line 1503 "test_characteristic_curves.pf"

  end subroutine testrpf_modified_kosugi_4param

! ************************************************************************** !

!  @Test
  subroutine testsf_constant(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: dsat_pres
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! saturation = f(capillary_pressure)
    capillary_pressure = 1.d5
    call this%cc_const%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'Constant saturation function liquid saturation'
    value = 0.5d0
#line 1530 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1530) )
  if (anyExceptions()) return
#line 1531 "test_characteristic_curves.pf"
    string = 'Constant saturation function liquid saturation derivative'
    value = 0.d0
#line 1533 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1533) )
  if (anyExceptions()) return
#line 1534 "test_characteristic_curves.pf"

    ! capillary_pressure = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_const%saturation_function% &
                                    CapillaryPressure(liquid_saturation, &
                                                      capillary_pressure, &
                                                      dpc_dsatl, &
                                                      this%option)
    string = 'Constant saturation function capillary pressure'
    value = 1.d5
#line 1544 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1544) )
  if (anyExceptions()) return
#line 1545 "test_characteristic_curves.pf"

  end subroutine testsf_constant
  
! ************************************************************************** !

!  @Test
  subroutine testsf_KRP1(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.22
    call this%cc_krp1%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP1 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 211645.56118006655d0
#line 1574 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1574) )
  if (anyExceptions()) return
#line 1575 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp1%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP1 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.74556316913866438d0
#line 1582 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1582) )
  if (anyExceptions()) return
#line 1583 "test_characteristic_curves.pf"
    call this%cc_krp1%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.10878895018640269d0
#line 1589 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1589) )
  if (anyExceptions()) return
#line 1590 "test_characteristic_curves.pf"
    call this%cc_krp1%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 5.9765155008000600d-4
#line 1596 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1596) )
  if (anyExceptions()) return
#line 1597 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp1%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP1 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 0.0d0
#line 1605 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1605) )
  if (anyExceptions()) return
#line 1606 "test_characteristic_curves.pf"
    select type(sf=>this%cc_krp1%saturation_function)
      class is(sat_func_KRP1_type)
        capillary_pressure = 10.d0/sf%alpha
    end select
    call this%cc_krp1%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP1 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.20000256369851988d0
#line 1616 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1616) )
  if (anyExceptions()) return
#line 1617 "test_characteristic_curves.pf"
    call this%cc_krp1%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 3.1269006686492169d-22
#line 1623 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1623) )
  if (anyExceptions()) return
#line 1624 "test_characteristic_curves.pf"
    call this%cc_krp1%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.0d0
#line 1630 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1630) )
  if (anyExceptions()) return
#line 1631 "test_characteristic_curves.pf"

  end subroutine testsf_KRP1
  
! ************************************************************************** !

!  @Test
  subroutine testsf_KRP2(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.05
    call this%cc_krp2%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP2 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 9990000.0d0
#line 1660 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1660) )
  if (anyExceptions()) return
#line 1661 "test_characteristic_curves.pf"
    capillary_pressure = 100.d0
    call this%cc_krp2%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP2 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 1668 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1668) )
  if (anyExceptions()) return
#line 1669 "test_characteristic_curves.pf"
    call this%cc_krp2%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1675 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1675) )
  if (anyExceptions()) return
#line 1676 "test_characteristic_curves.pf"
    call this%cc_krp2%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1682 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1682) )
  if (anyExceptions()) return
#line 1683 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp2%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP2 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 51637.669399930157d0
#line 1691 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1691) )
  if (anyExceptions()) return
#line 1692 "test_characteristic_curves.pf"
    select type(sf=>this%cc_krp2%saturation_function)
      class is(sat_func_KRP2_type)
        capillary_pressure = 15.d0/sf%alpha
    end select
    call this%cc_krp2%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP2 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.55731947333915333d0
#line 1701 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1701) )
  if (anyExceptions()) return
#line 1702 "test_characteristic_curves.pf"
    call this%cc_krp2%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 5.8310805074531403d-4
#line 1708 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1708) )
  if (anyExceptions()) return
#line 1709 "test_characteristic_curves.pf"
    call this%cc_krp2%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.24138701885980715d0
#line 1715 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1715) )
  if (anyExceptions()) return
#line 1716 "test_characteristic_curves.pf"

  end subroutine testsf_KRP2
  
! ************************************************************************** !

!  @Test
  subroutine testsf_KRP3(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.21
    call this%cc_krp3%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP3 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 143937.17953954436d0
#line 1745 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1745) )
  if (anyExceptions()) return
#line 1746 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp3%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP3 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 1752 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1752) )
  if (anyExceptions()) return
#line 1753 "test_characteristic_curves.pf"
    call this%cc_krp3%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1759 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1759) )
  if (anyExceptions()) return
#line 1760 "test_characteristic_curves.pf"
    call this%cc_krp3%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1766 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1766) )
  if (anyExceptions()) return
#line 1767 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp3%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP3 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 200.0d0
#line 1775 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1775) )
  if (anyExceptions()) return
#line 1776 "test_characteristic_curves.pf"
    capillary_pressure = 1.d5
    call this%cc_krp3%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP3 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.21815720128839766d0
#line 1782 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1782) )
  if (anyExceptions()) return
#line 1783 "test_characteristic_curves.pf"
    call this%cc_krp3%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 5.8632915084341688d-9
#line 1789 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1789) )
  if (anyExceptions()) return
#line 1790 "test_characteristic_curves.pf"
    call this%cc_krp3%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.78571287226820719d0
#line 1796 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1796) )
  if (anyExceptions()) return
#line 1797 "test_characteristic_curves.pf"

  end subroutine testsf_KRP3
  
! ************************************************************************** !

!  @Test
  subroutine testsf_KRP4(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.16
    call this%cc_krp4%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP4 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 9990000.0d0
#line 1826 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1826) )
  if (anyExceptions()) return
#line 1827 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp4%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP4 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 1833 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1833) )
  if (anyExceptions()) return
#line 1834 "test_characteristic_curves.pf"
    call this%cc_krp4%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1840 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1840) )
  if (anyExceptions()) return
#line 1841 "test_characteristic_curves.pf"
    call this%cc_krp4%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1847 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1847) )
  if (anyExceptions()) return
#line 1848 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp4%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP4 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 152.32289193581161d0
#line 1856 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1856) )
  if (anyExceptions()) return
#line 1857 "test_characteristic_curves.pf"
    capillary_pressure = 1.d5
    call this%cc_krp4%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP4 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.21815720128839766d0
#line 1863 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1863) )
  if (anyExceptions()) return
#line 1864 "test_characteristic_curves.pf"
    call this%cc_krp4%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 2.8180697920854255d-10
#line 1870 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1870) )
  if (anyExceptions()) return
#line 1871 "test_characteristic_curves.pf"
    call this%cc_krp4%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.78571287226820719d0
#line 1877 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1877) )
  if (anyExceptions()) return
#line 1878 "test_characteristic_curves.pf"

  end subroutine testsf_KRP4
  
! ************************************************************************** !

!  @Test
  subroutine testsf_KRP5(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.16
    call this%cc_krp5%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP5 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 8525093.3809598293d0
#line 1907 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1907) )
  if (anyExceptions()) return
#line 1908 "test_characteristic_curves.pf"
    capillary_pressure = 50.d0
    call this%cc_krp5%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP5 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.80014642571085304d0
#line 1915 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1915) )
  if (anyExceptions()) return
#line 1916 "test_characteristic_curves.pf"
    call this%cc_krp5%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1922 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1922) )
  if (anyExceptions()) return
#line 1923 "test_characteristic_curves.pf"
    call this%cc_krp5%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1929 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1929) )
  if (anyExceptions()) return
#line 1930 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp5%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP5 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 2000.0d0
#line 1938 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1938) )
  if (anyExceptions()) return
#line 1939 "test_characteristic_curves.pf"
    capillary_pressure = 8.d5
    call this%cc_krp5%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP5 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.74007809371245503d0
#line 1946 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1946) )
  if (anyExceptions()) return
#line 1947 "test_characteristic_curves.pf"
    call this%cc_krp5%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.92010412494993998d0
#line 1953 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1953) )
  if (anyExceptions()) return
#line 1954 "test_characteristic_curves.pf"
    call this%cc_krp5%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 7.9895875050060017d-2
#line 1960 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1960) )
  if (anyExceptions()) return
#line 1961 "test_characteristic_curves.pf"

  end subroutine testsf_KRP5
  
! ************************************************************************** !

!  @Test
  subroutine testsf_KRP8(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.16
    call this%cc_krp8%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP8 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 390287.05422819848d0
#line 1990 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1990) )
  if (anyExceptions()) return
#line 1991 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp8%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP8 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.99813934029313034d0
#line 1998 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1998) )
  if (anyExceptions()) return
#line 1999 "test_characteristic_curves.pf"
    call this%cc_krp8%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.82967396002016347d0
#line 2005 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2005) )
  if (anyExceptions()) return
#line 2006 "test_characteristic_curves.pf"
    call this%cc_krp8%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 3.5744473041460765d-4
#line 2012 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2012) )
  if (anyExceptions()) return
#line 2013 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp8%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP8 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 11673.176741409796d0
#line 2021 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2021) )
  if (anyExceptions()) return
#line 2022 "test_characteristic_curves.pf"
    select type(sf=>this%cc_krp8%saturation_function)
      class is(sat_func_KRP8_type)
        capillary_pressure = 15.d0/sf%alpha
    end select
    call this%cc_krp8%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP8 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.10220037076344601d0
#line 2032 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2032) )
  if (anyExceptions()) return
#line 2033 "test_characteristic_curves.pf"
    call this%cc_krp8%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 2.4705248145444649d-14
#line 2039 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2039) )
  if (anyExceptions()) return
#line 2040 "test_characteristic_curves.pf"
    call this%cc_krp8%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.99877541173464723d0
#line 2046 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2046) )
  if (anyExceptions()) return
#line 2047 "test_characteristic_curves.pf"

  end subroutine testsf_KRP8
  
! ************************************************************************** !

!  @Test
  subroutine testsf_KRP9(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.16
    call this%cc_krp9%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP9 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 6701.4503478003835d0
#line 2076 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2076) )
  if (anyExceptions()) return
#line 2077 "test_characteristic_curves.pf"
    capillary_pressure = 50.d0
    call this%cc_krp9%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP9 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.99999644138299326d0
#line 2084 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2084) )
  if (anyExceptions()) return
#line 2085 "test_characteristic_curves.pf"
    call this%cc_krp9%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.99999998839687010d0
#line 2091 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2091) )
  if (anyExceptions()) return
#line 2092 "test_characteristic_curves.pf"
    call this%cc_krp9%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.1603129901338605d-8
#line 2098 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2098) )
  if (anyExceptions()) return
#line 2099 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp9%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP9 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 2294.5062171637132d0
#line 2107 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2107) )
  if (anyExceptions()) return
#line 2108 "test_characteristic_curves.pf"
    select type(sf=>this%cc_krp9%saturation_function)
      class is(sat_func_KRP9_type)
        capillary_pressure = 8.d5
    end select
    call this%cc_krp9%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP9 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 1.8062138845808708d-7
#line 2118 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2118) )
  if (anyExceptions()) return
#line 2119 "test_characteristic_curves.pf"
    call this%cc_krp9%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.d0
#line 2125 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2125) )
  if (anyExceptions()) return
#line 2126 "test_characteristic_curves.pf"
    call this%cc_krp9%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.d0
#line 2132 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2132) )
  if (anyExceptions()) return
#line 2133 "test_characteristic_curves.pf"

  end subroutine testsf_KRP9
  
! ************************************************************************** !

!  @Test
  subroutine testsf_KRP11(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
 
    liquid_saturation = 0.15
    call this%cc_krp11%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP11 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 0.0d0
#line 2162 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2162) )
  if (anyExceptions()) return
#line 2163 "test_characteristic_curves.pf"
    capillary_pressure = 100.d0
    call this%cc_krp11%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP11 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 2170 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2170) )
  if (anyExceptions()) return
#line 2171 "test_characteristic_curves.pf"
    call this%cc_krp11%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 2177 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2177) )
  if (anyExceptions()) return
#line 2178 "test_characteristic_curves.pf"
    call this%cc_krp11%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.d0
#line 2184 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2184) )
  if (anyExceptions()) return
#line 2185 "test_characteristic_curves.pf"

    liquid_saturation = 0.85
    call this%cc_krp11%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP11 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 0.0d0
#line 2193 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2193) )
  if (anyExceptions()) return
#line 2194 "test_characteristic_curves.pf"
    capillary_pressure = 1.d6
    call this%cc_krp11%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP11 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 1.0d0
#line 2201 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2201) )
  if (anyExceptions()) return
#line 2202 "test_characteristic_curves.pf"
    call this%cc_krp11%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.0d0
#line 2208 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2208) )
  if (anyExceptions()) return
#line 2209 "test_characteristic_curves.pf"
    call this%cc_krp11%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.0d0
#line 2215 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2215) )
  if (anyExceptions()) return
#line 2216 "test_characteristic_curves.pf"

  end subroutine testsf_KRP11
  
! ************************************************************************** !

!  @Test
  subroutine testsf_KRP12(this)
    
    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
  
    liquid_saturation = 0.15
    call this%cc_krp12%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP12 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 35364351.440591194d0
#line 2245 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2245) )
  if (anyExceptions()) return
#line 2246 "test_characteristic_curves.pf"
    capillary_pressure = 100.d0
    call this%cc_krp12%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP12 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 2253 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2253) )
  if (anyExceptions()) return
#line 2254 "test_characteristic_curves.pf"
    call this%cc_krp12%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 2260 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2260) )
  if (anyExceptions()) return
#line 2261 "test_characteristic_curves.pf"
    call this%cc_krp12%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.d0
#line 2267 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2267) )
  if (anyExceptions()) return
#line 2268 "test_characteristic_curves.pf"

    liquid_saturation = 0.85
    call this%cc_krp12%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP12 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 4547.8850244305122d0
#line 2276 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2276) )
  if (anyExceptions()) return
#line 2277 "test_characteristic_curves.pf"
    capillary_pressure = 1.d6
    call this%cc_krp12%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP12 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.30276918155781385d0
#line 2284 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2284) )
  if (anyExceptions()) return
#line 2285 "test_characteristic_curves.pf"
    call this%cc_krp12%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.3380237859265290d-9
#line 2291 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2291) )
  if (anyExceptions()) return
#line 2292 "test_characteristic_curves.pf"
    call this%cc_krp12%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.53468502607262869d0
#line 2298 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2298) )
  if (anyExceptions()) return
#line 2299 "test_characteristic_curves.pf"

  end subroutine testsf_KRP12

! ************************************************************************** !

end module Test_Characteristic_Curves_module




function Test_Characteristic_Curves_module_suite() result(suite)
   use pFUnit_mod
   use Test_Characteristic_Curves_module
   implicit none

   type (TestSuite) :: suite
   suite = newTestSuite('Test_Characteristic_Curves_module_suite')

   call suite%addTest(Test_Characteristic_Curves('testSF_BC_SetupPolynomials',testSF_BC_SetupPolynomials))
   call suite%addTest(Test_Characteristic_Curves('testsf_Brooks_Corey',testsf_Brooks_Corey))
   call suite%addTest(Test_Characteristic_Curves('testcp_Brooks_Corey',testcp_Brooks_Corey))
   call suite%addTest(Test_Characteristic_Curves('testrpf_BC_Burdine',testrpf_BC_Burdine))
   call suite%addTest(Test_Characteristic_Curves('testrpf_BC_Mualem',testrpf_BC_Mualem))
   call suite%addTest(Test_Characteristic_Curves('testsf_van_Genuchten',testsf_van_Genuchten))
   call suite%addTest(Test_Characteristic_Curves('testcp_van_Genuchten',testcp_van_Genuchten))
   call suite%addTest(Test_Characteristic_Curves('testrpf_van_Genuchten_Mualem',testrpf_van_Genuchten_Mualem))
   call suite%addTest(Test_Characteristic_Curves('testrpf_van_Genuchten_Burdine',testrpf_van_Genuchten_Burdine))
   call suite%addTest(Test_Characteristic_Curves('testrpf_TOUGH2_IRP7_gas',testrpf_TOUGH2_IRP7_gas))
   call suite%addTest(Test_Characteristic_Curves('testsf_Linear',testsf_Linear))
   call suite%addTest(Test_Characteristic_Curves('testcp_Linear',testcp_Linear))
   call suite%addTest(Test_Characteristic_Curves('testrpf_Linear_Mualem',testrpf_Linear_Mualem))
   call suite%addTest(Test_Characteristic_Curves('testrpf_Linear_Burdine',testrpf_Linear_Burdine))
   call suite%addTest(Test_Characteristic_Curves('testsf_modified_kosugi_3param',testsf_modified_kosugi_3param))
   call suite%addTest(Test_Characteristic_Curves('testcp_modified_kosugi_3param',testcp_modified_kosugi_3param))
   call suite%addTest(Test_Characteristic_Curves('testrpf_modified_kosugi_3param',testrpf_modified_kosugi_3param))
   call suite%addTest(Test_Characteristic_Curves('testcp_modified_kosugi_4param',testcp_modified_kosugi_4param))
   call suite%addTest(Test_Characteristic_Curves('testrpf_modified_kosugi_4param',testrpf_modified_kosugi_4param))
   call suite%addTest(Test_Characteristic_Curves('testsf_constant',testsf_constant))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP1',testsf_KRP1))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP2',testsf_KRP2))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP3',testsf_KRP3))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP4',testsf_KRP4))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP5',testsf_KRP5))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP8',testsf_KRP8))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP9',testsf_KRP9))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP11',testsf_KRP11))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP12',testsf_KRP12))

end function Test_Characteristic_Curves_module_suite

