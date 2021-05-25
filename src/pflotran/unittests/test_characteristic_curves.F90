module Test_Characteristic_Curves_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use pFUnit_mod
  use Characteristic_Curves_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_Common_module
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
        sf%m = 0.527d0
        sf%alpha = 5.1054d-5
      class default
        print *, 'not VG type in setUp'
    end select
    this%cc_vgm%liq_rel_perm_function => RPFMualemVGLiqCreate()
    this%cc_vgm%liq_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgm%liq_rel_perm_function)
      class is(rpf_Mualem_VG_liq_type)
        sf%m = 0.527d0
      class default
        print *, 'not Mualem VG Liq type'
    end select
    this%cc_vgm%gas_rel_perm_function => RPFMualemVGGasCreate()
    this%cc_vgm%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgm%gas_rel_perm_function)
      class is(rpf_Mualem_VG_gas_type)
        sf%m = 0.527d0
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
        sf%m = 0.527d0
      class default
        print *, 'not Burdine VG Liq type'
    end select
    this%cc_vgb%gas_rel_perm_function => RPFBurdineVGGasCreate()
    this%cc_vgb%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgb%gas_rel_perm_function)
      class is(rpf_Burdine_VG_gas_type)
        sf%m = 0.527d0
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
#line 574 "test_characteristic_curves.pf"
  call assertEqual(values(i), this%cc_bcb%saturation_function%pres_poly%coefficients(i), dabs(values(i))*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 574) )
  if (anyExceptions()) return
#line 575 "test_characteristic_curves.pf"
    enddo

    ! saturation spline
    values = [-83508464.603000879d0, 173197055.36650354d0, &
              -89688590.763502657d0, 0.d0]
    do i = 1, 3
      write(string,*) i
      string = 'Brooks-Corey-Burdine saturation spline coefficient #' // &
               trim(adjustl(string))
#line 584 "test_characteristic_curves.pf"
  call assertEqual(values(i), this%cc_bcb%saturation_function%sat_poly%coefficients(i),  dabs(values(i))*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 584) )
  if (anyExceptions()) return
#line 585 "test_characteristic_curves.pf"
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
#line 625 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 625) )
  if (anyExceptions()) return
#line 626 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine relative permeability as a function of &
             &capillary pressure below polynomial fit'
    value = 1.d0
#line 629 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 629) )
  if (anyExceptions()) return
#line 630 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of saturation as a function &
             &of capillary pressure below polynomial fit'
    value = 0.d0
#line 633 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 633) )
  if (anyExceptions()) return
#line 634 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of relative permeability &
             &as a function of pressure below polynomial fit'
    value = 0.d0
#line 637 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 637) )
  if (anyExceptions()) return
#line 638 "test_characteristic_curves.pf"

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
#line 656 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 656) )
  if (anyExceptions()) return
#line 657 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine relative permeability as a function of &
             &capillary pressure within polynomial fit'
    value = 0.99789158871529349d0
#line 660 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 660) )
  if (anyExceptions()) return
#line 661 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of saturation as a function &
             &of capillary pressure within polynomial fit'
    value = 5.6675690490728353d-7
#line 664 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 664) )
  if (anyExceptions()) return
#line 665 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of relative permeability &
             &as a function of pressure within polynomial fit'
    value = 4.1422137957785640d-006
#line 668 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 668) )
  if (anyExceptions()) return
#line 669 "test_characteristic_curves.pf"

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
#line 687 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 687) )
  if (anyExceptions()) return
#line 688 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine relative permeability as a function of &
             &capillary pressure above polynomial fit'
    value = 0.78749164071142996d0
#line 691 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 691) )
  if (anyExceptions()) return
#line 692 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of saturation as a function &
             &of capillary pressure above polynomial fit'
    value = 5.0054278424773111d-6
#line 695 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 695) )
  if (anyExceptions()) return
#line 696 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of relative permeability &
             &as a function of pressure above polynomial fit'
    value = 3.0060561800889172d-005
#line 699 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 699) )
  if (anyExceptions()) return
#line 700 "test_characteristic_curves.pf"

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
#line 733 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 733) )
  if (anyExceptions()) return
#line 734 "test_characteristic_curves.pf"
    string = 'Brooks-Corey derivative of capillary pressure as a function of &
             &liquid saturation barely within polynomial fit'
    value = -7.7231957793806121d6
#line 737 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 737) )
  if (anyExceptions()) return
#line 738 "test_characteristic_curves.pf"

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
#line 752 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 752) )
  if (anyExceptions()) return
#line 753 "test_characteristic_curves.pf"
    string = 'Brooks-Corey derivative of capillary pressure as a function of &
             &liquid saturation well within polynomial fit'
    value = -1.9672548074671443d5
#line 756 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 756) )
  if (anyExceptions()) return
#line 757 "test_characteristic_curves.pf"

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
#line 771 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 771) )
  if (anyExceptions()) return
#line 772 "test_characteristic_curves.pf"
    string = 'Brooks-Corey derivative of capillary pressure as a function of &
             &liquid saturation above within polynomial fit'
    value = -2.0492439614025093d5
#line 775 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 775) )
  if (anyExceptions()) return
#line 776 "test_characteristic_curves.pf"

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
#line 803 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 803) )
  if (anyExceptions()) return
#line 804 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 6.2460411971762310d-2
#line 807 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 807) )
  if (anyExceptions()) return
#line 808 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_bcb%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Brooks-Corey-Burdine gas relative permeability as a &
             &function of liquid saturation'
    value = 0.38173220142506209d0
#line 816 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 816) )
  if (anyExceptions()) return
#line 817 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.6412199910843073d0
#line 820 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 820) )
  if (anyExceptions()) return
#line 821 "test_characteristic_curves.pf"

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
#line 848 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 848) )
  if (anyExceptions()) return
#line 849 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 9.3290328326124022d-2
#line 852 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 852) )
  if (anyExceptions()) return
#line 853 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_bcm%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Brooks-Corey-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.65126653365343257d0
#line 861 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 861) )
  if (anyExceptions()) return
#line 862 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Mualem derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.7243442005604193d0
#line 865 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 865) )
  if (anyExceptions()) return
#line 866 "test_characteristic_curves.pf"

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
#line 900 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 900) )
  if (anyExceptions()) return
#line 901 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.99957025105913566d0
#line 904 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 904) )
  if (anyExceptions()) return
#line 905 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of saturation as a function &
             &of capillary pressure at low capillary pressure'
    value = 1.0475199529417896d-8
#line 908 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 908) )
  if (anyExceptions()) return
#line 909 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of relative permeability &
             &as a function of capillary pressure at low capillary pressure'
    value = 4.7878857031474202d-005
#line 912 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 912) )
  if (anyExceptions()) return
#line 913 "test_characteristic_curves.pf"

    ! saturation = f(capillary_pressure) at high capillary pressure
    select type(sf=>this%cc_vgm%saturation_function)
      class is(sat_func_VG_type)
        capillary_pressure = 10.d0/sf%alpha
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
#line 931 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 931) )
  if (anyExceptions()) return
#line 932 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 4.4900562293186444d-6
#line 935 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 935) )
  if (anyExceptions()) return
#line 936 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of saturation as a function &
             &of capillary pressure at high capillary pressure'
    value = 3.7043838142841442d-7
#line 939 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 939) )
  if (anyExceptions()) return
#line 940 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of relative permeability &
             &as a function of capillary pressure at high capillary pressure'
    value = 1.0903614584398361d-010
#line 943 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 943) )
  if (anyExceptions()) return
#line 944 "test_characteristic_curves.pf"

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
#line 972 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 972) )
  if (anyExceptions()) return
#line 973 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of capillary pressure as a &
              &function of saturation'
    value = -1.2074621994359563d5
#line 976 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 976) )
  if (anyExceptions()) return
#line 977 "test_characteristic_curves.pf"

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
#line 1004 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1004) )
  if (anyExceptions()) return
#line 1005 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 8.9580035202641822d-2
#line 1008 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1008) )
  if (anyExceptions()) return
#line 1009 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_vgm%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'van Genuchten-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.61184154078016839d0
#line 1017 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1017) )
  if (anyExceptions()) return
#line 1018 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.4149310375033495d0
#line 1021 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1021) )
  if (anyExceptions()) return
#line 1022 "test_characteristic_curves.pf"

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
#line 1049 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1049) )
  if (anyExceptions()) return
#line 1050 "test_characteristic_curves.pf"
    string = 'van Genuchten-Burdine derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 0.20400586616752553d0
#line 1053 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1053) )
  if (anyExceptions()) return
#line 1054 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_vgb%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'van Genuchten-Burdine gas relative permeability as a &
             &function of liquid saturation'
    value = 0.29870096712333277d0
#line 1062 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1062) )
  if (anyExceptions()) return
#line 1063 "test_characteristic_curves.pf"
    string = 'van Genuchten-Burdine derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.4207000510364920d0
#line 1066 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1066) )
  if (anyExceptions()) return
#line 1067 "test_characteristic_curves.pf"

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
#line 1096 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1096) )
  if (anyExceptions()) return
#line 1097 "test_characteristic_curves.pf"
    string = 'TOUGH2 IRP7 derivative of gas relative permeability as a &
             &function of liquid saturation'
    value = -1.4564360410147879d0
#line 1100 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1100) )
  if (anyExceptions()) return
#line 1101 "test_characteristic_curves.pf"
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
#line 1134 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1134) )
  if (anyExceptions()) return
#line 1135 "test_characteristic_curves.pf"
    string = 'Linear-Mualem relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0000000000000000d0
#line 1138 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1138) )
  if (anyExceptions()) return
#line 1139 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of saturation as a function &
             &of capillary pressure at low capillary pressure'
    value = 8.5802608854958100d-9
#line 1142 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1142) )
  if (anyExceptions()) return
#line 1143 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of relative permeability as a &
             &function of capillary pressure at low capillary pressure'
    value = 0.d0
#line 1146 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1146) )
  if (anyExceptions()) return
#line 1147 "test_characteristic_curves.pf"

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
#line 1165 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1165) )
  if (anyExceptions()) return
#line 1166 "test_characteristic_curves.pf"
    string = 'Linear-Mualem relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.21040639641042236d0
#line 1169 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1169) )
  if (anyExceptions()) return
#line 1170 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of saturation as a function &
             &of capillary pressure at high capillary pressure'
    value = 8.5802608854958100d-9
#line 1173 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1173) )
  if (anyExceptions()) return
#line 1174 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of relative permeability as a &
             &function of capillary pressure at high capillary pressure'
    value = 3.7741597839123353d-007
#line 1177 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1177) )
  if (anyExceptions()) return
#line 1178 "test_characteristic_curves.pf"

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
#line 1206 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1206) )
  if (anyExceptions()) return
#line 1207 "test_characteristic_curves.pf"
    string = 'Linear derivative of capillary pressure as a function of &
             &saturation'
    value = -1.1654657280764198d8
#line 1210 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1210) )
  if (anyExceptions()) return
#line 1211 "test_characteristic_curves.pf"

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
#line 1238 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1238) )
  if (anyExceptions()) return
#line 1239 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 8.6657420154635095d-3
#line 1242 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1242) )
  if (anyExceptions()) return
#line 1243 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_lm%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'Linear-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.70258636467899827d0
#line 1251 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1251) )
  if (anyExceptions()) return
#line 1252 "test_characteristic_curves.pf"

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
#line 1279 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1279) )
  if (anyExceptions()) return
#line 1280 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 1.1668611435239207d0
#line 1283 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1283) )
  if (anyExceptions()) return
#line 1284 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_lb%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, & 
                              dkr_sat, this%option)
    string = 'Linear-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.57851239669421495d0
#line 1292 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1292) )
  if (anyExceptions()) return
#line 1293 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = -1.1806375442739079d0
#line 1296 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1296) )
  if (anyExceptions()) return
#line 1297 "test_characteristic_curves.pf"

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
#line 1332 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1332) )
  if (anyExceptions()) return
#line 1333 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi relative permeability as a '//&
         &'function of capillary pressure'
    value = 0.92508223462292838d0
#line 1336 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1336) )
  if (anyExceptions()) return
#line 1337 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of saturation as '//&
         &'a function of capillary pressure'
    value = 4.7607337955577978D-5
#line 1340 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1340) )
  if (anyExceptions()) return
#line 1341 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of relative perm as'// &
         &' a function of capillary pressure'
    value = 1.2551629198938072D-4
#line 1344 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1344) )
  if (anyExceptions()) return
#line 1345 "test_characteristic_curves.pf"

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
#line 1373 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1373) )
  if (anyExceptions()) return
#line 1374 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of capillary ' // &
         'pressure as a function of saturation'
    value = -24500.406539923177d0
#line 1377 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1377) )
  if (anyExceptions()) return
#line 1378 "test_characteristic_curves.pf"

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
#line 1405 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1405) )
  if (anyExceptions()) return
#line 1406 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of liquid relative'//&
         &' permeability as a function of liquid saturation'
    value = 0.92477440290277491d0
#line 1409 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1409) )
  if (anyExceptions()) return
#line 1410 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_mk3%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = '3-parameter modified Kosugi gas relative permeability as'//&
         &' a function of liquid saturation'
    value = 0.35040465818058608d0
#line 1418 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1418) )
  if (anyExceptions()) return
#line 1419 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of gas relative'//&
         &' permeability as a function of liquid saturation'
    value = -1.2770282015052450d0
#line 1422 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1422) )
  if (anyExceptions()) return
#line 1423 "test_characteristic_curves.pf"

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
#line 1451 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1451) )
  if (anyExceptions()) return
#line 1452 "test_characteristic_curves.pf"
    string = 'r-parameter modified Kosugi derivative of capillary ' // &
         'pressure as a function of saturation'
    value = -19192.362637570313d0
#line 1455 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1455) )
  if (anyExceptions()) return
#line 1456 "test_characteristic_curves.pf"

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
#line 1483 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1483) )
  if (anyExceptions()) return
#line 1484 "test_characteristic_curves.pf"
    string = '4-parameter modified Kosugi derivative of liquid relative'//&
         &' permeability as a function of liquid saturation'
    value = 0.92477440290277491d0  ! same as 3-param
#line 1487 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1487) )
  if (anyExceptions()) return
#line 1488 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_mk4%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = '4-parameter modified Kosugi gas relative permeability as'//&
         &' a function of liquid saturation'
    value = 0.35040465818058608d0  ! same as 3-param
#line 1496 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1496) )
  if (anyExceptions()) return
#line 1497 "test_characteristic_curves.pf"
    string = '4-parameter modified Kosugi derivative of gas relative'//&
         &' permeability as a function of liquid saturation'
    value = -1.2770282015052450d0   ! same as 3-param
#line 1500 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1500) )
  if (anyExceptions()) return
#line 1501 "test_characteristic_curves.pf"

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
#line 1528 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1528) )
  if (anyExceptions()) return
#line 1529 "test_characteristic_curves.pf"
    string = 'Constant saturation function liquid saturation derivative'
    value = 0.d0
#line 1531 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1531) )
  if (anyExceptions()) return
#line 1532 "test_characteristic_curves.pf"

    ! capillary_pressure = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_const%saturation_function% &
                                    CapillaryPressure(liquid_saturation, &
                                                      capillary_pressure, &
                                                      dpc_dsatl, &
                                                      this%option)
    string = 'Constant saturation function capillary pressure'
    value = 1.d5
#line 1542 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1542) )
  if (anyExceptions()) return
#line 1543 "test_characteristic_curves.pf"

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
#line 1572 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1572) )
  if (anyExceptions()) return
#line 1573 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp1%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP1 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.74556316913866438d0
#line 1580 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1580) )
  if (anyExceptions()) return
#line 1581 "test_characteristic_curves.pf"
    call this%cc_krp1%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.10878895018640269d0
#line 1587 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1587) )
  if (anyExceptions()) return
#line 1588 "test_characteristic_curves.pf"
    call this%cc_krp1%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 5.9765155008000600d-4
#line 1594 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1594) )
  if (anyExceptions()) return
#line 1595 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp1%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP1 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 0.0d0
#line 1603 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1603) )
  if (anyExceptions()) return
#line 1604 "test_characteristic_curves.pf"
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
#line 1614 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1614) )
  if (anyExceptions()) return
#line 1615 "test_characteristic_curves.pf"
    call this%cc_krp1%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 3.1269006686492169d-22
#line 1621 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1621) )
  if (anyExceptions()) return
#line 1622 "test_characteristic_curves.pf"
    call this%cc_krp1%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.0d0
#line 1628 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1628) )
  if (anyExceptions()) return
#line 1629 "test_characteristic_curves.pf"

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
#line 1658 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1658) )
  if (anyExceptions()) return
#line 1659 "test_characteristic_curves.pf"
    capillary_pressure = 100.d0
    call this%cc_krp2%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP2 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 1666 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1666) )
  if (anyExceptions()) return
#line 1667 "test_characteristic_curves.pf"
    call this%cc_krp2%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1673 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1673) )
  if (anyExceptions()) return
#line 1674 "test_characteristic_curves.pf"
    call this%cc_krp2%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1680 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1680) )
  if (anyExceptions()) return
#line 1681 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp2%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP2 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 51637.669399930157d0
#line 1689 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1689) )
  if (anyExceptions()) return
#line 1690 "test_characteristic_curves.pf"
    select type(sf=>this%cc_krp2%saturation_function)
      class is(sat_func_KRP2_type)
        capillary_pressure = 15.d0/sf%alpha
    end select
    call this%cc_krp2%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP2 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.55731947333915333d0
#line 1699 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1699) )
  if (anyExceptions()) return
#line 1700 "test_characteristic_curves.pf"
    call this%cc_krp2%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 5.8310805074531403d-4
#line 1706 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1706) )
  if (anyExceptions()) return
#line 1707 "test_characteristic_curves.pf"
    call this%cc_krp2%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.24138701885980715d0
#line 1713 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1713) )
  if (anyExceptions()) return
#line 1714 "test_characteristic_curves.pf"

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
#line 1743 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1743) )
  if (anyExceptions()) return
#line 1744 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp3%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP3 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 1750 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1750) )
  if (anyExceptions()) return
#line 1751 "test_characteristic_curves.pf"
    call this%cc_krp3%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1757 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1757) )
  if (anyExceptions()) return
#line 1758 "test_characteristic_curves.pf"
    call this%cc_krp3%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1764 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1764) )
  if (anyExceptions()) return
#line 1765 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp3%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP3 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 200.0d0
#line 1773 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1773) )
  if (anyExceptions()) return
#line 1774 "test_characteristic_curves.pf"
    capillary_pressure = 1.d5
    call this%cc_krp3%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP3 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.21815720128839766d0
#line 1780 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1780) )
  if (anyExceptions()) return
#line 1781 "test_characteristic_curves.pf"
    call this%cc_krp3%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 5.8632915084341688d-9
#line 1787 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1787) )
  if (anyExceptions()) return
#line 1788 "test_characteristic_curves.pf"
    call this%cc_krp3%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.78571287226820719d0
#line 1794 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1794) )
  if (anyExceptions()) return
#line 1795 "test_characteristic_curves.pf"

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
#line 1824 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1824) )
  if (anyExceptions()) return
#line 1825 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp4%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP4 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 1831 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1831) )
  if (anyExceptions()) return
#line 1832 "test_characteristic_curves.pf"
    call this%cc_krp4%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1838 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1838) )
  if (anyExceptions()) return
#line 1839 "test_characteristic_curves.pf"
    call this%cc_krp4%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1845 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1845) )
  if (anyExceptions()) return
#line 1846 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp4%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP4 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 152.32289193581161d0
#line 1854 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1854) )
  if (anyExceptions()) return
#line 1855 "test_characteristic_curves.pf"
    capillary_pressure = 1.d5
    call this%cc_krp4%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP4 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.21815720128839766d0
#line 1861 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1861) )
  if (anyExceptions()) return
#line 1862 "test_characteristic_curves.pf"
    call this%cc_krp4%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 2.8180697920854255d-10
#line 1868 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1868) )
  if (anyExceptions()) return
#line 1869 "test_characteristic_curves.pf"
    call this%cc_krp4%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.78571287226820719d0
#line 1875 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1875) )
  if (anyExceptions()) return
#line 1876 "test_characteristic_curves.pf"

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
#line 1905 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1905) )
  if (anyExceptions()) return
#line 1906 "test_characteristic_curves.pf"
    capillary_pressure = 50.d0
    call this%cc_krp5%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP5 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.80014642571085304d0
#line 1913 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1913) )
  if (anyExceptions()) return
#line 1914 "test_characteristic_curves.pf"
    call this%cc_krp5%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1920 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1920) )
  if (anyExceptions()) return
#line 1921 "test_characteristic_curves.pf"
    call this%cc_krp5%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1927 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1927) )
  if (anyExceptions()) return
#line 1928 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp5%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP5 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 2000.0d0
#line 1936 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1936) )
  if (anyExceptions()) return
#line 1937 "test_characteristic_curves.pf"
    capillary_pressure = 8.d5
    call this%cc_krp5%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP5 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.74007809371245503d0
#line 1944 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1944) )
  if (anyExceptions()) return
#line 1945 "test_characteristic_curves.pf"
    call this%cc_krp5%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.92010412494993998d0
#line 1951 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1951) )
  if (anyExceptions()) return
#line 1952 "test_characteristic_curves.pf"
    call this%cc_krp5%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 7.9895875050060017d-2
#line 1958 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1958) )
  if (anyExceptions()) return
#line 1959 "test_characteristic_curves.pf"

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
#line 1988 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1988) )
  if (anyExceptions()) return
#line 1989 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp8%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP8 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.99813934029313034d0
#line 1996 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1996) )
  if (anyExceptions()) return
#line 1997 "test_characteristic_curves.pf"
    call this%cc_krp8%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.82967396002016347d0
#line 2003 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2003) )
  if (anyExceptions()) return
#line 2004 "test_characteristic_curves.pf"
    call this%cc_krp8%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 3.5744473041460765d-4
#line 2010 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2010) )
  if (anyExceptions()) return
#line 2011 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp8%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP8 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 11673.176741409796d0
#line 2019 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2019) )
  if (anyExceptions()) return
#line 2020 "test_characteristic_curves.pf"
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
#line 2030 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2030) )
  if (anyExceptions()) return
#line 2031 "test_characteristic_curves.pf"
    call this%cc_krp8%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 2.4705248145444649d-14
#line 2037 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2037) )
  if (anyExceptions()) return
#line 2038 "test_characteristic_curves.pf"
    call this%cc_krp8%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.99877541173464723d0
#line 2044 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2044) )
  if (anyExceptions()) return
#line 2045 "test_characteristic_curves.pf"

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
#line 2074 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2074) )
  if (anyExceptions()) return
#line 2075 "test_characteristic_curves.pf"
    capillary_pressure = 50.d0
    call this%cc_krp9%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP9 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.99999644138299326d0
#line 2082 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2082) )
  if (anyExceptions()) return
#line 2083 "test_characteristic_curves.pf"
    call this%cc_krp9%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.99999998839687010d0
#line 2089 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2089) )
  if (anyExceptions()) return
#line 2090 "test_characteristic_curves.pf"
    call this%cc_krp9%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.1603129901338605d-8
#line 2096 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2096) )
  if (anyExceptions()) return
#line 2097 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp9%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP9 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 2294.5062171637132d0
#line 2105 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2105) )
  if (anyExceptions()) return
#line 2106 "test_characteristic_curves.pf"
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
#line 2116 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2116) )
  if (anyExceptions()) return
#line 2117 "test_characteristic_curves.pf"
    call this%cc_krp9%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.d0
#line 2123 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2123) )
  if (anyExceptions()) return
#line 2124 "test_characteristic_curves.pf"
    call this%cc_krp9%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.d0
#line 2130 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2130) )
  if (anyExceptions()) return
#line 2131 "test_characteristic_curves.pf"

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
#line 2160 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2160) )
  if (anyExceptions()) return
#line 2161 "test_characteristic_curves.pf"
    capillary_pressure = 100.d0
    call this%cc_krp11%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP11 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 2168 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2168) )
  if (anyExceptions()) return
#line 2169 "test_characteristic_curves.pf"
    call this%cc_krp11%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 2175 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2175) )
  if (anyExceptions()) return
#line 2176 "test_characteristic_curves.pf"
    call this%cc_krp11%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.d0
#line 2182 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2182) )
  if (anyExceptions()) return
#line 2183 "test_characteristic_curves.pf"

    liquid_saturation = 0.85
    call this%cc_krp11%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP11 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 0.0d0
#line 2191 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2191) )
  if (anyExceptions()) return
#line 2192 "test_characteristic_curves.pf"
    capillary_pressure = 1.d6
    call this%cc_krp11%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP11 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 1.0d0
#line 2199 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2199) )
  if (anyExceptions()) return
#line 2200 "test_characteristic_curves.pf"
    call this%cc_krp11%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.0d0
#line 2206 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2206) )
  if (anyExceptions()) return
#line 2207 "test_characteristic_curves.pf"
    call this%cc_krp11%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.0d0
#line 2213 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2213) )
  if (anyExceptions()) return
#line 2214 "test_characteristic_curves.pf"

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
#line 2243 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2243) )
  if (anyExceptions()) return
#line 2244 "test_characteristic_curves.pf"
    capillary_pressure = 100.d0
    call this%cc_krp12%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP12 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 2251 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2251) )
  if (anyExceptions()) return
#line 2252 "test_characteristic_curves.pf"
    call this%cc_krp12%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 2258 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2258) )
  if (anyExceptions()) return
#line 2259 "test_characteristic_curves.pf"
    call this%cc_krp12%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.d0
#line 2265 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2265) )
  if (anyExceptions()) return
#line 2266 "test_characteristic_curves.pf"

    liquid_saturation = 0.85
    call this%cc_krp12%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP12 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 4547.8850244305122d0
#line 2274 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2274) )
  if (anyExceptions()) return
#line 2275 "test_characteristic_curves.pf"
    capillary_pressure = 1.d6
    call this%cc_krp12%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP12 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.30276918155781385d0
#line 2282 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2282) )
  if (anyExceptions()) return
#line 2283 "test_characteristic_curves.pf"
    call this%cc_krp12%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.3380237859265290d-9
#line 2289 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2289) )
  if (anyExceptions()) return
#line 2290 "test_characteristic_curves.pf"
    call this%cc_krp12%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.53468502607262869d0
#line 2296 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2296) )
  if (anyExceptions()) return
#line 2297 "test_characteristic_curves.pf"

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

