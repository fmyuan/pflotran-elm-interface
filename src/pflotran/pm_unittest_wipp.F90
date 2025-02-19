module PM_Unit_Test_WIPP_class

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PM_Base_class
  use PM_Unit_Test_class
  use Realization_Subsurface_class
  use Communicator_Base_class
  use Material_module
  use Material_Aux_module
  use Option_module

  use PFLOTRAN_Constants_module

  private

  public :: PMUnitTestWIPPRead, &
            PMUnitTestWIPPRunFracture, &
            PMUnitTestWIPPRunKlinkenberg, &
            PMUnitTestWIPPRunPcSat, &
            PMUnitTestWIPPRunGasPermeability, &
            PMUnitTestWIPPRunLiqPermeability, &
            PMUnitTestWIPPRunCompressibility, &
            PMUnitTestWIPPRunEOSGasDensity, &
            PMUnitTestWIPPRunEOSWaterDensity
contains

! ************************************************************************** !

subroutine PMUnitTestWIPPRead(input,option,this)
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !

  use Input_Aux_module
  use Option_module
  use String_module

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_unittest_type), pointer :: this
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_str

  error_str = 'SIMULATION,PROCESS_MODELS,UNITTEST'
  input%ierr = INPUT_ERROR_NONE

  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_str)
    call StringToUpper(keyword)
    select case(trim(keyword))
      case('WIPP_FRACTURE')
        call InputReadWord(input,option,this%filename,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        this%RunUnitTest => PMUnitTestWIPPRunFracture
        call InputReadDouble(input,option,this%tolerance)
        allocate(this%next_unittest)
        this => this%next_unittest
      case('WIPP_KLINKENBERG')
        call InputReadWord(input,option,this%filename,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        this%RunUnitTest => PMUnitTestWIPPRunKlinkenberg
        call InputReadDouble(input,option,this%tolerance)
        allocate(this%next_unittest)
        this => this%next_unittest
      case('WIPP_COMPRESSIBILITY')
        call InputReadWord(input,option,this%filename,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        this%RunUnitTest => PMUnitTestWIPPRunCompressibility
        call InputReadDouble(input,option,this%tolerance)
        allocate(this%next_unittest)
        this => this%next_unittest
      case('EOS_GAS')
        call InputReadWord(input,option,this%filename,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        this%RunUnitTest => PMUnitTestWIPPRunEOSGasDensity
        call InputReadDouble(input,option,this%tolerance)
        allocate(this%next_unittest)
        this => this%next_unittest
      case('EOS_WATER')
        call InputReadWord(input,option,this%filename,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        this%RunUnitTest => PMUnitTestWIPPRunEOSWaterDensity
        call InputReadDouble(input,option,this%tolerance)
        allocate(this%next_unittest)
        this => this%next_unittest
      case('PC_SAT')
        call InputReadWord(input,option,this%filename,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        this%RunUnitTest => PMUnitTestWIPPRunPcSat
        call InputReadDouble(input,option,this%tolerance)
        allocate(this%next_unittest)
        this => this%next_unittest
      case('RELPERM_GAS')
        call InputReadWord(input,option,this%filename,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        this%RunUnitTest => PMUnitTestWIPPRunGasPermeability
        call InputReadDouble(input,option,this%tolerance)
        allocate(this%next_unittest)
        this => this%next_unittest
      case('RELPERM_LIQUID')
        call InputReadWord(input,option,this%filename,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        this%RunUnitTest => PMUnitTestWIPPRunLiqPermeability
        call InputReadDouble(input,option,this%tolerance)
        allocate(this%next_unittest)
        this => this%next_unittest
      case default
        call InputKeywordUnrecognized(input,keyword,error_str,option)
    end select
  enddo

  call InputPopBlock(input,option)

end subroutine PMUnitTestWIPPRead

! ************************************************************************** !

subroutine PMUnitTestWIPPRunFracture(this,realization)
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !
  use Fracture_module

  implicit none

  type(material_property_type), pointer :: cur_mat
  class(pm_unittest_type) :: this
  type(material_auxvar_type), pointer :: material_aux
  type(realization_subsurface_type) :: realization

  ! Fracture values
  PetscInt, allocatable  :: material_id(:)
  PetscInt               :: temp_material_id(150)
  PetscReal, allocatable :: init_liq_pressure(:)
  PetscReal              :: temp_init_liq_pressure(150)
  PetscReal, allocatable :: liq_pressure(:)
  PetscReal              :: temp_liq_pressure(150)
  PetscReal, allocatable :: correct_compressed_phi(:)
  PetscReal              :: temp_correct_compressed_phi(150)
  PetscReal, allocatable :: correct_scaling_factor(:)
  PetscReal              :: temp_correct_scaling_factor(150)
  PetscReal, allocatable :: compressed_porosity(:)
  PetscReal              :: dcompressed_porosity_dp
  PetscReal, allocatable :: scaling_factor(:)
  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscInt :: values(8)
  PetscReal :: tolerance
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time
  character(len=MAXWORDLENGTH) :: pass_fail, filename_out, input_filename

  PetscInt  :: i, j, k

  input_filename = trim(this%filename)
  ! Read input file
  realization%option%io_buffer = 'Running fracture unit test: '//trim(this%filename)
  call PrintMsg(realization%option)
  i = 1
  open(action='read', file=trim(input_filename), iostat=rc_in, &
       newunit=fu_in)
  if (rc_in /= 0) then
    realization%option%io_buffer = 'File read error: fracture unit test: '//trim(this%filename)
    call PrintErrMsg(realization%option)
  endif
  read(fu_in, *) ! skip header line
  do
    read (fu_in, *, iostat=rc_in) temp_material_id(i), &
                                  temp_init_liq_pressure(i), &
                                  temp_liq_pressure(i), &
                                  temp_correct_compressed_phi(i), &
                                  temp_correct_scaling_factor(i)
    if (rc_in /= 0) exit
    i = i + 1
    if (i > 150) exit
  enddo

  allocate(material_id(i-1))
  allocate(init_liq_pressure(i-1))
  allocate(liq_pressure(i-1))
  allocate(correct_compressed_phi(i-1))
  allocate(correct_scaling_factor(i-1))

  material_id(:) = temp_material_id(1:i-1)
  init_liq_pressure(:)      = temp_init_liq_pressure(1:i-1)
  liq_pressure(:)           = temp_liq_pressure(1:i-1)
  correct_compressed_phi(:) = temp_correct_compressed_phi(1:i-1)
  correct_scaling_factor(:) = temp_correct_scaling_factor(1:i-1)

  allocate(compressed_porosity(i-1))
  allocate(scaling_factor(i-1))

  ! Fracture unit test
  !   ! inputs:
  !   !   From PFLOTRAN:       Ci, P0, phi0, Pi, Pa, Phia, n
  !   !   From Unit test file:
  !   ! outputs:
  !   !   Altered porosity
  !   !   Permeability scaling factor
  !   ! functions tested:
  !   !   FracturePoroEvaluate :: compressed_porosity
  !   !   FracturePermScale    :: scaling_factor
  !   ! loop through each material, read input parameters from file, output
  !   ! material id, input parameters, altered porosity and permeability scaling factor

  cur_mat => realization%material_properties
  ! material_properties is a linked list of material properties.
  ! material_properties have an external (pflotran input deck) and internal id.
  ! The internal id associated with the external material id  is not known prior
  ! to traversing the material_properties linked list, so each material_property must
  ! be looped through, and the below functions are run for each property.
  do i = 1, size(material_id)
    if (cur_mat%external_id /= material_id(i)) then
      cur_mat => cur_mat%next
      if (.not.associated(cur_mat)) then
        cur_mat => realization%material_properties !loop back through to find the right material
      endif
    endif
    ! here, loop through pressures read in from the input deck
    do j = 1, size(realization%patch%aux%material%auxvars)
      if (realization%patch%aux%material%auxvars(j)%id == cur_mat%internal_id) then
        material_aux => realization%patch%aux%material%auxvars(j)
      endif
    end do

    material_aux%fracture%initial_pressure = init_liq_pressure(i)
    call FracturePoroEvaluate(material_aux,liq_pressure(i),&
                              compressed_porosity(i),dcompressed_porosity_dp)
    call FracturePermScale(material_aux,liq_pressure(i),&
                           compressed_porosity(i),scaling_factor(i))
 enddo
  tolerance = this%tolerance

  call PMUnitTestWIPPFilenameInToOut(input_filename, filename_out)
  filename_out = trim(filename_out)
  open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
  call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
  write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                  time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
  write(fu_out,*)
  write(fu_out,'(a)') 'NOTE: The input file provided was:'
  write(fu_out,'(a,a)') '      ', trim(input_filename)
  write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                               tolerance, '.'
  write(fu_out,*)
  i = 0
  do k=1,size(init_liq_pressure)
    write(fu_out,'(a,I3,a)') '||-----------TEST-#',k,'-----------------------&
                             &--------------||'
    write(fu_out,'(a)') '  [in]  init. liquid pressure [Pa]:'
    write(fu_out,'(d17.10)') init_liq_pressure(k)
    write(fu_out,'(a)') '  [in]  liquid pressure [Pa]:'
    write(fu_out,'(d17.10)') liq_pressure(k)

    write(fu_out,'(a)') '  [out]  altered porosity [-]:'
    write(fu_out,'(d17.10)') compressed_porosity(k)
    write(fu_out,'(a)') '  [correct]  altered porosity [-]:'
    write(fu_out,'(d17.10)') correct_compressed_phi(k)

    call CalcDiff(compressed_porosity(k),correct_compressed_phi(k),tolerance,pass_fail,i)
    write(fu_out,'(a)') trim(pass_fail)

    write(fu_out,'(a)') '  [out]  perm scaling factor [-]:'
    write(fu_out,'(d17.10)') scaling_factor(k)
    write(fu_out,'(a)') '  [correct]  perm scaling factor [-]:'
    write(fu_out,'(d17.10)') correct_scaling_factor(k)
    call CalcDiff(scaling_factor(k),correct_scaling_factor(k),tolerance,pass_fail,i)
    write(fu_out,'(a)') trim(pass_fail)
    write(fu_out,*)
  enddo

  write(fu_out,'(a)') 'TEST SUMMARY:'
  if (i == 0) then
    write(fu_out,'(a)') ' All tests passed!'
  else
    write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
  endif

  close(fu_out)
  close(fu_in)

end subroutine PMUnitTestWIPPRunFracture

! ************************************************************************** !

subroutine PMUnitTestWIPPRunKlinkenberg(this,realization)
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !

  use Klinkenberg_module

  implicit none

  class(pm_unittest_type) :: this
  type(realization_subsurface_type) :: realization

  PetscReal, allocatable :: perm_x(:)
  PetscReal, allocatable :: perm_y(:)
  PetscReal, allocatable :: perm_z(:)
  PetscReal, allocatable :: gas_pressure(:)
  PetscReal, allocatable :: correct_perm_out_x(:)
  PetscReal, allocatable :: correct_perm_out_y(:)
  PetscReal, allocatable :: correct_perm_out_z(:)
  PetscReal :: temp_perm_x(150)
  PetscReal :: temp_perm_y(150)
  PetscReal :: temp_perm_z(150)
  PetscReal :: temp_gas_pressure(150)
  PetscReal :: temp_correct_perm_out_x(150)
  PetscReal :: temp_correct_perm_out_y(150)
  PetscReal :: temp_correct_perm_out_z(150)

  PetscReal, allocatable :: permeability_scale(:,:)
  PetscReal, allocatable :: scaling_factor(:)

  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscInt :: values(8)
  PetscReal :: tolerance
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time
  character(len=MAXWORDLENGTH) :: pass_fail_x,pass_fail_y,pass_fail_z, filename_out, input_filename
  PetscInt  :: i,k

  klinkenberg => KlinkenbergGetPtr()

  input_filename = trim(this%filename)
  ! Read input file
  realization%option%io_buffer = 'Running Klinkenberg unit test, input ' // trim(this%filename)
  call PrintMsg(realization%option)
  i = 1
  open(action='read', file=trim(input_filename), iostat=rc_in, &
       newunit=fu_in)
  if (rc_in /= 0) then
    realization%option%io_buffer = 'File read error: Klinkenberg unit test: '//trim(this%filename)
    call PrintErrMsg(realization%option)
  endif
  read(fu_in, *) ! skip header line
  do
    read (fu_in, *, iostat=rc_in) temp_perm_x(i), &
                                  temp_perm_y(i), &
                                  temp_perm_z(i), &
                                  temp_gas_pressure(i), &
                                  temp_correct_perm_out_x(i), &
                                  temp_correct_perm_out_y(i), &
                                  temp_correct_perm_out_z(i)
    if (rc_in /= 0) exit
    i = i + 1
    if (i > 150) exit
  enddo


  allocate(perm_x(i-1))
  allocate(perm_y(i-1))
  allocate(perm_z(i-1))
  allocate(gas_pressure(i-1))
  allocate(correct_perm_out_x(i-1))
  allocate(correct_perm_out_y(i-1))
  allocate(correct_perm_out_z(i-1))
  allocate(permeability_scale(i-1,3))

  perm_x(:) = temp_perm_x(1:i-1)
  perm_y(:) = temp_perm_y(1:i-1)
  perm_z(:) = temp_perm_z(1:i-1)
  gas_pressure(:) = temp_gas_pressure(1:i-1)
  correct_perm_out_x(:) = temp_correct_perm_out_x(1:i-1)
  correct_perm_out_y(:) = temp_correct_perm_out_y(1:i-1)
  correct_perm_out_z(:) = temp_correct_perm_out_z(1:i-1)

  allocate(scaling_factor(i-1))

  do i = 1, size(perm_x)
    call klinkenberg%Scale((/perm_x(i),perm_y(i),perm_z(i)/),gas_pressure(i),permeability_scale(i,:))
  end do

  tolerance = this%tolerance

  call PMUnitTestWIPPFilenameInToOut(input_filename, filename_out)
  filename_out = trim(filename_out)
  open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
  call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
  write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                  time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
  write(fu_out,*)
  write(fu_out,'(a)') 'NOTE: The input file provided was:'
  write(fu_out,'(a,a)') '      ', trim(input_filename)
  write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                               tolerance, '.'
  write(fu_out,*)
  i = 0
  do k=1,size(perm_x)
    write(fu_out,'(a,I3,a)') '||-----------TEST-#',k,'-----------------------&
                             &--------------||'
    write(fu_out,'(a)') '[in]  liquid permeability (x,y,z) [m2]:'
    write(fu_out,'(" ("d17.10,", "d17.10,", "d17.10")")') (/perm_x(k),perm_y(k),perm_z(k)/)
    write(fu_out,'(a)') '[in]  gas pressure [Pa]:'
    write(fu_out,'(" ("d17.10,", "d17.10,", "d17.10")")') gas_pressure(k)

    write(fu_out,'(a)') '[out]  gas permeability (x,y,z) [m2]:'
    write(fu_out,'(" ("d17.10,", "d17.10,", "d17.10")")') (/perm_x(k)*permeability_scale(k,1),perm_y(k)*permeability_scale(k,2),perm_z(k)*permeability_scale(k,3)/)
    write(fu_out,'(a)') '[correct]  gas permeability (x,y,z) [m2]:'
    write(fu_out,'(" ("d17.10,", "d17.10,", "d17.10")")') (/correct_perm_out_x(k),correct_perm_out_y(k),correct_perm_out_z(k)/)

    call CalcDiff(perm_x(k)*permeability_scale(k,1),correct_perm_out_x(k),tolerance,pass_fail_x,i)
    call CalcDiff(perm_y(k)*permeability_scale(k,2),correct_perm_out_x(k),tolerance,pass_fail_y,i)
    call CalcDiff(perm_z(k)*permeability_scale(k,3),correct_perm_out_x(k),tolerance,pass_fail_z,i)
    write(fu_out,'(a,", "a,", "a)') trim(pass_fail_x), trim(pass_fail_y), trim(pass_fail_z)
    write(fu_out,*)
  enddo

  write(fu_out,'(a)') 'TEST SUMMARY:'
  if (i == 0) then
    write(fu_out,'(a)') ' All tests passed!'
  else
    write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
  endif

  close(fu_out)
  close(fu_in)

end subroutine PMUnitTestWIPPRunKlinkenberg

! ************************************************************************** !

subroutine PMUnitTestWIPPRunPcSat(this,realization)
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !

  use Characteristic_curves_module
  use Option_module
  use String_module

  implicit none

  class(pm_unittest_type) :: this
  type(characteristic_curves_type), pointer :: cc
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH), allocatable :: cc_name(:)
  character(len=MAXWORDLENGTH) :: temp_cc_name(150)
  PetscReal, allocatable :: pc(:)
  PetscReal :: temp_pc(150)
  PetscReal, allocatable :: correct_sat(:)
  PetscReal :: temp_correct_sat(150)
  PetscReal, allocatable :: liq_sat(:)
  PetscReal :: temp_liq_sat(150)
  PetscReal, allocatable :: correct_pc(:)
  PetscReal :: temp_correct_pc(150)

  PetscReal, allocatable :: calc_sat(:)
  PetscReal, allocatable :: calc_pc(:)
  PetscReal :: dummy

  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscInt :: values(8)
  PetscReal :: tolerance
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time
  character(len=MAXWORDLENGTH) :: pass_fail, filename_out, input_filename
  PetscBool :: cycled = PETSC_FALSE

  PetscInt  :: i,k

  input_filename = trim(this%filename)
  ! Read input file
  realization%option%io_buffer = 'Running capillary pressure-saturation unit test, input '//trim(this%filename)
  call PrintMsg(realization%option)
  i = 1
  open(action='read', file=trim(input_filename), iostat=rc_in, &
       newunit=fu_in)
  if (rc_in /= 0) then
    realization%option%io_buffer = 'File read error: capillary pressure-sauration unit test: '//trim(this%filename)
    call PrintErrMsg(realization%option)
  endif
  read(fu_in, *) ! skip header line
  do
    read (fu_in, *, iostat=rc_in) temp_cc_name(i), &
                                  temp_pc(i), &
                                  temp_correct_sat(i), &
                                  temp_liq_sat(i), &
                                  temp_correct_pc(i)
    if (rc_in /= 0) exit
    i = i + 1
    if (i > 150) exit
  enddo

  allocate(cc_name(i-1))
  allocate(pc(i-1))
  allocate(correct_sat(i-1))
  allocate(liq_sat(i-1))
  allocate(correct_pc(i-1))

  cc_name(:)     = temp_cc_name(1:i-1)
  pc(:)          = temp_pc(1:i-1)
  correct_sat(:) = temp_correct_sat(1:i-1)
  liq_sat(:)     = temp_liq_sat(1:i-1)
  correct_pc(:)  = temp_correct_pc(1:i-1)

  allocate(calc_sat(i-1))
  allocate(calc_pc(i-1))

  option => this%option
  cc => realization%characteristic_curves

  do i = 1, size(calc_pc)
     do while (.not.StringCompare(cc%name,trim(cc_name(i))))
       cc => cc%next
       if (.not.associated(cc)) then
         cc => realization%characteristic_curves
         if (cycled) then
           option%io_buffer = 'Characteristic curve not found.'
           call PrintErrMsg(option)
         endif
         cycled = PETSC_TRUE
       endif
     end do
     call cc%saturation_function%Saturation(pc(i),calc_sat(i),dummy,option)
     call cc%saturation_function%CapillaryPressure(liq_sat(i),calc_pc(i),dummy,option)
     cycled = PETSC_FALSE
  end do

  tolerance = this%tolerance

  call PMUnitTestWIPPFilenameInToOut(input_filename, filename_out)
  filename_out = trim(filename_out)
  open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
  call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
  write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                  time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
  write(fu_out,*)
  write(fu_out,'(a)') 'NOTE: The input file provided was:'
  write(fu_out,'(a,a)') '      ', trim(input_filename)
  write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                               tolerance, '.'
  write(fu_out,*)
  i = 0
  do k=1,size(pc)
    write(fu_out,'(a,I3,a)') '||-----------TEST-#',k,'-----------------------&
                             &--------------||'
    write(fu_out,'(a)') '  [in]  characteristic curves name:'
    write(fu_out,'(a)') trim(cc_name(k))
    write(fu_out,'(a)') '  [in]  capillary pressure [Pa]:'
    write(fu_out,'(d17.10)') pc(k)

    write(fu_out,'(a)') '  [out]  liquid saturation [-]:'
    write(fu_out,'(d17.10)') calc_sat(k)
    write(fu_out,'(a)') '  [correct]  liquid saturation [-]:'
    write(fu_out,'(d17.10)') correct_sat(k)

    call CalcDiff(calc_sat(k),correct_sat(k),tolerance,pass_fail,i)
    write(fu_out,'(a)') trim(pass_fail)

    write(fu_out,'(a)') '  [in]  liquid saturation [-]:'
    write(fu_out,'(d17.10)') liq_sat(k)
    write(fu_out,'(a)') '  [out]  capillary pressure [Pa]:'
    write(fu_out,'(d17.10)') calc_pc(k)
    write(fu_out,'(a)') '  [correct]  capillary pressure [Pa]:'
    write(fu_out,'(d17.10)') correct_pc(k)
    call CalcDiff(calc_pc(k),correct_pc(k),tolerance,pass_fail,i)
    write(fu_out,'(a)') trim(pass_fail)
    write(fu_out,*)
  enddo

  write(fu_out,'(a)') 'TEST SUMMARY:'
  if (i == 0) then
    write(fu_out,'(a)') ' All tests passed!'
  else
    write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
  endif

  close(fu_out)
  close(fu_in)

!  realization%characteristic_curves%name
  !call WIPPCharacteristicCurves()

end subroutine PMUnitTestWIPPRunPcSat

! ************************************************************************** !

subroutine PMUnitTestWIPPRunGasPermeability(this,realization)
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !

  use WIPP_Characteristic_Curve_module
  use Characteristic_curves_module
  use Option_module
  use String_module

  implicit none

  class(pm_unittest_type) :: this
  type(characteristic_curves_type), pointer :: cc
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: temp_cc_name(150)
  character(len=MAXWORDLENGTH), allocatable :: cc_name(:)
  PetscReal, allocatable :: gas_saturation(:)
  PetscReal :: temp_gas_saturation(150)
  PetscReal, allocatable :: correct_relperm(:)
  PetscReal :: temp_correct_relperm(150)
  PetscReal :: temp_krg(150)
  PetscReal, allocatable :: krg(:)
  PetscBool :: cycled = PETSC_TRUE
  PetscReal :: dummy
  PetscReal, parameter :: perm = 1.d-12
  PetscInt :: i,k
  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscInt :: values(8)
  PetscReal :: tolerance
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time
  character(len=MAXWORDLENGTH) :: pass_fail, filename_out, input_filename


  input_filename = this%filename
  ! Read input file
  realization%option%io_buffer = 'Running gas relative permeability unit test, input '//trim(this%filename)
  call PrintMsg(realization%option)
  i = 1
  open(action='read', file=trim(input_filename), iostat=rc_in, &
       newunit=fu_in)
  if (rc_in /= 0) then
    realization%option%io_buffer = 'File read error: gas relative permeability unit test: '//trim(this%filename)
    call PrintErrMsg(realization%option)
  endif
  read(fu_in, *) ! skip header line
  do
    read (fu_in, *, iostat=rc_in) temp_cc_name(i), &
                                  temp_gas_saturation(i), &
                                  temp_correct_relperm(i)
    if (rc_in /= 0) exit
    i = i + 1
    if (i > 150) exit
  enddo

  allocate(cc_name(i-1))
  allocate(gas_saturation(i-1))
  allocate(correct_relperm(i-1))
  allocate(krg(i-1))

  cc_name(:) = temp_cc_name(1:i-1)
  gas_saturation(:) = temp_gas_saturation(1:i-1)
  correct_relperm(:) = temp_correct_relperm(1:i-1)
  krg(:) = temp_krg(1:i-1)

  option => this%option
  cc => realization%characteristic_curves

  do i = 1, size(correct_relperm)
     do while (.not.StringCompare(cc%name,trim(cc_name(i))))
       cc => cc%next
       if (.not.associated(cc)) then
         cc => realization%characteristic_curves ! loop back through cc linked list
         if (cycled) then
           option%io_buffer = 'Characteristic curve not found.'
           call PrintErrMsg(option)
         endif
         cycled = PETSC_TRUE
       endif
     end do
     call cc%gas_rel_perm_function%RelativePermeability(gas_saturation(i),krg(i),dummy,option)
     cycled = PETSC_FALSE
  end do

  tolerance = this%tolerance

  call PMUnitTestWIPPFilenameInToOut(input_filename, filename_out)
  filename_out = trim(filename_out)
  open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
  call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
  write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                  time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
  write(fu_out,*)
  write(fu_out,'(a)') 'NOTE: The input file provided was:'
  write(fu_out,'(a,a)') '      ', trim(input_filename)
  write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                               tolerance, '.'
  write(fu_out,*)
  i = 0
  do k=1,size(krg)
    write(fu_out,'(a,I3,a)') '||-----------TEST-#',k,'-----------------------&
                             &--------------||'
    write(fu_out,'(a)') '  [in]  characteristic curves name:'
    write(fu_out,'(a)') trim(cc_name(k))
    write(fu_out,'(a)') '  [in]  saturation (gas) [Pa]:'
    write(fu_out,'(d17.10)') gas_saturation(k)

    write(fu_out,'(a)') '  [out]  relative permeability (gas) [m2]'
    write(fu_out,'(d17.10)') krg(k)
    write(fu_out,'(a)') '  [correct]  relative permeability (gas) [m2]'
    write(fu_out,'(d17.10)') correct_relperm(k)

    call CalcDiff(krg(k),correct_relperm(k),tolerance,pass_fail,i)
    write(fu_out,'(a)') trim(pass_fail)
    write(fu_out,*)
  enddo

  write(fu_out,'(a)') 'TEST SUMMARY:'
  if (i == 0) then
    write(fu_out,'(a)') ' All tests passed!'
  else
    write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
  endif

  close(fu_out)
  close(fu_in)

end subroutine PMUnitTestWIPPRunGasPermeability

! ************************************************************************** !

subroutine PMUnitTestWIPPRunLiqPermeability(this,realization)
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !

  use WIPP_Characteristic_Curve_module
  use Characteristic_curves_module
  use Option_module
  use String_module

  implicit none

  class(pm_unittest_type) :: this
  type(characteristic_curves_type), pointer :: cc
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: temp_cc_name(150)
  character(len=MAXWORDLENGTH), allocatable :: cc_name(:)
  PetscReal :: temp_liq_saturation(150)
  PetscReal, allocatable :: liq_saturation(:)
  PetscReal :: temp_correct_relperm(150)
  PetscReal, allocatable :: correct_relperm(:)

  PetscReal, allocatable :: krl(:)
  PetscReal :: temp_krl(150)
  PetscReal :: dummy
  PetscReal, parameter :: perm = 1.d-12
  PetscBool :: cycled = PETSC_FALSE
  PetscInt :: i,k
  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscInt :: values(8)
  PetscReal :: tolerance
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time
  character(len=MAXWORDLENGTH) :: pass_fail, filename_out, input_filename

  option => this%option

  input_filename = trim(this%filename)
  ! Read input file
  i = 1
  realization%option%io_buffer = 'Running liquid relative permeability unit test, input '//trim(this%filename)
  call PrintMsg(realization%option)
  open(action='read', file=trim(input_filename), iostat=rc_in, &
       newunit=fu_in)
  read(fu_in, *) ! skip header line
  do
    read (fu_in, *, iostat=rc_in) temp_cc_name(i), &
                                  temp_liq_saturation(i), &
                                  temp_correct_relperm(i)
    if (rc_in /= 0) exit
    i = i + 1
    if (i > 150) exit
  enddo

  allocate(cc_name(i-1))
  allocate(liq_saturation(i-1))
  allocate(correct_relperm(i-1))
  allocate(krl(i-1))

  cc_name(:) = temp_cc_name(1:i-1)
  liq_saturation(:) = temp_liq_saturation(1:i-1)
  correct_relperm(:) = temp_correct_relperm(1:i-1)
  krl(:) = temp_krl(1:i-1)

  cc => realization%characteristic_curves
  do i = 1, size(correct_relperm)
    do while (.not.StringCompare(cc%name,trim(cc_name(i))))
      cc => cc%next
      if (.not.associated(cc)) then
        cc => realization%characteristic_curves
         if (cycled) then
           option%io_buffer = 'Characteristic curve not found.'
           call PrintErrMsg(option)
         endif
         cycled = PETSC_TRUE
      endif
     end do
     call cc%liq_rel_perm_function%RelativePermeability(liq_saturation(i),krl(i),dummy,option)
     cycled = PETSC_FALSE
  end do

  tolerance = this%tolerance

  call PMUnitTestWIPPFilenameInToOut(input_filename, filename_out)
  filename_out = trim(filename_out)
  open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
  call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
  write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                  time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
  write(fu_out,*)
  write(fu_out,'(a)') 'NOTE: The input file provided was:'
  write(fu_out,'(a,a)') '      ', trim(input_filename)
  write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                               tolerance, '.'
  write(fu_out,*)
  i = 0
  do k=1,size(krl)
    write(fu_out,'(a,I3,a)') '||-----------TEST-#',k,'-----------------------&
                             &--------------||'
    write(fu_out,'(a)') '  [in]  characteristic curves name:'
    write(fu_out,'(a)') trim(cc_name(k))
    write(fu_out,'(a)') '  [in]  saturation (liq) [Pa]:'
    write(fu_out,'(d17.10)') liq_saturation(k)

    write(fu_out,'(a)') '  [out]  relative permeability (gas) [m2]'
    write(fu_out,'(d17.10)') krl(k)
    write(fu_out,'(a)') '  [correct]  relative permeability (gas0 [m2]'
    write(fu_out,'(d17.10)') correct_relperm(k)

    call CalcDiff(krl(k),correct_relperm(k),tolerance,pass_fail,i)
    write(fu_out,'(a)') trim(pass_fail)
    write(fu_out,*)
  enddo

  write(fu_out,'(a)') 'TEST SUMMARY:'
  if (i == 0) then
    write(fu_out,'(a)') ' All tests passed!'
  else
    write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
  endif

  close(fu_out)
  close(fu_in)

end subroutine PMUnitTestWIPPRunLiqPermeability

! ************************************************************************** !

subroutine PMUnitTestWIPPRunCompressibility(this,realization)
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !

  use Fracture_module

  implicit none

  class(pm_unittest_type) :: this
  type(realization_subsurface_type) :: realization
  character(len=MAXWORDLENGTH), allocatable :: material_name(:)
  character(len=MAXWORDLENGTH) :: temp_material_name(150)
  PetscReal, allocatable :: material_id(:)
  PetscReal :: temp_material_id(150)
  PetscReal, allocatable :: liquid_pressure(:)
  PetscReal :: temp_liquid_pressure(150)
  PetscReal, allocatable :: porosity(:)
  PetscReal, allocatable :: correct_altered_porosity(:)
  PetscReal :: temp_correct_altered_porosity(150)
  PetscReal, parameter :: Cr = 1.00374d-7    ! Pa^-1
  PetscReal, parameter :: P0 = 1.249d7       ! Pa
  PetscReal, parameter :: phi0 = 7.382392d-3 ! -
  PetscReal, parameter :: initiating_pressure = 1.d12    ! Pa, dummy value
  PetscReal, parameter :: fully_altered_pressure = 1.d12 ! Pa, dummy value
  PetscReal, parameter :: fully_altered_porosity = 1.d-1 ! -, dummy value
  PetscReal, parameter :: fully_alt_compressibility = 1.d-1 ! dummy value
  PetscInt :: i,k
  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscInt :: values(8)
  PetscReal :: tolerance
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time
  character(len=MAXWORDLENGTH) :: pass_fail, filename_out, input_filename

  input_filename = trim(this%filename)
  ! Read input file
  realization%option%io_buffer = 'Running compressibility unit test, input '//trim(this%filename)
  call PrintMsg(realization%option)
  i = 1
  open(action='read', file=trim(input_filename), iostat=rc_in, &
       newunit=fu_in)
  if (rc_in /= 0) then
    realization%option%io_buffer = 'File read error: compressibility unit test.'//trim(this%filename)
    call PrintErrMsg(realization%option)
  endif
  read(fu_in, *) ! skip header line
  do
    read (fu_in, *, iostat=rc_in) temp_material_name(i), &
                                  temp_material_id(i), &
                                  temp_liquid_pressure(i), &
                                  temp_correct_altered_porosity(i)
    if (rc_in /= 0) exit
    i = i + 1
    if (i > 150) exit
  enddo

  allocate(material_name(i-1))
  allocate(material_id(i-1))
  allocate(liquid_pressure(i-1))
  allocate(correct_altered_porosity(i-1))
  allocate(porosity(i-1))

  material_name(:) = temp_material_name(1:i-1)
  material_id(:) = temp_material_id(1:i-1)
  liquid_pressure(:)      = temp_liquid_pressure(1:i-1)
  correct_altered_porosity(:) = temp_correct_altered_porosity(1:i-1)

  do i = 1, size(liquid_pressure)
    call FractureCalcFracPorosity(Cr,phi0,initiating_pressure,fully_altered_pressure,P0, &
         liquid_pressure(i),fully_altered_porosity,fully_alt_compressibility,porosity(i))
  end do

  tolerance = this%tolerance

  call PMUnitTestWIPPFilenameInToOut(input_filename, filename_out)
  filename_out = trim(filename_out)
  open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
  call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
  write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                  time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
  write(fu_out,*)
  write(fu_out,'(a)') 'NOTE: The input file provided was:'
  write(fu_out,'(a,a)') '      ', trim(input_filename)
  write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                               tolerance, '.'
  write(fu_out,*)
  i = 0
  do k=1,size(liquid_pressure)
    write(fu_out,'(a,I3,a)') '||-----------TEST-#',k,'-----------------------&
                             &--------------||'
    write(fu_out,'(a)') '  [in]  material property name [-]:'
    write(fu_out,'(a)') material_name(k)
    write(fu_out,'(a)') '  [in]  liquid pressure [Pa]:'
    write(fu_out,'(d17.10)') liquid_pressure(k)

    write(fu_out,'(a)') '  [out]  altered porosity [-]:'
    write(fu_out,'(d17.10)') porosity(k)
    write(fu_out,'(a)') '  [correct]  altered porosity [-]:'
    write(fu_out,'(d17.10)') correct_altered_porosity(k)

    call CalcDiff(porosity(k),correct_altered_porosity(k),tolerance,pass_fail,i)
    write(fu_out,'(a)') trim(pass_fail)
    write(fu_out,*)
  enddo

  write(fu_out,'(a)') 'TEST SUMMARY:'
  if (i == 0) then
    write(fu_out,'(a)') ' All tests passed!'
  else
    write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
  endif

  close(fu_out)
  close(fu_in)

end subroutine PMUnitTestWIPPRunCompressibility

! ************************************************************************** !

subroutine PMUnitTestWIPPRunEOSGasDensity(this,realization)
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !

  use EOS_Gas_module

  implicit none

  class(pm_unittest_type) :: this
  type(realization_subsurface_type) :: realization

  PetscReal, allocatable :: temperature(:)
  PetscReal :: temp_temperature(150)
  PetscReal, allocatable :: pressure(:)
  PetscReal :: temp_pressure(150)
  PetscReal, allocatable :: den_kg(:)
  PetscReal, allocatable :: correct_den_kg(:)
  PetscReal :: temp_correct_den_kg(150)
  PetscErrorCode :: ierr
  PetscInt :: i,k
  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscInt :: values(8)
  PetscReal :: tolerance
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time
  character(len=MAXWORDLENGTH) :: pass_fail, filename_out, input_filename
  ierr = 0

  input_filename = trim(this%filename)
  ! Read input file
  realization%option%io_buffer = 'Running gas density unit test, input '//trim(this%filename)
  call PrintMsg(realization%option)
  i = 1
  open(action='read', file=trim(input_filename), iostat=rc_in, &
       newunit=fu_in)
  if (rc_in /= 0) then
    realization%option%io_buffer = 'File read error: gas density unit test.'//trim(this%filename)
    call PrintErrMsg(realization%option)
  endif
  read(fu_in, *) ! skip header line
  do
    read (fu_in, *, iostat=rc_in) temp_temperature(i), &
                                  temp_pressure(i), &
                                  temp_correct_den_kg(i)
    if (rc_in /= 0) exit
    i = i + 1
    if (i > 150) exit
  enddo

  allocate(temperature(i-1))
  allocate(pressure(i-1))
  allocate(correct_den_kg(i-1))
  allocate(den_kg(i-1))

  temperature(:) = temp_temperature(1:i-1)
  pressure(:) = temp_pressure(1:i-1)
  correct_den_kg(:)      = temp_correct_den_kg(1:i-1)

  do i = 1, size(pressure)
     call EOSGasDensity(temperature(i),pressure(i),den_kg(i),ierr)
     den_kg(i) = den_kg(i)
  end do

  tolerance = this%tolerance

  call PMUnitTestWIPPFilenameInToOut(input_filename, filename_out)
  filename_out = trim(filename_out)
  open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
  call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
  write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                  time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
  write(fu_out,*)
  write(fu_out,'(a)') 'NOTE: The input file provided was:'
  write(fu_out,'(a,a)') '      ', trim(input_filename)
  write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                               tolerance, '.'
  write(fu_out,*)
  i = 0
  do k=1,size(den_kg)
    write(fu_out,'(a,I3,a)') '||-----------TEST-#',k,'-----------------------&
                             &--------------||'
    write(fu_out,'(a)') '  [in]  temperature [C]:'
    write(fu_out,'(d17.10)') temperature(k)
    write(fu_out,'(a)') '  [out]  density [kg/m3] (BRAGFLO)'
    write(fu_out,'(d17.10)') den_kg(k)

    write(fu_out,'(a)') '  [correct]  density [kg/m3] (BRAGFLO):'
    write(fu_out,'(d17.10)') correct_den_kg(k)
    call CalcDiff(den_kg(k),correct_den_kg(k),tolerance,pass_fail,i)
    write(fu_out,'(a)') trim(pass_fail)
    write(fu_out,*)
  enddo

  write(fu_out,'(a)') 'TEST SUMMARY:'
  if (i == 0) then
    write(fu_out,'(a)') ' All tests passed!'
  else
    write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
  endif

  close(fu_out)
  close(fu_in)

end subroutine PMUnitTestWIPPRunEOSGasDensity

! ************************************************************************** !

subroutine PMUnitTestWIPPRunEOSWaterDensity(this,realization)
  !
  ! Tests water density (kg/m^3) for WIPP unit tests.
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/20/25
  !

  use EOS_Water_module

  implicit none
  class(pm_unittest_type) :: this
  type(realization_subsurface_type) :: realization

  PetscReal, allocatable :: temperature(:)
  PetscReal :: temp_temperature(150)
  PetscReal, allocatable :: pressure(:)
  PetscReal :: temp_pressure(150)
  PetscReal, allocatable :: den_kg(:)
  PetscReal :: den
  PetscReal, allocatable :: correct_den_kg(:)
  PetscReal :: temp_correct_den_kg(150)
  PetscErrorCode :: ierr
  PetscInt :: i,k
  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscInt :: values(8)
  PetscReal :: tolerance
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time
  character(len=MAXWORDLENGTH) :: pass_fail, filename_out, input_filename

  input_filename = trim(this%filename)
  ! Read input file
  realization%option%io_buffer = 'Running water density unit test, input '//trim(this%filename)
  call PrintMsg(realization%option)
  i = 1
  open(action='read', file=trim(input_filename), iostat=rc_in, &
       newunit=fu_in)
  if (rc_in /= 0) then
    realization%option%io_buffer = 'File read error: fracture unit test.'//trim(this%filename)
    call PrintErrMsg(realization%option)
  endif
  read(fu_in, *) ! skip header line
  do
    read (fu_in, *, iostat=rc_in) temp_temperature(i), &
                                  temp_pressure(i), &
                                  temp_correct_den_kg(i)
    if (rc_in /= 0) exit
    i = i + 1
    if (i > 150) exit
  enddo

  allocate(temperature(i-1))
  allocate(pressure(i-1))
  allocate(correct_den_kg(i-1))
  allocate(den_kg(i-1))

  temperature(:) = temp_temperature(1:i-1)
  pressure(:) = temp_pressure(1:i-1)
  correct_den_kg(:)      = temp_correct_den_kg(1:i-1)
  do i = 1, size(pressure)
    call EOSWaterDensity(temperature(i),pressure(i),den_kg(i),den,ierr)
  end do

  tolerance = this%tolerance

  call PMUnitTestWIPPFilenameInToOut(input_filename, filename_out)
  filename_out = trim(filename_out)
  open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
  call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
  write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                  time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
  write(fu_out,*)
  write(fu_out,'(a)') 'NOTE: The input file provided was:'
  write(fu_out,'(a,a)') '      ', trim(input_filename)
  write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                               tolerance, '.'
  write(fu_out,*)
  i = 0
  do k=1,size(den_kg)
    write(fu_out,'(a,I3,a)') '||-----------TEST-#',k,'-----------------------&
                             &--------------||'
    write(fu_out,'(a)') '  [in]  temperature [C]:'
    write(fu_out,'(d17.10)') temperature(k)
    write(fu_out,'(a)') '  [out]  density [kg/m3] (BRAGFLO)'
    write(fu_out,'(d17.10)') den_kg(k)

    write(fu_out,'(a)') '  [correct]  density [kg/m3] (BRAGFLO):'
    write(fu_out,'(d17.10)') correct_den_kg(k)
    call CalcDiff(den_kg(k),correct_den_kg(k),tolerance,pass_fail,i)
    write(fu_out,'(a)') trim(pass_fail)
    write(fu_out,*)
  enddo

  write(fu_out,'(a)') 'TEST SUMMARY:'
  if (i == 0) then
    write(fu_out,'(a)') ' All tests passed!'
  else
    write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
  endif

  close(fu_out)
  close(fu_in)

end subroutine PMUnitTestWIPPRunEOSWaterDensity

! ************************************************************************** !

subroutine CalcDiff(val, correct_val, tolerance,pass_fail,i)
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !

  implicit none

  PetscReal :: val
  PetscReal :: correct_val
  PetscReal :: tolerance, diff
  PetscInt :: i
  character(len=MAXWORDLENGTH) :: pass_fail

  diff = abs(correct_val-val)

  if (diff > tolerance*correct_val) then
    pass_fail = 'FAIL!'
    i = i + 1
  else
    pass_fail = 'pass'
  endif

end subroutine

! ************************************************************************** !

subroutine PMUnitTestWIPPFilenameInToOut(filename_in, filename_out)
  !
  ! Creates .out file from a .in file
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/20/25
  !

  implicit none

  character(len=MAXWORDLENGTH), intent(in) :: filename_in
  character(len=MAXWORDLENGTH), intent(out) :: filename_out

  integer :: dot_position

  ! Find the position of the last dot
  dot_position = index(filename_in, '.')

  ! Check if a dot was found
  if (dot_position > 0) then
    ! Extract the substring before the dot
    filename_out = trim(adjustl(filename_in(1:dot_position-1)))
  else
    ! If no dot is found, return the original filename
    filename_out = trim(adjustl(filename_in))
  end if
  filename_out = trim(filename_out) // '.out'

end subroutine PMUnitTestWIPPFilenameInToOut

! ************************************************************************** !

end module PM_Unit_Test_WIPP_class
