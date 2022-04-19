module Lookup_Table_module
  
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none
  
  private

  ! Variable extrapolation types
  PetscInt, parameter, public :: VAR_EXTRAP_CONST_VAL = 1
  PetscInt, parameter, public :: VAR_EXTRAP_CONST_GRAD = 2
  ! Variable interpolation types
  PetscInt, parameter :: VAR_INTERP_LINEAR = 1
  PetscInt, parameter :: VAR_INTERP_X_LINLOG = 2
  ! 3D Interpolation options
  PetscInt, parameter :: INTERP_3D_LP = 1 ! Lagrange polynomials
  PetscInt, parameter :: INTERP_3D_TL = 2 ! Trilinear
  
  type, abstract, public :: lookup_table_base_type
    PetscInt :: dim
    PetscInt :: dims(3)
    PetscReal, pointer :: data(:)
    PetscReal, pointer :: var_data(:,:)
    type(lookup_table_var_ptr_type), pointer :: var_array(:)
    type(lookup_table_var_list_type), pointer :: vars    
    class(lookup_table_axis_type), pointer :: axis1
  contains
    procedure(LookupTableEvaluateDummy), deferred, public :: Sample
    procedure(LookupTableValAndGradDummy),deferred,public :: SampleAndGradient
    procedure(LookupTableAxesAreSMIncDummy), deferred,public :: AxesAreSMInc
    procedure, public :: LookupTableVarConvFactors  
    procedure, public :: LookupTableVarsInit
    procedure, public :: LookupTableVarIsPresent
    procedure, public :: LookupTableVarIsSMInc
    procedure, public :: CreateAddLookupTableVar
    procedure, public :: VarDataRead
    procedure, public :: VarDataReverse
    procedure, public :: VarPointAndUnitConv
    procedure, public :: SetupConstGradExtrap
    procedure, public :: SetupConstValExtrap
    procedure, public :: LookupTableVarInitGradients
    procedure, public :: SetupVarLinLogInterp
    procedure, public :: SetupVarUserUnits
  end type lookup_table_base_type
  
  type, public, extends(lookup_table_base_type) :: lookup_table_uniform_type
    class(lookup_table_axis_type), pointer :: axis2
    class(lookup_table_axis_type), pointer :: axis3
  contains
    procedure, public :: Sample => LookupTableEvaluateUniform
    procedure, public :: SampleAndGradient => ValAndGradUniform
    procedure, public :: AxesAreSMInc => AxesAreSMIncUniform
  end type lookup_table_uniform_type

  type, public :: irregular_array_type
    ! implements an irregular (or jagged) array
    PetscReal, pointer :: data(:)
  end type irregular_array_type
  
  type, public, extends(lookup_table_base_type) :: lookup_table_general_type
    class(lookup_table_axis2_general_type), pointer :: axis2
    class(lookup_table_axis3_general_type), pointer :: axis3
    PetscInt :: mode ! interpolation mode
    PetscReal, pointer :: data_references(:,:,:) ! reference data values from interpolation indices
    type(irregular_array_type), allocatable :: partition(:) ! if allocated, data is partitioned
  contains
    procedure, public :: Sample => LookupTableEvaluateGeneral
    procedure, public :: SampleAndGradient => ValAndGradGeneral
    procedure, public :: AxesAreSMInc => AxesAreSMIncGeneral
  end type lookup_table_general_type
  
  type, public :: lookup_table_axis_type
    PetscInt :: itype
    PetscInt :: saved_index
    PetscInt :: saved_index1
    PetscReal, pointer :: values(:)
  end type lookup_table_axis_type
  
  type, public, extends(lookup_table_axis_type) :: lookup_table_axis2_general_type
    PetscInt :: saved_index2
  end type lookup_table_axis2_general_type

  type, public, extends(lookup_table_axis_type) :: lookup_table_axis3_general_type
    PetscBool :: extrapolate ! extrapolation from axis3 required
    PetscInt :: num_partitions ! number of axis3 partitions
    PetscInt, allocatable :: bounds(:) ! bounds of axis3 partitions
    PetscInt, allocatable :: saved_index_partition(:,:) ! index of partition per i, j coordinate
    PetscInt, allocatable :: saved_indices3(:,:,:) ! left/right index k per i, j coordinate
    type(irregular_array_type), allocatable :: partition(:) ! if allocated, axis 3 is partitioned
  end type lookup_table_axis3_general_type

  type, public :: lookup_table_var_type
    PetscInt :: id
    PetscInt :: iname
    PetscInt :: data_idx
    PetscInt :: extrapolation_itype
    PetscInt :: interp_type
    character(len=MAXWORDLENGTH) :: internal_units
    character(len=MAXWORDLENGTH) :: user_units
    PetscReal :: conversion_factor
    PetscReal, pointer :: data(:)
    PetscReal :: sample
    PetscReal, pointer :: sample_grad(:)
    type(lookup_table_var_type), pointer :: next
  end type  lookup_table_var_type

  type, public :: lookup_table_var_ptr_type
    type(lookup_table_var_type), pointer :: ptr
  end type lookup_table_var_ptr_type

  type, public :: lookup_table_var_list_type
    PetscInt :: num_lookup_table_vars
    type(lookup_table_var_type), pointer :: first
    type(lookup_table_var_type), pointer :: last
    !type(lookup_table_var_type), pointer :: array(:)
  end type lookup_table_var_list_type
  
  abstract interface
    function LookupTableEvaluateDummy(this,lookup1,lookup2,lookup3)
      import lookup_table_base_type
      implicit none
      class(lookup_table_base_type) :: this
      PetscReal :: lookup1
      PetscReal, optional :: lookup2
      PetscReal, optional :: lookup3
      PetscReal :: LookupTableEvaluateDummy
    end function LookupTableEvaluateDummy
    
    subroutine LookupTableValAndGradDummy(this,var_iname,lookup1,lookup2,lookup3)
      import lookup_table_base_type
      implicit none
      class(lookup_table_base_type) :: this
      PetscInt, intent(in) :: var_iname
      PetscReal :: lookup1
      PetscReal, optional :: lookup2
      PetscReal, optional :: lookup3
    end subroutine LookupTableValAndGradDummy

    subroutine LookupTableAxesAreSMIncDummy(this,AxisIsSMInc)
      import lookup_table_base_type
      implicit none
      class(lookup_table_base_type) :: this
      PetscBool, intent(out) :: AxisIsSMInc(this%dim)
    end subroutine LookupTableAxesAreSMIncDummy

  end interface
  
  interface LookupTableTest
    module procedure LookupTableTest1D
    module procedure LookupTableTest2D
  end interface

  interface LookupTableCreateGeneral
    module procedure LookupTableCreateGeneralDim
    module procedure LookupTableCreateGeneralNoDim  
  end interface

  interface CreateLookupTableVar
    module procedure CreateLookupTableVar1
    module procedure CreateLookupTableVar2
  end interface

  interface LookupTableDestroy
    module procedure LookupTableUniformDestroy
    module procedure LookupTableGeneralDestroy
  end interface

  public :: LookupTableCreateUniform, &
            LookupTableCreateGeneral, &
            InverseLookupTableCreateGen, &
            LookupTableAxisInit, &
            LookupTableDestroy, &
            !InverseLookupTableGenDestroy, &
            LookupTableTest, &
            LookupTableVarInitList, &
            CreateLookupTableVar, &
            LookupTableVarAddToList, &
            LookupTableVarListDestroy            
  
contains

! ************************************************************************** !

subroutine LookupTableBaseInit(lookup_table)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 

  implicit none
  
  class(lookup_table_base_type) :: lookup_table
  
  lookup_table%dim = 0
  lookup_table%dims = 0
  nullify(lookup_table%data)
  nullify(lookup_table%var_data)
  nullify(lookup_table%var_array)
  nullify(lookup_table%vars)
  nullify(lookup_table%axis1)

end subroutine LookupTableBaseInit

! ************************************************************************** !

function LookupTableCreateUniform(dim)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 

  implicit none
  
  PetscInt :: dim
  
  class(lookup_table_uniform_type), pointer :: LookupTableCreateUniform

  class(lookup_table_uniform_type), pointer :: lookup_table
  
  allocate(lookup_table)
  call LookupTableBaseInit(lookup_table)
  lookup_table%dim = dim
  nullify(lookup_table%axis2)
  nullify(lookup_table%axis3)

  allocate(lookup_table%axis1)
  call LookupTableAxisInit(lookup_table%axis1)
  if (dim > 1) then
    allocate(lookup_table%axis2)
    call LookupTableAxisInit(lookup_table%axis2)
  endif
  if (dim > 2) then
    allocate(lookup_table%axis3)
    call LookupTableAxisInit(lookup_table%axis3)
  endif
  
  LookupTableCreateUniform => lookup_table

end function LookupTableCreateUniform

! ************************************************************************** !

function LookupTableCreateGeneralDim(dim)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14 - 
  ! Modified by Paolo Orsini: change of function name to add new constructor
  ! Date: 05/02/18
  ! Modified by Alex Salazar: accommodate 3D interpolation
  ! Date: 02/16/2022
  ! 

  implicit none
  
  PetscInt :: dim
  
  class(lookup_table_general_type), pointer :: LookupTableCreateGeneralDim

  class(lookup_table_general_type), pointer :: lookup_table
  
  allocate(lookup_table)
  call LookupTableBaseInit(lookup_table)
  lookup_table%dim = dim
  nullify(lookup_table%axis2)
  nullify(lookup_table%axis3)
  allocate(lookup_table%axis1)
  call LookupTableAxisInit(lookup_table%axis1)
  if (dim > 1) then
    allocate(lookup_table%axis2)
    call LookupTableAxisInit(lookup_table%axis2)
    lookup_table%axis2%saved_index2 = 1
  endif
  if (dim > 2) then
    allocate(lookup_table%axis3)
    call LookupTableAxisInit(lookup_table%axis3)
    lookup_table%mode = INTERP_3D_LP
    lookup_table%axis3%num_partitions = UNINITIALIZED_INTEGER
    lookup_table%axis3%extrapolate = PETSC_FALSE
    nullify(lookup_table%data_references)
  endif
  
  LookupTableCreateGeneralDim => lookup_table

end function LookupTableCreateGeneralDim

! ************************************************************************** !

function LookupTableCreateGeneralNoDim()
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 05/02/18
  ! 
  ! contructor before knowing the lookup table dims

  implicit none
  
  class(lookup_table_general_type), pointer :: LookupTableCreateGeneralNoDim

  class(lookup_table_general_type), pointer :: lookup_table
  
  allocate(lookup_table)
  call LookupTableBaseInit(lookup_table)
  nullify(lookup_table%axis2)
  nullify(lookup_table%axis3)
  
  LookupTableCreateGeneralNoDim => lookup_table

end function LookupTableCreateGeneralNoDim

! ************************************************************************** !

function InverseLookupTableCreateGen(lookup_table,new_axis_var_iname,option)
  ! 
  ! Author: Paolo Orsini
   
  ! Create inverse lookup given: 
  !   - a general lookup table 
  !   - the new axis var, assumed to be stored in var_data and with monotonic
  !     beahaviour
  ! Note only 1D inverse lookup are currently supported
  !
  ! Date: 08/15/18
  ! 

  use Option_module

  implicit none
  
  class(lookup_table_general_type), pointer :: InverseLookupTableCreateGen
  class(lookup_table_general_type), intent(in) :: lookup_table
  PetscInt, intent(in) ::  new_axis_var_iname
  type (option_type) :: option
  
  class(lookup_table_general_type), pointer :: inv_lookup_table
  type(lookup_table_var_type), pointer :: NewInvLookupVar
  PetscBool :: reverse
  PetscInt :: data_idx, i_tmp, var_idx_max
  
  allocate(inv_lookup_table)
  call LookupTableBaseInit(inv_lookup_table)
  Inv_lookup_table%dim = lookup_table%dim
  Inv_lookup_table%dims = lookup_table%dims
  nullify(inv_lookup_table%axis2)
  allocate(inv_lookup_table%axis1)
  call LookupTableAxisInit(inv_lookup_table%axis1)
  if (Inv_lookup_table%dim > 1) then
    option%io_buffer = "only 1D inverse lookup are supported"
    call PrintErrMsg(option)
  endif
  
  nullify(inv_lookup_table%data)
  nullify(inv_lookup_table%var_data)
  
  allocate(inv_lookup_table%var_data( size(lookup_table%var_data,1), &
                                      size(lookup_table%var_data,2)  ) )
                                      
  !determine if table inversion is required or not
  data_idx = lookup_table%var_array(new_axis_var_iname)%ptr%data_idx
  var_idx_max = size(lookup_table%axis1%values(:))
  if (lookup_table%var_data(data_idx,1) > &
      lookup_table%var_data(data_idx,var_idx_max)) then
    reverse = PETSC_TRUE
  else if (lookup_table%var_data(data_idx,1) < &
      lookup_table%var_data(data_idx,var_idx_max)) then
    reverse = PETSC_FALSE  
  else
    option%io_buffer = "InverseLookupTableCreateGen: cannot create " // &
      "inverse gradient lookup as the new axis var selected is not monotonic "
    call PrintErrMsg(option)
  end if
    
  !assumes that any unit conversion on lookup_table%var_data is already done
  if (reverse) then
    data_idx =  1
    do i_tmp = size(lookup_table%var_data,2),1,-1
      inv_lookup_table%var_data(:,data_idx) = lookup_table%var_data(:,i_tmp)
      data_idx = data_idx + 1    
    end do 
  else
    inv_lookup_table%var_data = lookup_table%var_data
  end if
    
  call inv_lookup_table%LookupTableVarsInit(size(lookup_table%var_array(:)))
  do i_tmp = 1,size(lookup_table%var_array(:))
    if ( associated(lookup_table%var_array(i_tmp)%ptr) ) then
      NewInvLookupVar => CreateLookupTableVar( &
                                          lookup_table%var_array(i_tmp)%ptr)
      nullify(NewInvLookupVar%data)
      data_idx = lookup_table%var_array(i_tmp)%ptr%data_idx
      NewInvLookupVar%data_idx = data_idx
      NewInvLookupVar%data => inv_lookup_table%var_data(data_idx,:)
      inv_lookup_table%var_array(i_tmp)%ptr => NewInvLookupVar
      call LookupTableVarAddToList(NewInvLookupVar,inv_lookup_table%vars)
    end if
  end do

  
  allocate(inv_lookup_table%axis1)
  data_idx = inv_lookup_table%var_array(new_axis_var_iname)%ptr%data_idx
  !copying axis1%values (not pointing) for consistency with all other tables
  allocate(inv_lookup_table%axis1%values( &
                          size(inv_lookup_table%var_data(data_idx,:))) )
  inv_lookup_table%axis1%values = inv_lookup_table%var_data(data_idx,:)
  
  InverseLookupTableCreateGen => inv_lookup_table

end function InverseLookupTableCreateGen

! ************************************************************************** !

subroutine LookupTableAxisInit(axis)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 

  implicit none
  
  class(lookup_table_axis_type) :: axis
  
  axis%itype = 0
  axis%saved_index = 1
  axis%saved_index1 = 1
  nullify(axis%values)
  
end subroutine LookupTableAxisInit

! ************************************************************************** !

function LookupTableEvaluateUniform(this,lookup1,lookup2,lookup3)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3 
  
  PetscReal :: LookupTableEvaluateUniform

  call LookupTableIndexUniform(this,lookup1,lookup2,lookup3)
  if (present(lookup3)) then
!    call LookupTableInterpolate3DUniform(this,lookup1,lookup2,lookup3,LookupTableEvaluateUniform)
  else if (present(lookup2)) then
    call LookupTableInterpolate2DUniform(this,lookup1,lookup2,LookupTableEvaluateUniform)
  else
    call LookupTableInterpolate1D(this,lookup1,LookupTableEvaluateUniform)
  endif
  
end function LookupTableEvaluateUniform

! ************************************************************************** !

subroutine ValAndGradUniform(this,var_iname,lookup1,lookup2,lookup3)
  ! 
  ! Computes value and gradient for given coordinates
  !
  ! Author: Paolo Orsini
  ! Date: 05/11/18
  ! 
  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3

  call LookupTableIndexUniform(this,lookup1,lookup2,lookup3)
  if (present(lookup3)) then
    !3D interpolation and gradient computation not yet supported
  else if (present(lookup2)) then
    call InterpExtrapGrad2DUniform(this,var_iname,lookup1,lookup2)
  else
    call InterpExtrapGrad1D(this,var_iname,lookup1)
  endif
  
end subroutine ValAndGradUniform

! ************************************************************************** !

subroutine AxesAreSMIncUniform(this,AxisIsSMInc)

  use Utility_module, only : ArrayIsSMonotonicInc

  implicit none

  class(lookup_table_uniform_type) :: this
  PetscBool, intent(out) :: AxisIsSMInc(this%dim)

  AxisIsSMInc(:) = PETSC_FALSE
  if (associated(this%axis1%values)) then
    AxisIsSMInc(ONE_INTEGER) = ArrayIsSMonotonicInc(this%axis1%values)
  end if

  if (this%dim == 2) then
    if (associated(this%axis2%values)) then
      AxisIsSMInc(TWO_INTEGER) = ArrayIsSMonotonicInc(this%axis2%values)
    end if
  end if

end subroutine AxesAreSMIncUniform

! ************************************************************************** !

function LookupTableEvaluateGeneral(this,lookup1,lookup2,lookup3)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  class(lookup_table_general_type) :: this
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3 
  
  PetscReal :: LookupTableEvaluateGeneral

  call LookupTableIndexGeneral(this,lookup1,lookup2,lookup3)
  if (present(lookup3)) then
    select case(this%mode)
    !-------------------------------------
      case(INTERP_3D_LP)
        call LookupTableInterpolate3DLP(this,lookup1,lookup2,lookup3, &
                                          LookupTableEvaluateGeneral)
    !-------------------------------------
      case(INTERP_3D_TL)
        call LookupTableInterpolate3DTrilinear(this,lookup1,lookup2,lookup3, &
                                                 LookupTableEvaluateGeneral)
    !-------------------------------------
      case default 
!    call LookupTableInterpolate3DGeneral(this,lookup1,lookup2,lookup3,LookupTableEvaluateGeneral)
    end select
  else if (present(lookup2)) then
    call LookupTableInterpolate2DGeneral(this,lookup1,lookup2,LookupTableEvaluateGeneral)
  else
    call LookupTableInterpolate1D(this,lookup1,LookupTableEvaluateGeneral)
  endif
  
end function LookupTableEvaluateGeneral

! ************************************************************************** !

subroutine ValAndGradGeneral(this,var_iname,lookup1,lookup2,lookup3)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/21/18
  ! 
  implicit none
  
  class(lookup_table_general_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3  
  
  call LookupTableIndexGeneral(this,lookup1,lookup2,lookup3)
  
  if (present(lookup3)) then
    ! 3D general lookup table not supported
  else if (present(lookup2)) then
    call InterpExtrapGradGeneral2D(this,var_iname,lookup1,lookup2)
  else
    call InterpExtrapGrad1D(this,var_iname,lookup1)
  end if    
  
  
end subroutine ValAndGradGeneral

! ************************************************************************** !

subroutine AxesAreSMIncGeneral(this,AxisIsSMInc)

  use Utility_module, only : ArrayIsSMonotonicInc

  implicit none

  class(lookup_table_general_type) :: this
  PetscBool, intent(out) :: AxisIsSMInc(this%dim)

  PetscInt :: i_tmp
  PetscInt :: i1, i2

  AxisIsSMInc(:) = PETSC_FALSE

  if (associated(this%axis1%values)) then
    AxisIsSMInc(ONE_INTEGER) = ArrayIsSMonotonicInc(this%axis1%values)
  end if

  if (this%dim == 2) then
    if (associated(this%axis2%values)) then
      AxisIsSMInc(TWO_INTEGER) = PETSC_TRUE
      do i_tmp = 1,this%dims(1)
        i1 = (i_tmp - 1) * this%dims(2) + 1
        i2 = i_tmp * this%dims(2)
        AxisIsSMInc(TWO_INTEGER) = &
            ArrayIsSMonotonicInc(this%axis2%values(i1:i2))
      end do
    end if
  end if

end subroutine AxesAreSMIncGeneral

! ************************************************************************** !

subroutine LookupTableIndexUniform(this,lookup1,lookup2,lookup3)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3 
  
  call LookupTableAxisIndexUniform(this%axis1,lookup1)
  if (associated(this%axis2)) then
    call LookupTableAxisIndexUniform(this%axis2,lookup2)
  endif
  if (associated(this%axis3)) then
    call LookupTableAxisIndexUniform(this%axis3,lookup3)
  endif

end subroutine LookupTableIndexUniform

! ************************************************************************** !

subroutine LookupTableIndexGeneral(this,lookup1,lookup2,lookup3)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  ! Modified by: Alex Salazar III
  ! Date: 02/18/2022
  !
  implicit none
  
  class(lookup_table_general_type) :: this
  PetscReal :: lookup1
  PetscReal, optional :: lookup2
  PetscReal, optional :: lookup3 
  
  PetscInt :: ja, jb
  PetscInt :: istart, iend
  class(lookup_table_axis2_general_type), pointer :: axis2
  class(lookup_table_axis3_general_type), pointer :: axis3
  
  ! axis 1 corresponds to the j dim below
  call LookupTableAxisIndexGeneral(lookup1,this%axis1%values, &
                                   this%axis1%saved_index)
  
  if (associated(this%axis3)) then

    call LookupTableIndexAxis3(this,lookup1,lookup2,lookup3)

  elseif (associated(this%axis2)) then

    axis2 => this%axis2

    ja = this%axis1%saved_index
    if (ja > 0) then
      jb = max(min(ja+1,this%dims(1)),1)
    else
      ja = 1
      jb = 1
    endif

    iend = ja*this%dims(2)
    istart = iend - this%dims(2) + 1
    call LookupTableAxisIndexGeneral(lookup2,axis2%values(istart:iend), &
                                     axis2%saved_index)
    if (ja /= jb) then
      iend = jb*this%dims(2)
      istart = iend - this%dims(2) + 1
      call LookupTableAxisIndexGeneral(lookup2,axis2%values(istart:iend), &
                                       axis2%saved_index2)
    else 
      axis2%saved_index2 = axis2%saved_index
    endif
  endif

end subroutine LookupTableIndexGeneral

! ************************************************************************** !

subroutine LookupTableIndexAxis3(this,lookup1,lookup2,lookup3)
  !
  ! Identifies indices for 3D interpolation for given lookups
  !
  ! Author: Alex Salazar III
  ! Date: 02/21/2022
  !
  implicit none
  ! ----------------------------------
  class(lookup_table_general_type) :: this
  PetscReal :: lookup1
  PetscReal :: lookup2
  PetscReal :: lookup3
  ! ----------------------------------
  class(lookup_table_axis_type), pointer :: axis1
  class(lookup_table_axis2_general_type), pointer :: axis2
  class(lookup_table_axis3_general_type), pointer :: axis3
  PetscInt :: i, j, k, m ! iterators
  PetscInt :: li, lj, lk ! array lengths
  PetscInt :: i1, j1, k1 ! left bounds
  PetscInt :: i2, j2, k2 ! right bounds
  PetscBool :: iextrp ! lookup1 requires extrapolation from axis 1
  PetscBool :: jextrp ! lookup2 requires extrapolation from axis 2
  PetscBool :: kextrp ! lookup3 requires extrapolation from axis 3
  PetscInt :: pselect ! relevant partition of axis3
  PetscReal, pointer :: v3(:) ! subset of axis3 values
  PetscInt :: kstart, kend    ! start and end of axis3 array to interpolate
  PetscInt :: indices1(4,2) ! first set of indices used for axis3 indexing
  PetscInt :: indices2(8,3) ! second set of indices used for data indexing
  ! ----------------------------------

  axis1 => this%axis1
  axis2 => this%axis2
  axis3 => this%axis3

  iextrp = PETSC_FALSE
  jextrp = PETSC_FALSE
  kextrp = PETSC_FALSE

  li = size(axis1%values)
  lj = size(axis2%values)
  lk = size(axis3%values)

  ! axis1 indices (pivot variable)
  if (li == 1) then
    ! only one value given for axis1
    i1 = 1
    i2 = 1
  else
    if(lookup1 <= axis1%values(1)) then
      i2 = 2   ! extrapolation - use first and second points
    elseif (lookup1 > axis1%values(li)) then
      i2 = li  ! extrapolation - use final and penultimate points
      iextrp = PETSC_TRUE
    else
      ! interpolation
      do i = 2, li
        if (axis1%values(i - 1) < lookup1 .and. axis1%values(i) >= lookup1) then
          i2 = i
          exit
        endif
      enddo
    endif
    i1 = i2 - 1
  endif

  ! axis2 indices (pivot variable)
  if (lj == 1) then
    ! only one value given for axis2
    j1 = 1
    j2 = 1
  else
    if(lookup2 <= axis2%values(1)) then
      j2 = 2   ! extrapolation - first and second points will be used
    elseif (lookup2 > axis2%values(lj)) then
      j2 = lj  ! extrapolation - final and penultimate points will be used
      jextrp = PETSC_TRUE
    else
      ! interpolation
      do j = 2, lj
        if (axis2%values(j - 1) < lookup2 .and. axis2%values(j) >= lookup2) then
          j2 = j
          exit
        endif
      enddo
    endif
    j1 = j2 - 1
  endif

  ! save i, j indices
  axis1%saved_index  = i1
  axis1%saved_index1 = i2
  axis2%saved_index  = j1
  axis2%saved_index2 = j2

  ! list i, j combinations
  indices1(1,:) = (/i1, j1/)
  indices1(2,:) = (/i1, j2/)
  indices1(3,:) = (/i2, j1/)
  indices1(4,:) = (/i2, j2/)

  ! axis3 (independent variables) 
  !   indexing is based on i, j combinations
  if (.not. allocated(axis3%saved_indices3)) then
    allocate(axis3%saved_indices3(li,lj,2))
    axis3%saved_indices3 = 0
  endif
  kstart = 0
  kend = 0
  pselect = 0
  if (allocated(axis3%bounds)) then
    ! ---> axis3 is has defined partitions (non-rectangular)
    ! allocate saved partition indices if needed
    if (.not. allocated(axis3%saved_index_partition)) then
      allocate(axis3%saved_index_partition(li,lj))
      axis3%saved_index_partition = 0
    endif

    ! loop through i, j coordinates
    do m = 1, 4
      i = indices1(m,1)
      j = indices1(m,2)

      pselect = this%dims(2)*(i - 1) + j

      axis3%saved_index_partition(i, j) = pselect

      v3 => axis3%partition(pselect)%data
      lk = size(v3) ! redfine length

      ! axis3 indices
      if (lk == 1) then
        ! only one value given for axis3
        k1 = 1
        k2 = 1
      else
        if(lookup3 <= v3(1)) then
          k2 = 2   ! extrapolation - first and second points will be used
        elseif (lookup3 > v3(lk)) then
          k2 = lk  ! extrapolation - final and penultimate points will be used
          kextrp = PETSC_TRUE
        else
          ! interpolation
          do k = 2, lk
            if (v3(k - 1) < lookup3 .and. v3(k) >= lookup3) then
              k2 = k
              exit
            endif
          enddo
        endif
        k1 = k2 - 1
      endif

      indices2(m,:)   = (/i, j, k1/)
      indices2(m+4,:) = (/i, j, k2/)

      axis3%saved_indices3(i, j, 1) = k1
      axis3%saved_indices3(i, j, 2) = k2

      if (associated(v3)) nullify(v3)
    enddo

  else
    ! ---> axis3 is described by the dim(3) value (rectangular)
    ! loop through i, j coordinates
    do m = 1, 4
      i = indices1(m,1)
      j = indices1(m,2)

      kstart = (((this%dims(2)*(i - 1) + j) - 1) * this%dims(3)) + 1
      kend = ((this%dims(2)*(i - 1) + j)) * this%dims(3)

      ! axis3 indices
      if (lk == 1) then
        ! only one value given for axis3
        k1 = 1
        k2 = 1
      else
        if(lookup3 <= axis3%values(1)) then
          k2 = 2   ! extrapolation - first and second points will be used
        elseif (lookup3 > axis3%values(lk)) then
          k2 = lk  ! extrapolation - final and penultimate points will be used
          kextrp = PETSC_TRUE
        else
          ! interpolation
          do k = kstart + 1, kend
            if (axis3%values(k - 1) < lookup3 .and. &
                axis3%values(k) >= lookup3) then
              k2 = k
              exit
            endif
          enddo
        endif
        k1 = k2 - 1
      endif

      indices2(m,:)   = (/i, j, k1/)
      indices2(m+4,:) = (/i, j, k2/)

      axis3%saved_indices3(i, j, 1) = k1
      axis3%saved_indices3(i, j, 2) = k2

    enddo

  endif

  ! assign reference data points to the indices
  nullify(this%data_references)
  allocate(this%data_references(li,lj,lk))
  this%data_references = 0.0d0
  if (allocated(this%partition)) then
    do m = 1, 8
      i = indices2(m,1)
      j = indices2(m,2)
      k = indices2(m,3)
      pselect = axis3%saved_index_partition(i, j)
      this%data_references(i, j, k) = this%partition(pselect)%data(k)
    enddo
  else
    do m = 1, 8
      i = indices2(m,1)
      j = indices2(m,2)
      k = indices2(m,3)
      this%data_references(i, j, k) = this%data(k)
    enddo
  endif

  if (kextrp) this%axis3%extrapolate = PETSC_TRUE

end subroutine LookupTableIndexAxis3

! ************************************************************************** !

subroutine LookupTableAxisIndexUniform(this,lookup1)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  class(lookup_table_axis_type) :: this
  PetscReal :: lookup1
  
  PetscInt :: i1
  PetscInt :: size1
  PetscReal :: begin1
  
  size1 = size(this%values)
  begin1 = this%values(1)
  i1 = int((lookup1 - begin1) / (this%values(size1) - begin1) * (size1-1) + 1)
  ! truncate i1 to zero indicating the value is below the range specified
  i1 = max(min(i1,size1),0)
  this%saved_index = i1

end subroutine LookupTableAxisIndexUniform

! ************************************************************************** !

subroutine LookupTableAxisIndexGeneral(lookup1,values,saved_index)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  implicit none
  
  PetscReal :: lookup1
  PetscInt :: saved_index
  PetscReal :: values(:)
  
  PetscInt :: i1
  PetscInt :: j1
  PetscInt :: mid1
  PetscInt :: sizei
 
  sizei = size(values)
  if (lookup1 < values(1)) then
    saved_index = 0
    return
  else if (lookup1 > values(sizei)) then
    saved_index = sizei
    return
  endif
  i1 = max(min(saved_index,sizei-1),1)
  if (lookup1 > values(i1+1) .or. &
      lookup1 < values(i1)) then
    ! move either up or down array
    if (lookup1 > values(i1+1)) then
      i1 = i1+1
      j1 = sizei
    else 
      j1 = i1
      i1 = 1
    endif
    do
      mid1 = (j1+i1) / 2
      if (lookup1 > values(mid1)) then
        i1 = mid1
      else
        j1 = mid1
      endif
      if (j1-i1 <= 1) exit
    enddo
  endif
  saved_index = i1

end subroutine LookupTableAxisIndexGeneral

! ************************************************************************** !

subroutine LookupTableInterpolate1D(this,lookup1,result)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module, only : Interpolate
  
  implicit none
  
  class(lookup_table_base_type) :: this
  PetscReal :: lookup1
  PetscReal :: result
  
  PetscInt :: i1, j1
  
  i1 = this%axis1%saved_index
  if (i1 > 0) then
    j1 = max(min(i1+1,this%dims(1)),1)
    call Interpolate(this%axis1%values(j1),this%axis1%values(i1),lookup1, &
                     this%data(j1),this%data(i1),result)
  else ! catch values below axis range
    result = this%data(1)
  endif
  
end subroutine LookupTableInterpolate1D

! ************************************************************************** !

subroutine InterpExtrapGrad1D(this,var_iname,lookup1)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/12/18
  ! 
  use Utility_module, only : Interpolate, GradientLinear
  
  implicit none
  
  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscReal, intent(in) :: lookup1
  
  PetscReal :: x1, x2, x1_grad, x2_grad
  PetscInt :: i1, i2, i1_grad, i2_grad
  PetscInt :: var_idx, sizei
  
  i1 = this%axis1%saved_index
  sizei = this%dims(1)
  
  if (i1 > 0) then
    i2 = max(min(i1+1,this%dims(1)),1)
  else
    i1=1
    i2=1
  end if    
  
  if (i1 == i2) then !end points
    if (i1 == 1) then
      i1_grad = 1
      i2_grad = 2
    else if (i1 == sizei) then
      i1_grad = sizei - 1
      i2_grad = sizei
    end if
    !x1 = this%axis1%values(i1)
    x1_grad = this%axis1%values(i1_grad)
    x2_grad = this%axis1%values(i2_grad)
    !x_fract_grad = (lookup1 - x1_grad) / ( x2_grad - x1_grad )
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then
          call Var1DExtrapolate(this,var_idx,i1,i1_grad,i2_grad, &
                                      x1_grad,x2_grad,lookup1)
        end if     
      end do
    else 
      call Var1DExtrapolate(this,var_iname,i1,i1_grad,i2_grad, &
                                  x1_grad,x2_grad,lookup1)
    end if
  else ! away from end points
    x1 = this%axis1%values(i1)
    x2 = this%axis1%values(i2)
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then
          call Var1DInterpolate(this,var_idx,i1,i2,x1,x2,lookup1)
        end if     
      end do
    else
      call Var1DInterpolate(this,var_iname,i1,i2,x1,x2,lookup1)
    end if    
  end if  
  
end subroutine InterpExtrapGrad1D

! ************************************************************************** !

subroutine Var1DExtrapolate(this,var_iname,i1,i1_grad,i2_grad, &
                            x1_grad,x2_grad,lookup)
  ! 
  ! Author: Paolo Orsini
  ! Date: 02/06/18
  !   
  implicit none

  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: i1,i1_grad,i2_grad
  PetscReal, intent(in) :: x1_grad,x2_grad
  PetscReal, intent(in) :: lookup
          
  if ( this%var_array(var_iname)%ptr%extrapolation_itype == &          
       VAR_EXTRAP_CONST_GRAD ) then
    call Var1DInterpolate(this,var_iname,i1_grad,i2_grad, &
                                   x1_grad,x2_grad,lookup)
  else if( this%var_array(var_iname)%ptr%extrapolation_itype == &
         VAR_EXTRAP_CONST_VAL) then
    this%var_array(var_iname)%ptr%sample = &
                         this%var_array(var_iname)%ptr%data(i1)
    this%var_array(var_iname)%ptr%sample_grad(1) = 0.0d0
      
  end if

end subroutine Var1DExtrapolate

! ************************************************************************** !

subroutine Var1DInterpolate(this,var_iname,i1,i2,x1,x2,lookup)
  ! 
  ! Author: Paolo Orsini
  ! Date: 02/06/18
  ! 
  use Utility_module, only : Interpolate, GradientLinear
  
  implicit none
    
  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: i1,i2
  PetscReal, intent(in) :: x1,x2
  PetscReal, intent(in) :: lookup

  PetscReal :: val1, val2

  val1 = this%var_array(var_iname)%ptr%data(i1)
  val2 = this%var_array(var_iname)%ptr%data(i2)
  if (this%var_array(var_iname)%ptr%interp_type  == VAR_INTERP_X_LINLOG) then
     val1 = dlog(val1)
     val2 = dlog(val2)
  end if   
  call GradientLinear(x2,x1,val2,val1, &
                      this%var_array(var_iname)%ptr%sample_grad(1))
  call Interpolate(x2,x1,lookup,val2,val1, &
                   this%var_array(var_iname)%ptr%sample)
  if (this%var_array(var_iname)%ptr%interp_type  == VAR_INTERP_X_LINLOG) then
    this%var_array(var_iname)%ptr%sample = &
        dexp(this%var_array(var_iname)%ptr%sample)
    this%var_array(var_iname)%ptr%sample_grad(1) = &
        this%var_array(var_iname)%ptr%sample_grad(1) * &
        this%var_array(var_iname)%ptr%sample
  end if                 
  
end subroutine Var1DInterpolate

! ************************************************************************** !

subroutine LookupTableInterpolate2DUniform(this,lookup1,lookup2,result)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module, only : Interpolate, InterpolateBilinear

  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscReal :: lookup1
  PetscReal :: lookup2
  PetscReal :: result
  
  PetscInt :: i1, i2, j1, j2
  PetscReal :: x1, x2, y1, y2, z1, z2, z3, z4
  PetscInt :: sizei, sizej

  !  x1,y2,z3 ------ x2,y2,z4
  !     |               |
  !     |               |
  !     |   x,y         |
  !     |               |
  !  x1,y1,z1 ------ x2,y1,z2                   
  
  result = UNINITIALIZED_DOUBLE
  sizei = this%dims(1)
  sizej = this%dims(2)
  i1 = this%axis1%saved_index
  j1 = this%axis2%saved_index
  ! index axes
  if (i1 > 0) then
    i2 = max(min(i1+1,sizei),1)
  else
    i1 = 1
    i2 = 1
  endif
  if (j1 > 0) then
    j2 = max(min(j1+1,sizej),1)
  else
    j1 = 1
    j2 = 1
  endif
  if (i2 == i1) then
    if (j2 == j1) then
      ! corner of domain
      result = this%data(i1+(j1-1)*sizei)
    else
      y1 = this%axis2%values(j1)
      y2 = this%axis2%values(j2)
      z1 = this%data(i1+(j1-1)*sizei)
      z3 = this%data(i2+(j2-1)*sizei)
      call Interpolate(y2,y1,lookup2,z3,z1,result)
    endif
  else if (j2 == j1) then
    x1 = this%axis1%values(i1)
    x2 = this%axis1%values(i2)
    z1 = this%data(i1+(j1-1)*sizei)
    z2 = this%data(i2+(j2-1)*sizei)
    call Interpolate(x2,x1,lookup1,z2,z1,result)
  else
    x1 = this%axis1%values(i1)
    x2 = this%axis1%values(i2)
    y1 = this%axis2%values(j1)
    y2 = this%axis2%values(j2)
    z1 = this%data(i1+(j1-1)*sizei)
    z2 = this%data(i2+(j1-1)*sizei)
    z3 = this%data(i1+(j2-1)*sizei)
    z4 = this%data(i2+(j2-1)*sizei)
    result = InterpolateBilinear(lookup1,lookup2,x1,x2,y1,y2,z1,z2,z3,z4)
  endif
                   
end subroutine LookupTableInterpolate2DUniform

! ************************************************************************** !

subroutine InterpExtrapGrad2DUniform(this,var_iname,lookup1,lookup2)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/25/18
  ! 

  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscReal :: lookup1
  PetscReal :: lookup2
  class(lookup_table_axis_type), pointer :: axis1 => null()
  class(lookup_table_axis_type), pointer :: axis2 => null()

  PetscInt :: i1, i2, j1, j2
  PetscInt :: i1_grad, i2_grad, j1_grad, j2_grad
  PetscReal :: x1, x2
  PetscReal :: x1_grad, x2_grad
  PetscReal :: x_frac
  PetscInt :: sizei, sizej
  PetscInt :: var_idx

  !  x1,y2,z3 ------ x2,y2,z4
  !     |               |
  !     |               |
  !     |   x,y         |
  !     |               |
  !  x1,y1,z1 ------ x2,y1,z2


  axis1 => this%axis1
  axis2 => this%axis2
  
  sizei = this%dims(1)
  sizej = this%dims(2)
  i1 = axis1%saved_index
  j1 = axis2%saved_index
  ! index axes
  if (i1 > 0) then
    i2 = max(min(i1+1,sizei),1)
  else
    i1 = 1
    i2 = 1
  endif
  if (j1 > 0) then
    j2 = max(min(j1+1,sizej),1)
  else
    j1 = 1
    j2 = 1
  endif
  
  if (j2 == j1) then
    if ( j1 == 1 ) then
      j1_grad = 1
      j2_grad = 2
    else if ( j1 == sizej ) then
      j1_grad = sizej -1
      j2_grad = sizej
    end if
  end if

  !compute i indices for interpolation/extrapolation
  if (i2 == i1) then
    !perform extrapolation in the x-direction
    if ( i1 == 1) then
      i1_grad = 1
      i2_grad = 1 + 1
    else if ( i1 == sizei ) then
      i1_grad = sizei - 1
      i2_grad = sizei
    end if
    !extrapolation in the x-direction
    x1_grad = axis1%values(i1_grad)
    x2_grad = axis1%values(i2_grad)
    x_frac = (lookup1 - x1_grad) / ( x2_grad - x1_grad )    
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then
          call VarUniformXExtrapolate(this,var_idx,i1,i1_grad,i2_grad,sizei, &
                                      j1,j1_grad,j2,j2_grad,&
                                      x1_grad,x2_grad,x_frac,lookup1,lookup2)        
        end if  
      end do
    else
      call VarUniformXExtrapolate(this,var_iname,i1,i1_grad,i2_grad,sizei, &
                                  j1,j1_grad,j2,j2_grad,&
                                  x1_grad,x2_grad,x_frac,lookup1,lookup2)      
    end if    
  else !away from end point in the x-direction
    !inteprolation in the x-direction
    x1 = axis1%values(i1)
    x2 = axis1%values(i2)
    x_frac = (lookup1 - x1) / ( x2 - x1 )    
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then    
          call VarUniformXInterpolate(this,var_idx,j1,j1_grad,j2,j2_grad, &
                                     i1,i2,sizei,x1,x2,x_frac,lookup1,lookup2)
        end if
      end do
    else
      call VarUniformXInterpolate(this,var_iname,j1,j1_grad,j2,j2_grad, &
                                 i1,i2,sizei,x1,x2,x_frac,lookup1,lookup2)      
    end if                                   
  end if
  
  nullify(axis1)
  nullify(axis2)
                   
end subroutine InterpExtrapGrad2DUniform

! ************************************************************************** !

subroutine VarUniformXExtrapolate(this,var_iname,i_col,i1_grad,i2_grad,sizei, &
                                  j1,j1_grad,j2,j2_grad,&
                                  x1_grad,x2_grad,x_frac_grad,lookup1,lookup2)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/25/18
  !     

  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: i_col,i1_grad,i2_grad,sizei
  PetscInt, intent(in) :: j1,j1_grad,j2,j2_grad
  PetscReal, intent(in) :: x1_grad,x2_grad,x_frac_grad,lookup1,lookup2
  
                           
  if (this%var_array(var_iname)%ptr%extrapolation_itype == &
      VAR_EXTRAP_CONST_GRAD ) then
    ! extrapolate  
    call VarUniformXInterpolate(this,var_iname,j1,j1_grad,j2,j2_grad, &
                                i1_grad,i2_grad,sizei,x1_grad,x2_grad, &
                                x_frac_grad,lookup1,lookup2)
  else if ( this%var_array(var_iname)%ptr%extrapolation_itype == & 
            VAR_EXTRAP_CONST_VAL ) then
    !value and gradient in the edge at x = const
    call UniformYValAndGrad(this,var_iname,j1,j2,j1_grad,j2_grad,i_col,sizei, &
                            lookup2,this%var_array(var_iname)%ptr%sample, &              
                            this%var_array(var_iname)%ptr%sample_grad(2))
    this%var_array(var_iname)%ptr%sample_grad(1) = 0.0d0        
    if (this%var_array(var_iname)%ptr%interp_type == &
                                 VAR_INTERP_X_LINLOG ) then
      this%var_array(var_iname)%ptr%sample = &
           dexp(this%var_array(var_iname)%ptr%sample)
      this%var_array(var_iname)%ptr%sample_grad(2) = &
            this%var_array(var_iname)%ptr%sample_grad(2) * &
            this%var_array(var_iname)%ptr%sample
    end if        
  end if

end subroutine VarUniformXExtrapolate

! ************************************************************************** !

subroutine VarUniformXInterpolate(this,var_iname,j1,j1_grad,j2,j2_grad, &
                                 i1,i2,sizei,x1,x2,x_frac,lookup1,lookup2)
! 
! Author: Paolo Orsini
! Date: 05/17/18
!     
  use Utility_module, only : Interpolate, GradientLinear

  implicit none

  class(lookup_table_uniform_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: j1,j1_grad
  PetscInt, intent(in) :: j2,j2_grad
  PetscInt, intent(in) :: i1,i2,sizei
  PetscReal, intent(in) :: x1,x2,x_frac,lookup1,lookup2

  PetscReal :: val_i1,grad_i1,val_i2,grad_i2

  !left
  call UniformYValAndGrad(this,var_iname,j1,j2,j1_grad,j2_grad,i1, &
                                        sizei,lookup2,val_i1,grad_i1)
  !right
  call UniformYValAndGrad(this,var_iname,j1,j2,j1_grad,j2_grad,i2, &
                                          sizei,lookup2,val_i2,grad_i2)
  call GradientLinear(x2,x1,val_i2,val_i1, &
                      this%var_array(var_iname)%ptr%sample_grad(1))
  call Interpolate(x2,x1,lookup1,val_i2,val_i1, &
                   this%var_array(var_iname)%ptr%sample)
  this%var_array(var_iname)%ptr%sample_grad(2) = &
                       grad_i1 * (1.0 - x_frac) + grad_i2 * x_frac
  if (this%var_array(var_iname)%ptr%interp_type == &
                                      VAR_INTERP_X_LINLOG ) then
    this%var_array(var_iname)%ptr%sample = &
       dexp(this%var_array(var_iname)%ptr%sample)
    this%var_array(var_iname)%ptr%sample_grad(1:2) =  &
       this%var_array(var_iname)%ptr%sample_grad(1:2) * &
       this%var_array(var_iname)%ptr%sample
  end if

end subroutine VarUniformXInterpolate

! ************************************************************************** !

subroutine UniformYValAndGrad(this,var_iname,j1,j2,j1_grad,j2_grad,i_col, &
                              sizei,lookup2,val,grad_val)                              
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/25/18
  !   
  ! Given a data set of point in y this routine computes either:
  ! 1) z and dz_dy
  !    OR
  ! 2) ln(z) and d(ln(z))/dy
  !
  ! depending on the interpolation method chosen for the variable
  !
  use Utility_module, only : Interpolate, GradientLinear

  implicit none
  
  class(lookup_table_uniform_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: j1,j2
  PetscInt, intent(in) :: j1_grad,j2_grad
  PetscInt, intent(in) :: i_col, sizei
  PetscReal,intent(in) :: lookup2
  PetscReal, intent(out) :: val, grad_val
  
  PetscReal :: y1, y2, z1, z2
  PetscReal :: y1_grad, y2_grad
  
  
  if (j2 == j1) then
    if ( this%var_array(var_iname)%ptr%extrapolation_itype == &
         VAR_EXTRAP_CONST_GRAD ) then
      y1_grad = this%axis2%values(j1_grad)
      y2_grad = this%axis2%values(j2_grad)
      z1 = this%var_array(var_iname)%ptr%data(i_col+(j1_grad-1)*sizei)
      z2 = this%var_array(var_iname)%ptr%data(i_col+(j2_grad-1)*sizei)
      call GradientLinear(y2_grad,y1_grad,z2,z1,grad_val)
      call Interpolate(y2_grad,y1_grad,lookup2,z2,z1,val)
       if ( this%var_array(var_iname)%ptr%interp_type == &
            VAR_INTERP_X_LINLOG ) then
          grad_val = ( 1.0 / val ) * grad_val
          val = dlog(val)
       end if    
    else if ( this%var_array(var_iname)%ptr%extrapolation_itype == &
              VAR_EXTRAP_CONST_VAL ) then
      val = this%var_array(var_iname)%ptr%data(i_col+(j1-1)*sizei)
      if ( this%var_array(var_iname)%ptr%interp_type == &
           VAR_INTERP_X_LINLOG ) then
        val = dlog(val)
      end if              
      grad_val = 0.0
    end if  
  else !away from end points
    y1 = this%axis2%values(j1)
    y2 = this%axis2%values(j2)
    z1 = this%var_array(var_iname)%ptr%data(i_col+(j1-1)*sizei)
    z2 = this%var_array(var_iname)%ptr%data(i_col+(j2-1)*sizei)
    call Interpolate(y2,y1,lookup2,z2,z1,val)
    call GradientLinear(y2,y1,z2,z1,grad_val)
     if ( this%var_array(var_iname)%ptr%interp_type == &
          VAR_INTERP_X_LINLOG ) then
       grad_val = 1.0 / val * grad_val
       val = dlog(val)
     end if    
  endif
  
end subroutine UniformYValAndGrad

! ************************************************************************** !

subroutine LookupTableInterpolate2DGeneral(this,lookup1,lookup2,result)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module, only : Interpolate

  implicit none
  
  class(lookup_table_general_type) :: this
  PetscInt :: i1a
  PetscInt :: i1b
  PetscInt :: ja
  PetscReal :: lookup1
  PetscReal :: lookup2
  PetscReal :: result
  
  PetscInt :: i2a, i2b
  PetscInt :: jb
  PetscInt :: ii1, ii2, iia, iib
  PetscReal :: x1, x2, z1, z2, xa, xb
  PetscReal :: interp_a, interp_b
  PetscInt :: sizei, sizej

  
  !         x2,y2,z4
  !           /|
  !          / |
  !         /  |
  !        /   |
  !       /    |
  !      /     |
  !  x1,y2,z3  |
  !     |      |
  !     |  x,y |
  !     |      |
  !  x1,y1,z1 -x2,y1,z2                   
 
  result = UNINITIALIZED_DOUBLE
  sizei = this%dims(2)
  sizej = this%dims(1)
  ! index axes
  ja = this%axis1%saved_index
  i1a = this%axis2%saved_index
  i1b = this%axis2%saved_index2
  if (ja > 0) then
    jb = max(min(ja+1,sizej),1)
  else
    ja = 1
    jb = 1
  endif
  if (i1a > 0) then
    i2a = max(min(i1a+1,sizei),1)
  else
    i1a = 1
    i2a = 1
  endif
  if (i1b > 0) then
    i2b = max(min(i1b+1,sizei),1)
  else
    i1b = 1
    i2b = 1
  endif
  if (jb == ja) then
    ! only use ja/i1a/i2a
    if (i2a == i1a) then
      ! corner of domain
      result = this%data(i1a+(ja-1)*sizei)
    else
      ii1 = i1a+(ja-1)*sizei
      ii2 = i2a+(ja-1)*sizei
      x1 = this%axis2%values(ii1)
      x2 = this%axis2%values(ii2)
      z1 = this%data(ii1)
      z2 = this%data(ii2)
      call Interpolate(x2,x1,lookup2,z2,z1,result)
    endif
  else
    ! ja / i*a
    if (i2a == i1a) then
      interp_a = this%data(i1a+(ja-1)*sizei)
    else
      iia = i1a+(ja-1)*sizei
      iib = i2a+(ja-1)*sizei
      x1 = this%axis2%values(iia)
      x2 = this%axis2%values(iib)
      z1 = this%data(iia)
      z2 = this%data(iib)
      call Interpolate(x2,x1,lookup2,z2,z1,interp_a)
    endif
    ! jb / i*b
    if (i2b == i1b) then
      interp_b = this%data(i1b+(jb-1)*sizei)
    else
      iia = i1b+(jb-1)*sizei
      iib = i2b+(jb-1)*sizei
      x1 = this%axis2%values(iia)
      x2 = this%axis2%values(iib)
      z1 = this%data(iia)
      z2 = this%data(iib)
      call Interpolate(x2,x1,lookup2,z2,z1,interp_b)
    endif
    xa = this%axis1%values(ja)
    xb = this%axis1%values(jb)
    call Interpolate(xb,xa,lookup1,interp_b,interp_a,result)
  endif
                   
end subroutine LookupTableInterpolate2DGeneral

! ************************************************************************** !

subroutine LookupTableInterpolate3DLP(this, lookup1, lookup2, lookup3, result)
  !
  ! Interpolation in three dimensions via Lagrange polynomials
  !
  ! Author: Alex Salazar III
  ! Date: 02/22/2022
  !
  implicit none
  ! ----------------------------------
  class(lookup_table_general_type) :: this   ! lookup table
  PetscReal, intent(in) :: lookup1           ! lookup1 value
  PetscReal, intent(in) :: lookup2           ! lookup2 value
  PetscReal, intent(in) :: lookup3           ! lookup3 value
  PetscReal, intent(out) :: result           ! interpolated value
  ! ----------------------------------
  PetscReal, pointer :: x(:)  ! lookup1 table
  PetscReal, pointer :: y(:)  ! lookup2 table
  PetscReal, pointer :: z0(:) ! lookup3 table - full
  type(irregular_array_type) :: za(4)  ! axis3 of interest, needed for partitions
  PetscReal :: F(8)  ! reference data
  PetscInt  :: i, j, k, m
  PetscInt  :: i1, i2
  PetscInt  :: j1, j2
  PetscInt  :: k1, k2
  PetscInt  :: pselect
  PetscInt  :: k11, k12, k13, k14, k21, k22, k23, k24
  PetscReal :: a1, a2, a3, a4, a51, a52, a53, a54, a61, a62, a63, a64
  PetscInt :: indices1(4,2) ! first set of indices used for axis3 indexing
  PetscInt :: indices2(8,3) ! second set of indices used for data indexing
  ! ----------------------------------
  !  METHOD
  !
  !  The lookup value (x, y, z) is bounded from x1 to x2, y1 to y2, and z1 to z2
  !
  !  Eight points form the basis of polynomial interpolation
  !    p1: (x1,y1,z1,u1)
  !    p2: (x1,y2,z1,u2)
  !    p3: (x2,y1,z1,u3)
  !    p4: (x2,y2,z1,u4)
  !    p5: (x1,y1,z2,u5)
  !    p6: (x1,y2,z2,u6)
  !    p7: (x2,y1,z2,u7)
  !    p8: (x2,y2,z2,u8)
  !
  !  Polynomial P(x,y,z) is defined that has roots at each point
  !    P(x,y,z) = (x - x1)(x - x2)(y - y1)(y - y2)(z - z1)(z - z2)
  !
  !  Constituent polynomials fi(x,y,z) are defined such that
  !    f1(x,y,z) = P(x,y,z)/(x - x1)(y - y1)(z - z1)
  !    f2(x,y,z) = P(x,y,z)/(x - x1)(y - y2)(z - z1)
  !    f3(x,y,z) = P(x,y,z)/(x - x2)(y - y1)(z - z1)
  !    f4(x,y,z) = P(x,y,z)/(x - x2)(y - y2)(z - z1)
  !    f5(x,y,z) = P(x,y,z)/(x - x1)(y - y1)(z - z2)
  !    f6(x,y,z) = P(x,y,z)/(x - x1)(y - y2)(z - z2)
  !    f7(x,y,z) = P(x,y,z)/(x - x2)(y - y1)(z - z2)
  !    f8(x,y,z) = P(x,y,z)/(x - x2)(y - y2)(z - z2)
  !
  !  Function U(x,y,z) defined that uses constituent polynomials to go through
  !    all eight points
  !    U(x,y,z) = u1*f1(x,y,z)/f1(x1,y1,z1) +
  !               u2*f2(x,y,z)/f2(x1,y2,z1) +
  !               u3*f3(x,y,z)/f3(x2,y1,z1) +
  !               u4*f4(x,y,z)/f4(x2,y2,z1) +
  !               u5*f5(x,y,z)/f5(x1,y1,z2) +
  !               u6*f6(x,y,z)/f6(x1,y2,z2) +
  !               u7*f7(x,y,z)/f7(x2,y1,z2) +
  !               u8*f8(x,y,z)/f8(x2,y2,z2)

  x  => this%axis1%values ! axis1 values
  y  => this%axis2%values ! axis2 values
  z0 => this%axis3%values ! axis3 values (full)

  ! retrieve i and j indices (pivot variables)
  i1 = this%axis1%saved_index
  i2 = this%axis1%saved_index1
  j1 = this%axis2%saved_index
  j2 = this%axis2%saved_index2

  ! (i, j) combinations
  indices1(1,:) = (/i1, j1/) ! a*1 coefficients
  indices1(2,:) = (/i1, j2/) ! a*2 coefficients
  indices1(3,:) = (/i2, j1/) ! a*3 coefficients
  indices1(4,:) = (/i2, j2/) ! a*4 coefficients

  ! retrieve k indices (independent variables)
  do m = 1, 4
    i = indices1(m,1)
    j = indices1(m,2)
    k1 = this%axis3%saved_indices3(i, j, 1)
    k2 = this%axis3%saved_indices3(i, j, 2)
    indices2(m,:)   = (/i, j, k1/) ! a5* coefficients
    indices2(m+4,:) = (/i, j, k2/) ! a6* coefficients
    ! redefine z-axis as-needed
    if (allocated(this%axis3%partition)) then
      pselect = this%axis3%saved_index_partition(i, j)
      za(m)%data => this%axis3%partition(pselect)%data
    else
      za(m)%data => this%axis3%values
    endif
  enddo

  ! retrieve data
  do m = 1, 8
    i = indices2(m,1)
    j = indices2(m,2)
    k = indices2(m,3)
    F(m) = this%data_references(i, j, k)
  enddo

  ! define coefficients of polymomial based on lookups
  !   lookups are compared to opposite i, j, and k indices
  a1  = (lookup1 - x(i2)) / ( x(i1) - x(i2) ) ! a1 -> i1
  a2  = (lookup1 - x(i1)) / ( x(i2) - x(i1) ) ! a2 -> i2
  a3  = (lookup2 - y(j2)) / ( y(j1) - y(j2) ) ! a3 -> j1
  a4  = (lookup2 - y(j1)) / ( y(j2) - y(j1) ) ! a4 -> j2

  k11 = this%axis3%saved_indices3(i1, j1, 1) ! k11 -> i1,j1,k1
  k12 = this%axis3%saved_indices3(i1, j2, 1) ! k12 -> i1,j2,k1
  k13 = this%axis3%saved_indices3(i2, j1, 1) ! k13 -> i2,j1,k1
  k14 = this%axis3%saved_indices3(i2, j2, 1) ! k14 -> i2,j2,k1
  k21 = this%axis3%saved_indices3(i1, j1, 2) ! k21 -> i1,j1,k2
  k22 = this%axis3%saved_indices3(i1, j2, 2) ! k22 -> i1,j2,k2
  k23 = this%axis3%saved_indices3(i2, j1, 2) ! k23 -> i2,j1,k2
  k24 = this%axis3%saved_indices3(i2, j2, 2) ! k24 -> i2,j2,k2

  a51 = (lookup3 - za(1)%data(k21)) / (za(1)%data(k11) - za(1)%data(k21)) ! a51 -> i1,j1,k1
  a52 = (lookup3 - za(2)%data(k22)) / (za(2)%data(k12) - za(2)%data(k22)) ! a52 -> i1,j2,k1
  a53 = (lookup3 - za(3)%data(k23)) / (za(3)%data(k13) - za(3)%data(k23)) ! a53 -> i2,j1,k1
  a54 = (lookup3 - za(4)%data(k24)) / (za(4)%data(k14) - za(4)%data(k24)) ! a54 -> i2,j2,k1
  a61 = (lookup3 - za(1)%data(k11)) / (za(1)%data(k21) - za(1)%data(k11)) ! a61 -> i1,j1,k2
  a62 = (lookup3 - za(2)%data(k12)) / (za(2)%data(k22) - za(2)%data(k12)) ! a62 -> i1,j2,k2
  a63 = (lookup3 - za(3)%data(k13)) / (za(3)%data(k23) - za(3)%data(k13)) ! a63 -> i2,j1,k2
  a64 = (lookup3 - za(4)%data(k14)) / (za(4)%data(k24) - za(4)%data(k14)) ! a64 -> i2,j2,k2
  
  ! product operator bypasses identical indices --> pick one term to cancel
  if (i1 == i2) then
    a1 = 1
    a2 = 0
  endif
  if (j1 == j2) then
    a3 = 1
    a4 = 0
  endif
  if (k11 == k21) then
    a51 = 1
    a61 = 0
  endif
  if (k12 == k22) then
    a52 = 1
    a62 = 0
  endif
  if (k13 == k23) then
    a53 = 1
    a63 = 0
  endif
  if (k14 == k24) then
    a54 = 1
    a64 = 0
  endif

  ! obtain result from interpolating polynomial
  result = a1*a3*a51*F(1) + & ! F and a51 -> i1,j1,k1 ! a1 -> i1 ! a3 -> j1
           a1*a4*a52*F(2) + & ! F and a52 -> i1,j2,k1 ! a1 -> i1 ! a4 -> j2
           a2*a3*a53*F(3) + & ! F and a53 -> i2,j1,k1 ! a2 -> i2 ! a3 -> j1
           a2*a4*a54*F(4) + & ! F and a54 -> i2,j2,k1 ! a2 -> i2 ! a4 -> j2
           a1*a3*a61*F(5) + & ! F and a61 -> i1,j1,k2 ! a1 -> i1 ! a3 -> j1
           a1*a4*a62*F(6) + & ! F and a62 -> i1,j2,k2 ! a1 -> i1 ! a4 -> j2
           a2*a3*a63*F(7) + & ! F and a63 -> i2,j1,k2 ! a2 -> i2 ! a3 -> j1
           a2*a4*a64*F(8)     ! F and a64 -> i2,j2,k2 ! a2 -> i2 ! a4 -> j2

  if (associated(x)) nullify(x)
  if (associated(y)) nullify(y)
  if (associated(z0)) nullify(z0)

end subroutine LookupTableInterpolate3DLP

! ************************************************************************** !

subroutine LookupTableInterpolate3DTrilinear(this, lookup1, lookup2, lookup3, &
                                             result)
  !
  ! Trilinear interpolation of lookup values
  !
  ! Author: Alex Salazar III
  ! Date: 02/22/2022
  !
  implicit none
  ! ----------------------------------
  class(lookup_table_general_type) :: this   ! lookup table
  PetscReal, intent(in) :: lookup1           ! lookup1 value
  PetscReal, intent(in) :: lookup2           ! lookup2 value
  PetscReal, intent(in) :: lookup3           ! lookup3 value
  PetscReal, intent(out) :: result ! interpolated value
  ! ----------------------------------
  PetscReal, pointer :: x(:)  ! lookup1 table
  PetscReal, pointer :: y(:)  ! lookup2 table
  PetscReal, pointer :: z0(:) ! lookup3 table - full
  type(irregular_array_type) :: za(4)  ! axis3 of interest, needed for partitions
  PetscReal :: F(8)  ! reference data
  PetscInt  :: i, j, k, m ! iterators
  PetscInt  :: i1, i2
  PetscInt  :: j1, j2
  PetscInt  :: k1, k2
  PetscInt  :: pselect
  PetscInt :: indices1(4,2) ! first set of indices used for axis3 indexing
  PetscInt :: indices2(8,3) ! second set of indices used for data indexing
  PetscReal :: xd, yd, zd   ! axial differences
  PetscReal :: zmin, zmax
  PetscReal :: zlog(8)
  PetscReal :: c000   ! data point cube vertex (x1,y1,z1)
  PetscReal :: c010   ! data point cube vertex (x1,y2,z1)
  PetscReal :: c100   ! data point cube vertex (x2,y1,z1)
  PetscReal :: c110   ! data point cube vertex (x2,y2,z1)
  PetscReal :: c001   ! data point cube vertex (x1,y1,z2)
  PetscReal :: c011   ! data point cube vertex (x1,y2,z2)
  PetscReal :: c101   ! data point cube vertex (x2,y1,z2)
  PetscReal :: c111   ! data point cube vertex (x2,y2,z2)
  PetscReal :: c00    ! data point plane vertex
  PetscReal :: c01    ! data point plane vertex
  PetscReal :: c10    ! data point plane vertex
  PetscReal :: c11    ! data point plane vertex
  PetscReal :: c0     ! data point line vertex
  PetscReal :: c1     ! data point line vertex
  ! ----------------------------------
  !  METHOD
  !
  !  The lookup value (x, y, z) is bounded from x1 to x2, y1 to y2, and z1 to z2
  !
  !  Eight points form vertices of cube around the lookup value
  !    p1: (x1,y1,z1,u1)
  !    p2: (x1,y2,z1,u2)
  !    p3: (x2,y1,z1,u3)
  !    p4: (x2,y2,z1,u4)
  !    p5: (x1,y1,z2,u5)
  !    p6: (x1,y2,z2,u6)
  !    p7: (x2,y1,z2,u7)
  !    p8: (x2,y2,z2,u8)
  !
  !  Determine differences along x, y, and z between the lookup values and
  !    vertices, where the minimum-bounding vertex serves as the origin.
  !
  !  Marching along x direction beginning at origin, find intermediate
  !    interpolated data values using data from all vertices
  !    -- this forms a yz plane.
  !
  !  Marching along y direction on the yz plane defined earlier, find
  !    intermediate interpolated data values using four points of the yz plane
  !    -- this forms a line parallel with the z axis.
  !
  !  Marching along the z line defined earlier, find final interpolated data
  !    value.

  x  => this%axis1%values ! axis1 values
  y  => this%axis2%values ! axis2 values
  z0 => this%axis3%values ! axis3 values (full)

  ! retrieve i and j indices (pivot variables)
  i1 = this%axis1%saved_index
  i2 = this%axis1%saved_index1
  j1 = this%axis2%saved_index
  j2 = this%axis2%saved_index2

  ! x and y differences
  xd = (lookup1 - x(i1))/(x(i2) - x(i1))
  yd = (lookup2 - y(j1))/(y(j2) - y(j1))

  if (i1 == i2) then
    xd = 1
  endif
  if (j1 == j2) then
    yd = 1
  endif

  ! (i, j) combinations
  indices1(1,:) = (/i1, j1/) ! x1, y1 per z
  indices1(2,:) = (/i1, j2/) ! x1, y2 per z
  indices1(3,:) = (/i2, j1/) ! x2, y1 per z
  indices1(4,:) = (/i2, j2/) ! x2, y2 per z

  ! retrieve k indices (independent variables)
  do m = 1, 4
    i = indices1(m,1)
    j = indices1(m,2)
    k1 = this%axis3%saved_indices3(i, j, 1)
    k2 = this%axis3%saved_indices3(i, j, 2)
    indices2(m,:)   = (/i, j, k1/) ! z1
    indices2(m+4,:) = (/i, j, k2/) ! z2
    ! redefine z-axis as-needed
    if (allocated(this%axis3%partition)) then
      pselect = this%axis3%saved_index_partition(i, j)
      za(m)%data => this%axis3%partition(pselect)%data
    else
      za(m)%data => this%axis3%values
    endif
    zlog(m)   = za(m)%data(k1)
    zlog(m+4) = za(m)%data(k2)
  enddo

  ! retrieve data
  do m = 1, 8
    i = indices2(m,1)
    j = indices2(m,2)
    k = indices2(m,3)
    F(m) = this%data_references(i, j, k)
  enddo

  ! z difference
  zmin = minval(zlog)
  zmax = maxval(zlog)
  zd = (lookup3 - zmin)/(zmax - zmin)
  if (zmax == zmin) then
    zd = 1
  endif

  ! define cube vertex data values
  c000 = F(1)
  c010 = F(2)
  c100 = F(3)
  c110 = F(4)
  c001 = F(5)
  c011 = F(6)
  c101 = F(7)
  c111 = F(8)

  ! interpolate along x-axis
  c00 = c000*(1 - xd) + c100*xd
  c01 = c001*(1 - xd) + c101*xd
  c10 = c010*(1 - xd) + c110*xd
  c11 = c011*(1 - xd) + c111*xd

  ! interpolate along y-axis
  c0 = c00*(1 - yd) + c10*yd
  c1 = c01*(1 - yd) + c11*yd

  ! interpolate along z-axis
  result = c0*(1- zd) + c1*zd
  
  if (associated(x)) nullify(x)
  if (associated(y)) nullify(y)
  if (associated(z0)) nullify(z0)

end subroutine LookupTableInterpolate3DTrilinear

! ************************************************************************** !

subroutine InterpExtrapGradGeneral2D(this,var_iname,lookup1,lookup2)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/17/18
  ! 
  use Utility_module, only : Interpolate, GradientLinear

  implicit none
  
  class(lookup_table_general_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscReal, intent(in) :: lookup1
  PetscReal, intent(in) :: lookup2
  
  PetscInt :: var_idx
  PetscInt :: ja
  PetscInt :: i1a, i1b
  PetscInt :: jb
  PetscInt :: i2a, i2b
  PetscReal :: xa, xb
  PetscInt :: sizei, sizej
  PetscInt :: ja_grad, jb_grad
  PetscInt :: i1a_grad, i2a_grad, i1b_grad, i2b_grad
  PetscReal :: xa_grad, xb_grad, x_frac

  
  !         x2,y2,z4
  !           /|
  !          / |
  !         /  |
  !        /   |
  !       /    |
  !      /     |
  !  x1,y2,z3  |
  !     |      |
  !     |  x,y |
  !     |      |
  !  x1,y1,z1 -x2,y1,z2   
  !
  !    ja     jb
  ! note that j refer to x
 
  ! x used for both i and j to indicate all operations are 1D
  sizei = this%dims(2)  ! i referes to y 
  sizej = this%dims(1)  ! j refers to x
  x_frac = 0.0
  ! index axes
  ja = this%axis1%saved_index
  i1a = this%axis2%saved_index
  i1b = this%axis2%saved_index2
  if (ja > 0) then
    jb = max(min(ja+1,sizej),1)
  else
    ja = 1
    jb = 1
  endif
  if (i1a > 0) then
    i2a = max(min(i1a+1,sizei),1)
  else
    i1a = 1
    i2a = 1
  endif
  if (i1b > 0) then
    i2b = max(min(i1b+1,sizei),1)
  else
    i1b = 1
    i2b = 1
  endif
  
  if ( i1a == i2a ) then
    if ( i1a == 1 ) then
      i1a_grad = 1
      i2a_grad = 2
    else if (i1a == sizei) then
      i1a_grad = sizei - 1
      i2a_grad = sizei
    end if  
  end if  

  if ( i1b == i2b ) then
    if ( i1b == 1 ) then
      i1b_grad = 1
      i2b_grad = 2
    else if (i1b == sizei) then
      i1b_grad = sizei - 1
      i2b_grad = sizei
    end if  
  end if
  
  if (jb == ja) then
  !left and right table edges (included the corners that are extrapolated)
    if ( ja == 1 ) then
      ja_grad = 1
      jb_grad = 2
      !new look up to find i1b_grad, i2b_grad 
      call GeneralFindIndicesAxis2(this,jb_grad,lookup2,i1b,i2b, &
                                   i1b_grad,i2b_grad)
    else if ( ja == sizej ) then
      ja_grad = sizej - 1
      jb_grad = sizej
      !new look up to find i1a_grad, i2a_grad 
      call GeneralFindIndicesAxis2(this,ja_grad,lookup2,i1b,i2b, &
                                   i1a_grad,i2a_grad)
    end if
    xa_grad = this%axis1%values(ja_grad)
    xb_grad = this%axis1%values(jb_grad)
    x_frac = (lookup1 - xa_grad) / (xb_grad - xa_grad)
    !interpolation/extrapolation left and right columns             
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then
          call VarGeneralXExtrapolate(this,var_idx,ja,i1a,i2a,sizei, &
                                      ja_grad,i1a_grad,i2a_grad, &
                                      jb_grad,i1b_grad,i2b_grad, &
                                      xa_grad,xb_grad,x_frac,lookup1,lookup2)
        end if
      end do                             
    else
      ! to call function var_extrapolate function
      call VarGeneralXExtrapolate(this,var_iname,ja,i1a,i2a,sizei, &
                                        ja_grad,i1a_grad,i2a_grad, &
                                        jb_grad,i1b_grad,i2b_grad, &
                                        xa_grad,xb_grad,x_frac,lookup1,lookup2) 
    end if
  else !away from left/right table edges - no extrapolation in x directions
    xa = this%axis1%values(ja)
    xb = this%axis1%values(jb)
    x_frac = ( lookup1 - xa ) / ( xb - xa )
    if (var_iname == -1 ) then
      do var_idx = 1,size(this%var_array(:))
        if ( associated(this%var_array(var_idx)%ptr) ) then
          call VarGeneralXInterpolate(this,var_idx,ja,i1a,i2a, &
                     i1a_grad,i2a_grad,sizei,jb,i1b,i2b,i1b_grad,i2b_grad, &
                     xa,xb,x_frac,lookup1,lookup2) 
        end if
      end do          
    else
      call VarGeneralXInterpolate(this,var_iname,ja,i1a,i2a, &
                 i1a_grad,i2a_grad,sizei,jb,i1b,i2b,i1b_grad,i2b_grad, &
                 xa,xb,x_frac,lookup1,lookup2)
    end if
  end if        
                   
end subroutine InterpExtrapGradGeneral2D

! ************************************************************************** !

subroutine VarGeneralXExtrapolate(this,var_iname,ja,i1a,i2a,sizei, &
                                  ja_grad,i1a_grad,i2a_grad, &
                                  jb_grad,i1b_grad,i2b_grad, &
                                  xa_grad,xb_grad,x_frac,lookup1,lookup2)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/17/18
  !     
  use Utility_module, only : Interpolate, GradientLinear

  implicit none
  
  class(lookup_table_general_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: ja,i1a,i2a,sizei
  PetscInt, intent(in) :: ja_grad,i1a_grad,i2a_grad
  PetscInt, intent(in) :: jb_grad,i1b_grad,i2b_grad
  PetscReal, intent(in) :: xa_grad,xb_grad,x_frac,lookup1,lookup2
  
  !PetsReal :: val_a,val_b,grad_a,grad_b
                           
  if (this%var_array(var_iname)%ptr%extrapolation_itype == &
      VAR_EXTRAP_CONST_GRAD ) then

    call VarGeneralXInterpolate(this,var_iname,ja_grad,i1a,i2a,&
                  i1a_grad,i2a_grad,sizei,jb_grad,i1a,i2a,i1b_grad,i2b_grad, &
                              xa_grad,xb_grad,x_frac,lookup1,lookup2)
  else if ( this%var_array(var_iname)%ptr%extrapolation_itype == & 
            VAR_EXTRAP_CONST_VAL ) then
    !value and gradient in the edge at x = const
    call GeneralYValAndGrad(this,var_iname,ja,i1a,i2a,i1a_grad,i2a_grad, &
                             sizei,lookup2, &
                             this%var_array(var_iname)%ptr%sample, &
                             this%var_array(var_iname)%ptr%sample_grad(2))
    this%var_array(var_iname)%ptr%sample_grad(1) = 0.0d0                         
    if (this%var_array(var_iname)%ptr%interp_type == &
                                 VAR_INTERP_X_LINLOG ) then
      this%var_array(var_iname)%ptr%sample = &
           dexp(this%var_array(var_iname)%ptr%sample)
      this%var_array(var_iname)%ptr%sample_grad(2) = &
            this%var_array(var_iname)%ptr%sample_grad(2) * &
            this%var_array(var_iname)%ptr%sample
    end if        
  end if

end subroutine VarGeneralXExtrapolate

! ************************************************************************** !

subroutine VarGeneralXInterpolate(this,var_iname,ja,i1a,i2a,i1a_grad,i2a_grad, &
                                 sizei,jb,i1b,i2b,i1b_grad,i2b_grad,&
                                 xa,xb,x_frac,lookup1,lookup2)
! 
! Author: Paolo Orsini
! Date: 05/17/18
!     
  use Utility_module, only : Interpolate, GradientLinear

  implicit none

  class(lookup_table_general_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: ja,i1a,i2a,i1a_grad,i2a_grad,sizei
  PetscInt, intent(in) :: jb,i1b,i2b,i1b_grad,i2b_grad
  PetscReal, intent(in) :: xa,xb,x_frac,lookup1,lookup2

  PetscReal :: val_a,grad_a,val_b,grad_b

  !left
  call GeneralYValAndGrad(this,var_iname,ja,i1a,i2a,i1a_grad,i2a_grad, &
                          sizei,lookup2,val_a,grad_a)
  !right  
  call GeneralYValAndGrad(this,var_iname,jb,i1b,i2b,i1b_grad,i2b_grad, &
                          sizei,lookup2,val_b,grad_b)
  call GradientLinear(xb,xa,val_b,val_a, &
                      this%var_array(var_iname)%ptr%sample_grad(1))
  call Interpolate(xb,xa,lookup1,val_b,val_a, &
                   this%var_array(var_iname)%ptr%sample)
  this%var_array(var_iname)%ptr%sample_grad(2) = &
                       grad_a * (1.0 - x_frac) + grad_b * x_frac
  if (this%var_array(var_iname)%ptr%interp_type == &
                                      VAR_INTERP_X_LINLOG ) then
    this%var_array(var_iname)%ptr%sample = &
       dexp(this%var_array(var_iname)%ptr%sample)
    this%var_array(var_iname)%ptr%sample_grad(1:2) =  &
       this%var_array(var_iname)%ptr%sample_grad(1:2) * &
       this%var_array(var_iname)%ptr%sample
  end if

end subroutine VarGeneralXInterpolate

! ************************************************************************** !

subroutine GeneralYValAndGrad(this,var_iname,j_col,i1,i2,i1_grad,i2_grad, &
                              sizei,lookup,val,grad_val)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/17/18
  !   
  ! Given a data set of point in y this routine computes either:
  ! 1) z and dz_dy
  !    OR
  ! 2) ln(z) and d(ln(z))/dy
  !
  ! depending on the interpolation method chosen for the variable
  !
  use Utility_module, only : Interpolate, GradientLinear

  implicit none
  
  class(lookup_table_general_type) :: this
  PetscInt, intent(in) :: var_iname
  PetscInt, intent(in) :: j_col
  PetscInt, intent(in) :: i1,i2
  PetscInt, intent(in) :: i1_grad, i2_grad
  PetscInt, intent(in) :: sizei
  PetscReal,intent(in) :: lookup
  PetscReal, intent(out) :: val, grad_val
  
  !PetscInt :: i1_grad, i2_grad
  PetscInt :: iia, iib
  PetscReal :: x1, x2, z1, z2
  PetscReal :: x1_grad, x2_grad

  
  if (i2 == i1) then
    if ( this%var_array(var_iname)%ptr%extrapolation_itype == &
         VAR_EXTRAP_CONST_GRAD ) then
      iia = i1+(j_col-1)*sizei
      x1 = this%axis2%values(iia)
      iia = i1_grad+(j_col-1)*sizei
      iib = i2_grad+(j_col-1)*sizei
      x1_grad = this%axis2%values(iia)
      x2_grad = this%axis2%values(iib)
      z1 = this%var_array(var_iname)%ptr%data(iia)
      z2 = this%var_array(var_iname)%ptr%data(iib)
      !log numerical derivatives for testing
      ! if ( this%var_array(var_iname)%ptr%interp_type == &
      !      VAR_INTERP_X_LINLOG ) then
      !   z1 = dlog(z1)
      !   z2 = dlog(z2)
      ! end if     
      call GradientLinear(x2_grad,x1_grad,z2,z1,grad_val)
      call Interpolate(x2_grad,x1_grad,lookup,z2,z1,val)
      !log analytical derivatives
      if ( this%var_array(var_iname)%ptr%interp_type == &
           VAR_INTERP_X_LINLOG ) then
        grad_val = ( 1.0 / val) * grad_val
        val = dlog(val)
      end if           
    else if ( this%var_array(var_iname)%ptr%extrapolation_itype == &
              VAR_EXTRAP_CONST_VAL ) then
      val = this%var_array(var_iname)%ptr%data(i1+(j_col-1)*sizei)
      if ( this%var_array(var_iname)%ptr%interp_type == &
           VAR_INTERP_X_LINLOG ) then
        val = dlog(val)   
      end if           
      grad_val = 0.0
    end if  
  else !away from end points
    iia = i1+(j_col-1)*sizei
    iib = i2+(j_col-1)*sizei
    x1 = this%axis2%values(iia)
    x2 = this%axis2%values(iib)
    z1 = this%var_array(var_iname)%ptr%data(iia)
    z2 = this%var_array(var_iname)%ptr%data(iib)
    !log numerical derivatives for testing
    ! if ( this%var_array(var_iname)%ptr%interp_type == &
    !      VAR_INTERP_X_LINLOG ) then
    !   z1 = dlog(z1)
    !   z2 = dlog(z2)
    ! end if         
    call Interpolate(x2,x1,lookup,z2,z1,val)
    call GradientLinear(x2,x1,z2,z1,grad_val)
    if ( this%var_array(var_iname)%ptr%interp_type == &
          VAR_INTERP_X_LINLOG ) then
      grad_val = 1.0 / val * grad_val 
      val = dlog(val)
    end if      
  endif
  
end subroutine GeneralYValAndGrad

! ************************************************************************** !

subroutine GeneralFindIndicesAxis2(this,j_col,lookup,i1,i2,i1_grad,i2_grad)
  ! 
  ! Author: Paolo Orsini
  ! Date: 05/21/18
  !   
  implicit none  
  
  class(lookup_table_general_type) :: this
  PetscInt, intent(in)  :: j_col
  PetscReal, intent(in)  :: lookup
  PetscInt, intent(out) :: i1
  PetscInt, intent(out) :: i2
  PetscInt, intent(out) :: i1_grad
  PetscInt, intent(out) :: i2_grad

  PetscInt :: iend, istart, sizei
  
  sizei = this%dims(2)
  iend = j_col*this%dims(2)
  istart = iend - this%dims(2) + 1
  call LookupTableAxisIndexGeneral(lookup,this%axis2%values(istart:iend),i1)
  if (i1 > 0) then
    i2 = max(min(i1+1,sizei),1)
  else
    i1 = 1
    i2 = 1
  end if
  
  if ( i1 == i2 ) then
    if ( i1 == 1 ) then
      i1_grad = 1
      i2_grad = 2
    else if ( i1 == sizei) then
      i1_grad = sizei - 1
      i2_grad = sizei
    end if    
  end if  
  
end subroutine GeneralFindIndicesAxis2

! ************************************************************************** !

subroutine LookupTableTest1D(lookup_table,lookup,desired_result)
  ! 
  ! Deallocates any allocated pointers in axis
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_base_type) :: lookup_table
  PetscReal :: lookup
  PetscReal :: desired_result
  
  PetscReal :: result
  
100 format(3(f12.6),l)

  result = lookup_table%Sample(lookup)
  write(*,100) lookup, desired_result, result, Equal(result,desired_result)

end subroutine LookupTableTest1D

! ************************************************************************** !

subroutine LookupTableTest2D(lookup_table,lookup1,lookup2,desired_result)
  ! 
  ! Deallocates any allocated pointers in axis
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_base_type) :: lookup_table
  PetscReal :: lookup1
  PetscReal :: lookup2
  PetscReal :: desired_result
  
  PetscReal :: result
  
100 format(4(f12.6),l)

  result = lookup_table%Sample(lookup1,lookup2)
  write(*,100) lookup1, lookup2, desired_result, result, &
    Equal(result,desired_result)

end subroutine LookupTableTest2D

! ************************************************************************** !

subroutine LookupTableAxisDestroy(axis)
  ! 
  ! Deallocates any allocated pointers in axis
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_axis_type), pointer :: axis
  
  if (.not.associated(axis)) return
  
  call DeallocateArray(axis%values)
  deallocate(axis)
  nullify(axis)

end subroutine LookupTableAxisDestroy

! ************************************************************************** !

subroutine LookupTableBaseDestroy(lookup_table)
  !
  ! Deallocates any allocated pointers in lookup table base type
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/17
  !
  use Utility_module

  implicit none

  class(lookup_table_base_type) :: lookup_table

  PetscInt :: i

  call LookupTableAxisDestroy(lookup_table%axis1)

  if (associated(lookup_table%var_array)) then
    do i =1,size(lookup_table%var_array(:))
      nullify(lookup_table%var_array(i)%ptr)
    end do
    deallocate(lookup_table%var_array)
    nullify(lookup_table%var_array)
  end if

  call LookupTableVarListDestroy(lookup_table%vars)

  !deallocate data arrays at last, as lookup variables point to it
  call DeallocateArray(lookup_table%data)
  call DeallocateArray(lookup_table%var_data)

end subroutine LookupTableBaseDestroy

! ************************************************************************** !

subroutine LookupTableUniformDestroy(lookup_table)
  ! 
  ! Deallocates any allocated pointers in lookup table
  ! 
  ! Author: Glenn Hammond, Paolo Orsini
  ! Date: 10/15/14, 04/18/2018
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_uniform_type), pointer :: lookup_table
  
  call LookupTableBaseDestroy(lookup_table)
  call LookupTableAxisDestroy(lookup_table%axis2)
  call LookupTableAxisDestroy(lookup_table%axis3)
  deallocate(lookup_table)
  nullify(lookup_table)

end subroutine LookupTableUniformDestroy

! ************************************************************************** !

subroutine LookupTableGeneralDestroy(lookup_table)
  ! 
  ! Deallocates any allocated pointers in lookup table
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/15/14
  ! 
  use Utility_module
  
  implicit none
  
  class(lookup_table_general_type), pointer :: lookup_table
  
  if (.not.associated(lookup_table)) return

  call LookupTableBaseDestroy(lookup_table)
  ! axis2 is different type
  if (associated(lookup_table%axis2)) then
    call DeallocateArray(lookup_table%axis2%values)
    deallocate(lookup_table%axis2)
    nullify(lookup_table%axis2)
  endif
  ! axis3 is different type
  if (associated(lookup_table%axis3)) then
    call DeallocateArray(lookup_table%axis3%values)
    deallocate(lookup_table%axis3)
    nullify(lookup_table%axis3)
  endif
  deallocate(lookup_table)
  nullify(lookup_table)

end subroutine LookupTableGeneralDestroy


! ************************************************************************** !
! ** lookup table variable procedures
! ************************************************************************** !

subroutine LookupTableVarsInit(this,n_var_max)
  !
  ! Initializes a lookup table variables
  !
  ! Author: Paolo Orsini
  ! Date: 05/04/18
  !

  implicit none

  class(lookup_table_base_type) :: this
  PetscInt :: n_var_max
  
  PetscInt :: i_var
  
  allocate(this%vars)
  call LookupTableVarInitList(this%vars)
  allocate(this%var_array(n_var_max))
  do i_var = 1, n_var_max
    nullify(this%var_array(i_var)%ptr)
  end do

end subroutine LookupTableVarsInit

! ************************************************************************** !

function LookupTableVarIsPresent(this,var_iname)
  !
  ! Check is a lookup table var is present in an existing lookup table
  !
  ! Author: Paolo Orsini
  ! Date: 05/04/18
  !

  implicit none

  PetscBool :: LookupTableVarIsPresent
  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname

  LookupTableVarIsPresent = associated(this%var_array(var_iname)%ptr)

end function LookupTableVarIsPresent

! ************************************************************************** !

function LookupTableVarIsSMInc(this,var_iname)
  !
  ! Check if a lookup table var is strictly monotonically increasing
  !
  ! Author: Paolo Orsini
  ! Date: 02/11/19
  !

  implicit none

  PetscBool :: LookupTableVarIsSMInc
  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname

  PetscInt :: i_data

  LookupTableVarIsSMInc = PETSC_TRUE

  do i_data = 2, size(this%var_array(var_iname)%ptr%data(:))
    if ( this%var_array(var_iname)%ptr%data(i_data) <= &
         this%var_array(var_iname)%ptr%data(i_data - 1 ) ) then
      LookupTableVarIsSMInc = PETSC_FALSE
    end if
  end do

end function LookupTableVarIsSMInc

! ************************************************************************** !

subroutine LookupTableVarInitList(list)
  !
  ! Initializes a lookup table var list
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/17
  !

  implicit none

  type(lookup_table_var_list_type) :: list

  nullify(list%first)
  nullify(list%last)
  list%num_lookup_table_vars = 0

end subroutine LookupTableVarInitList

! ************************************************************************** !

subroutine LookupTableVarAddToList(new_var,list)
  !
  ! Adds a new lookup table var to a lookup table vars list
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/2017
  !

  implicit none

  type(lookup_table_var_type), pointer :: new_var
  type(lookup_table_var_list_type) :: list

  list%num_lookup_table_vars = list%num_lookup_table_vars + 1
  new_var%id = list%num_lookup_table_vars
  if (.not.associated(list%first)) list%first => new_var
  if (associated(list%last)) list%last%next => new_var
  list%last => new_var

end subroutine LookupTableVarAddToList

! ************************************************************************** !

function CreateLookupTableVar1(var_iname,internal_units,user_units,data_idx)
  !
  ! Creates a new lookup table var
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/2017 - modfied 04/17/2018
  !

  implicit none

  type(lookup_table_var_type), pointer :: CreateLookupTableVar1
  PetscInt, intent(in) :: var_iname
  character(len=MAXWORDLENGTH), intent(in) :: internal_units
  character(len=MAXWORDLENGTH), intent(in) :: user_units
  PetscInt, intent(in) :: data_idx

  allocate(CreateLookupTableVar1)

  CreateLookupTableVar1%id = UNINITIALIZED_INTEGER
  CreateLookupTableVar1%iname = var_iname
  CreateLookupTableVar1%data_idx = data_idx
  CreateLookupTableVar1%extrapolation_itype = UNINITIALIZED_INTEGER
  CreateLookupTableVar1%interp_type = VAR_INTERP_LINEAR
  CreateLookupTableVar1%internal_units = trim(internal_units)
  CreateLookupTableVar1%user_units = trim(user_units)
  CreateLookupTableVar1%conversion_factor = UNINITIALIZED_DOUBLE
  nullify(CreateLookupTableVar1%data)
  CreateLookupTableVar1%sample = UNINITIALIZED_DOUBLE
  nullify(CreateLookupTableVar1%sample_grad)
  nullify(CreateLookupTableVar1%next)

end function CreateLookupTableVar1

! ************************************************************************** !

function CreateLookupTableVar2(lookup_var)
  !
  ! Creates a new lookup table var from an exisitng lookup table var
  !
  ! Author: Paolo Orsini
  ! Date: 08/20/2018
  !

  implicit none

  type(lookup_table_var_type), pointer :: CreateLookupTableVar2
  type(lookup_table_var_type), intent(in) :: lookup_var

  allocate(CreateLookupTableVar2)

  CreateLookupTableVar2%id = lookup_var%id
  CreateLookupTableVar2%iname = lookup_var%iname
  CreateLookupTableVar2%data_idx = lookup_var%data_idx
  CreateLookupTableVar2%extrapolation_itype = lookup_var%extrapolation_itype
  CreateLookupTableVar2%interp_type = lookup_var%interp_type
  CreateLookupTableVar2%internal_units = lookup_var%internal_units
  CreateLookupTableVar2%user_units = lookup_var%user_units
  CreateLookupTableVar2%conversion_factor = lookup_var%conversion_factor
  CreateLookupTableVar2%data => lookup_var%data
  CreateLookupTableVar2%sample = lookup_var%sample
  allocate(CreateLookupTableVar2%sample_grad(size(lookup_var%sample_grad(:))))
  CreateLookupTableVar2%sample_grad = UNINITIALIZED_DOUBLE
  nullify(CreateLookupTableVar2%next)

end function CreateLookupTableVar2

! ************************************************************************** !

subroutine CreateAddLookupTableVar(this,var_iname,internal_units,user_units, &
                                                            data_idx, option)

  ! create and dd a lokup variable to a lookup table
  !
  ! Author: Paolo Orsini
  ! Date: 01/19/2019
  !

  use Option_module

  implicit none

  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname
  character(len=MAXWORDLENGTH), intent(in) :: internal_units
  character(len=MAXWORDLENGTH), intent(in) :: user_units
  PetscInt, intent(in) :: data_idx
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: lookup_var => null()

  if ( .not.associated(this%var_array) ) then
    option%io_buffer = 'CreateAddLookupTableVar: Cannot add a lookup variable &
                        &without intialising the lookup variable list'
    call PrintErrMsg(option)
  end if

  if ( var_iname > size (this%var_array(:)) .or. &
       var_iname <= 0 ) then
    option%io_buffer = 'var_iname must be larger than zero and not larger &
                        &than the maximum number of lookup variables'
    call PrintErrMsg(option)
  end if

  lookup_var => CreateLookupTableVar(var_iname,internal_units, &
                                            user_units,data_idx)
  call LookupTableVarAddToList(lookup_var,this%vars)
  this%var_array(lookup_var%iname)%ptr => lookup_var
  nullify(lookup_var)

end subroutine CreateAddLookupTableVar

! ************************************************************************** !

subroutine LookupTableVarInitGradients(this,option)
  !
  ! allocate lookup variable gradients
  !
  ! Author: Paolo Orsini
  ! Date: 06/04/2018
  !
  
  use Option_module

  implicit none
  
  class(lookup_table_base_type) :: this
  type(option_type) :: option

  PetscInt :: prop_idx

  if ( this%dim > 0 ) then
    do prop_idx = 1,size(this%var_array(:))
      if ( associated(this%var_array(prop_idx)%ptr) ) then
        allocate(this%var_array(prop_idx)%ptr%sample_grad(this%dim))
        this%var_array(prop_idx)%ptr%sample_grad = UNINITIALIZED_DOUBLE
      end if  
    end do
  else
    option%io_buffer = "LookupTableVarInitGradients: cannot initialise " // &
                        "var gradient before defining table dims"
    call PrintErrMsg(option)
  end if  

end subroutine LookupTableVarInitGradients

! ************************************************************************** !

subroutine VarDataRead(this,input,num_fields,min_entries,error_string,option)
  !
  ! Reads in table data from input
  !
  ! Author: Paolo Orsini
  ! Date: 08/15/2018
  !

  use Option_module
  use Units_module
  use Input_Aux_module
  use Utility_module
  
  implicit none

  class(lookup_table_base_type) :: this
  type(input_type), pointer, intent(inout) :: input
  PetscInt, intent(in) :: num_fields
  PetscInt, intent(in) :: min_entries
  character(len=MAXSTRINGLENGTH), intent(in) :: error_string
  type(option_type), intent(inout) :: option

  PetscInt :: tmp_array_size 
  PetscReal, pointer :: tmp_data_array(:,:)
  PetscInt :: i_data, num_entries, i_entry
  character(len=MAXWORDLENGTH) :: word1

  tmp_array_size = 1000 !estimate max table points
  allocate(tmp_data_array(num_fields,tmp_array_size))
  tmp_data_array = UNINITIALIZED_DOUBLE
  
  num_entries = 0
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    if (InputError(input)) then
      option%io_buffer = 'Lookup table var_data reading'
      call PrintErrMsg(option)
    end if
    num_entries = num_entries + 1
    !press_idx = press_idx + 1
    if ( num_entries > tmp_array_size ) then
      !each time doubles the size of tmp_array_size
      !tmp_array_size overwritten by new size
      call ReallocateArray(tmp_data_array,tmp_array_size)
    end if
    do i_data = 1, num_fields
      call InputReadDouble(input,option,tmp_data_array(i_data,num_entries))
      call InputErrorMsg(input,option,'VALUE',error_string)
    end do
  end do

  if ( num_entries < min_entries ) then
    write(word1,*) min_entries
    option%io_buffer = trim(error_string) // &
                       ', number of entries less than = ' // trim(word1)
    call PrintErrMsg(option)
  end if

  allocate(this%var_data(num_fields,num_entries))
  do i_entry = 1,num_entries
    do i_data = 1, num_fields
      this%var_data(i_data,i_entry) = tmp_data_array(i_data,i_entry)
    end do
  end do
  
  call DeallocateArray(tmp_data_array)
  
end subroutine VarDataRead

! ************************************************************************** !

subroutine VarDataReverse(this,option)

  use Option_module

  implicit none

  class(lookup_table_base_type) :: this
  type(option_type), intent(inout) :: option

  PetscInt :: num_entries

  if ( .not. associated(this%var_data) ) then
    option%io_buffer = 'Cannot reverse VarData as not allocated'
    call PrintErrMsg(option)
  end if
  num_entries = size(this%var_data,2)

  this%var_data = this%var_data( :, num_entries:1:-1 )

end subroutine VarDataReverse

! ************************************************************************** !

subroutine LookupTableVarConvFactors(this,option)
  !
  ! Compute conversion varctors for all lookup table variables
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/2017
  !

  use Option_module
  use Units_module

  implicit none

  class(lookup_table_base_type) :: this
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: var

  var => this%vars%first

  do
    if (.not.associated(var)) exit
      var%conversion_factor = &
          UnitsConvertToInternal(var%user_units,var%internal_units,option)
      var => var%next
  enddo

end subroutine LookupTableVarConvFactors

! ************************************************************************** !

subroutine VarPointAndUnitConv(this,option)
  !
  ! Points variable data arrays to their table columns given the locations
  ! Convert variable data units given the convrsion factors 
  ! previously computed
  !
  ! Author: Paolo Orsini
  ! Date: 05/04/2018
  !

  use Option_module
  use Units_module

  implicit none

  class(lookup_table_base_type) :: this
  type(option_type) :: option

  PetscInt :: prop_idx, data_idx

  do prop_idx = 1,size(this%var_array(:))
    if ( associated(this%var_array(prop_idx)%ptr) ) then
      data_idx = this%var_array(prop_idx)%ptr%data_idx
      this%var_array(prop_idx)%ptr%data => this%var_data(data_idx,:)
      this%var_data(data_idx,:) = this%var_data(data_idx,:) * &
                              this%var_array(prop_idx)%ptr%conversion_factor
    end if
  end do

end subroutine VarPointAndUnitConv

! ************************************************************************** !

subroutine SetupConstGradExtrap(this,option)
  !
  ! Setup estrapolation method to constant gradient for all lookup vars
  !
  ! Author: Paolo Orsini
  ! Date: 05/04/2018
  !

  use Option_module

  implicit none

  class(lookup_table_base_type) :: this
  type(option_type) :: option

  PetscInt :: prop_idx
  
  do prop_idx = 1,size(this%var_array(:))
    if ( associated(this%var_array(prop_idx)%ptr) ) then
      this%var_array(prop_idx)%ptr%extrapolation_itype = VAR_EXTRAP_CONST_GRAD
    end if
  end do  
  
end subroutine SetupConstGradExtrap

! ************************************************************************** !

subroutine SetupConstValExtrap(this,option)
  !
  ! Setup extrapolation method to constant gradient for all lookup vars
  !
  ! Author: Paolo Orsini
  ! Date: 08/16/2018
  !
  use Option_module

  implicit none

  class(lookup_table_base_type) :: this
  type(option_type) :: option

  PetscInt :: prop_idx
  
  do prop_idx = 1,size(this%var_array(:))
    if ( associated(this%var_array(prop_idx)%ptr) ) then
      this%var_array(prop_idx)%ptr%extrapolation_itype = VAR_EXTRAP_CONST_VAL
    end if
  end do  
  
end subroutine SetupConstValExtrap
    
! ************************************************************************** !

subroutine SetupVarLinLogInterp(this,var_iname,option)
  ! 
  ! Author: Paolo Orsini
  ! Date: 06/06/18
  ! 
  
  use Option_module
  
  implicit none
  
  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname
  type(option_type) :: option
  
  if ( associated(this%var_array) ) then
    if ( associated(this%var_array(var_iname)%ptr) ) then
      this%var_array(var_iname)%ptr%interp_type = VAR_INTERP_X_LINLOG
    else
      option%io_buffer = "SetupVarLinLogInterp: cannot setup " // &
            "LinLog inteprolation method for a var not defined as lookupvar"
      call PrintErrMsg(option)
    end if  
  end if

end subroutine SetupVarLinLogInterp

! ************************************************************************** !

subroutine SetupVarUserUnits(this,var_iname,var_user_units,option)
  !
  ! Author: Paolo Orsini
  ! Date: 04/09/19
  !

  use Option_module

  implicit none

  class(lookup_table_base_type) :: this
  PetscInt, intent(in) :: var_iname
  character(len=MAXWORDLENGTH), intent(in) :: var_user_units
  type(option_type) :: option

  if ( associated(this%var_array) ) then
    if ( associated(this%var_array(var_iname)%ptr) ) then
      this%var_array(var_iname)%ptr%user_units = var_user_units
    else
      option%io_buffer = "SetupVarUserUnits: cannot setup " // &
            "User Units for a var not defined as lookupvar"
      call PrintErrMsg(option)
    end if
  end if

end subroutine SetupVarUserUnits

! ************************************************************************** !


subroutine LookupTableVarListDestroy(lookup_table_vars_list)
  !
  ! Deallocates a list of lookup table vars
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/2017
  !

  implicit none

  type(lookup_table_var_list_type), pointer :: lookup_table_vars_list

  type(lookup_table_var_type), pointer :: var, prev_var

  if (.not.associated(lookup_table_vars_list)) return

  var => lookup_table_vars_list%first
  do
    if (.not.associated(var)) exit
    prev_var => var
    var => var%next
    call LookupTableVarDestroy(prev_var)
  enddo

  lookup_table_vars_list%num_lookup_table_vars = 0
  nullify(lookup_table_vars_list%first)
  nullify(lookup_table_vars_list%last)

  deallocate(lookup_table_vars_list)
  nullify(lookup_table_vars_list)

end subroutine LookupTableVarListDestroy

! ************************************************************************** !

subroutine LookupTableVarDestroy(lookup_table_var)
  !
  ! Deallocates all members of a lookup table var
  !
  ! Author: Paolo Orsini
  ! Date: 12/01/2017 - mod 04/17/2018
  !
  use Utility_module

  implicit none

  type(lookup_table_var_type), pointer :: lookup_table_var

  if (associated(lookup_table_var%data)) then
    !data is only a pointer to a slice of var_data
    nullify(lookup_table_var%data)
  end if
  call DeallocateArray(lookup_table_var%sample_grad)
  nullify(lookup_table_var%next)

  deallocate(lookup_table_var)
  nullify(lookup_table_var)

end subroutine LookupTableVarDestroy

! ************************************************************************** !


end module Lookup_Table_module
