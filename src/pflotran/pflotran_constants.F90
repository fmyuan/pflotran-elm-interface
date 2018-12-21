module PFLOTRAN_Constants_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!

#include "petsc/finclude/petscsys.h"
  use petscsys

  use, intrinsic :: iso_fortran_env, only : stdout=>Output_Unit

  implicit none

  private

  PetscBool, parameter :: PFLOTRAN_RELEASE = PETSC_TRUE
  PetscInt, parameter :: PFLOTRAN_VERSION_MAJOR = 4
  PetscInt, parameter :: PFLOTRAN_VERSION_MINOR = 0
  PetscInt, parameter :: PFLOTRAN_VERSION_PATCH = 0 ! (alpha < -1; beta = -1)

#define VMAJOR 3
#define VMINOR 16
#define VSUBMINOR 2
#if (PETSC_VERSION_MAJOR < VMAJOR ||                    \
     (PETSC_VERSION_MAJOR == VMAJOR &&                  \
      (PETSC_VERSION_MINOR < VMINOR ||                  \
       (PETSC_VERSION_MINOR == VMINOR &&                \
        (PETSC_VERSION_SUBMINOR < VSUBMINOR)))))
!#error "Please use PETSc version 3.16.2 or later: 'git checkout v3.16.2' in $PETSC_DIR"
#endif
  ! MUST INCREMENT THIS NUMBER EVERYTIME A CHECKPOINT FILE IS
  ! MODIFIED TO PREVENT COMPATIBILITY ISSUES - geh.
  PetscInt, parameter, public :: CHECKPOINT_REVISION_NUMBER = 1

  PetscInt, parameter, public :: MAXSTRINGLENGTH = 512
  PetscInt, parameter, public :: MAXWORDLENGTH = 32
  PetscInt, parameter, public :: STDOUT_UNIT = stdout
  PetscInt, parameter, public :: DRIVER_OUT_UNIT = 14
  PetscInt, parameter, public :: FORWARD_OUT_UNIT = 15
  PetscInt, parameter, public :: OUTPUT_UNIT = 16
  PetscInt, parameter, public :: IN_UNIT = 17 ! 17-MAX_IN_UNIT are reserved
  ! If you increase MAX_IN_UNIT, you MUST ensure that no other units #
  ! lie between IN_UNIT and MAX_IN_UNIT, as these units are reserved
  ! for embedded input files.
  PetscInt, parameter, public :: MAX_IN_UNIT = 25
  PetscInt, parameter, public :: IUNIT_TEMP = 86
  ! Units 50-59 are reserved for reservoir engineering format files
  ! EKG_UNIT = 87
  PetscInt, parameter, public :: INPUT_RECORD_UNIT = 88
  PetscInt, parameter, public :: HHISTORY_LENGTH = 1000
  ! HHISTORY_LENGTH is the length of the array used to store the differencing
  ! values h.

  ! EXIT codes
  PetscInt, parameter, public :: EXIT_SUCCESS = 0
  PetscInt, parameter, public :: EXIT_USER_ERROR = 87
  PetscInt, parameter, public :: EXIT_FAILURE = 88

  ! formula weights
  PetscReal, parameter, public :: FMWNACL = 58.44277d0
  PetscReal, parameter, public :: FMWH2O = 18.01534d0  ! kg/kmol h2o
  PetscReal, parameter, public :: FMWCO2 = 44.0098d0
  PetscReal, parameter, public :: FMWAIR = 28.96d0

  ! constants
  PetscReal, parameter, public :: DAYS_PER_YEAR = 365.d0
!geh: for bragflo year
!  PetscReal, parameter, public :: DAYS_PER_YEAR = 365.24224537d0
  PetscReal, parameter, public :: H2O_CRITICAL_TEMPERATURE = 647.3d0  ! K
#if defined(MATCH_TOUGH2)
  PetscReal, parameter, public :: H2O_CRITICAL_PRESSURE = 22.12d6 ! Pa
#else
  PetscReal, parameter, public :: H2O_CRITICAL_PRESSURE = 22.064d6 ! Pa
#endif

  ! conversion factors
  PetscReal, parameter, public :: LOG_TO_LN = 2.30258509299d0
  PetscReal, parameter, public :: LN_TO_LOG = 0.434294481904d0

  ! constants
                             ! from http://physics.nist.gov/cgi-bin/cuu/Value?r
  PetscReal, parameter, public :: IDEAL_GAS_CONSTANT = 8.31446d0 ! J/mol-K
!to match BRAGFLO  PetscReal, parameter, public :: IDEAL_GAS_CONSTANT = 8.31451 ! J/mol-K
  PetscReal, parameter, public :: HEAT_OF_FUSION = 3.34d5  ! J/kg
  PetscReal, parameter, public :: PI = 3.14159265359d0
  PetscReal, parameter, public :: FARADAY = 96485.3365d0 ! C/mol
  PetscReal, parameter, public :: EARTH_GRAVITY = 9.8068d0 ! m/s^2

  PetscInt, parameter, public :: ZERO_INTEGER = 0
  PetscInt, parameter, public :: ONE_INTEGER = 1
  PetscInt, parameter, public :: TWO_INTEGER = 2
  PetscInt, parameter, public :: THREE_INTEGER = 3
  PetscInt, parameter, public :: FOUR_INTEGER = 4
  PetscInt, parameter, public :: FIVE_INTEGER = 5
  PetscInt, parameter, public :: SIX_INTEGER = 6
  PetscInt, parameter, public :: SEVEN_INTEGER = 7
  PetscInt, parameter, public :: EIGHT_INTEGER = 8
  PetscInt, parameter, public :: NINE_INTEGER = 9
  PetscInt, parameter, public :: TEN_INTEGER = 10
  PetscInt, parameter, public :: ELEVEN_INTEGER = 11
  PetscInt, parameter, public :: TWELVE_INTEGER = 12
  PetscInt, parameter, public :: NEG_ONE_INTEGER = -1

  PetscMPIInt, parameter, public :: ZERO_INTEGER_MPI = ZERO_INTEGER
  PetscMPIInt, parameter, public :: ONE_INTEGER_MPI = ONE_INTEGER
  PetscMPIInt, parameter, public :: TWO_INTEGER_MPI = TWO_INTEGER
  PetscMPIInt, parameter, public :: THREE_INTEGER_MPI = THREE_INTEGER
  PetscMPIInt, parameter, public :: FOUR_INTEGER_MPI = FOUR_INTEGER
  PetscMPIInt, parameter, public :: FIVE_INTEGER_MPI = FIVE_INTEGER
  PetscMPIInt, parameter, public :: SIX_INTEGER_MPI = SIX_INTEGER
  PetscMPIInt, parameter, public :: SEVEN_INTEGER_MPI = SEVEN_INTEGER
  PetscMPIInt, parameter, public :: TEN_INTEGER_MPI = TEN_INTEGER
  PetscMPIInt, parameter, public :: ELEVEN_INTEGER_MPI = ELEVEN_INTEGER
  PetscMPIInt, parameter, public :: TWELVE_INTEGER_MPI = TWELVE_INTEGER
  PetscMPIInt, parameter, public :: MAXSTRINGLENGTH_MPI = MAXSTRINGLENGTH

  PetscInt, parameter, public :: X_DIRECTION = 1
  PetscInt, parameter, public :: Y_DIRECTION = 2
  PetscInt, parameter, public :: Z_DIRECTION = 3
  PetscInt, parameter, public :: XY_DIRECTION = 4
  PetscInt, parameter, public :: XZ_DIRECTION = 5
  PetscInt, parameter, public :: YZ_DIRECTION = 6
  PetscInt, parameter, public :: LOWER = 1
  PetscInt, parameter, public :: UPPER = 2

  PetscInt, parameter, public :: TIME_NULL = 0
  PetscInt, parameter, public :: TIME_T = 1
  PetscInt, parameter, public :: TIME_TpDT = 2

  PetscInt, parameter, public :: SORPTION_LINEAR = 1
  PetscInt, parameter, public :: SORPTION_LANGMUIR = 2
  PetscInt, parameter, public :: SORPTION_FREUNDLICH  = 3

  ! Classes
  PetscInt, parameter, public :: NULL_CLASS = 0
  PetscInt, parameter, public :: FLOW_CLASS = 1
  PetscInt, parameter, public :: TRANSPORT_CLASS = 2
  PetscInt, parameter, public :: GEOPHYSICS_CLASS = 3

  ! Macros that are used as 'dm_index' values.  --RTM
  PetscInt, parameter, public :: ONEDOF = 1
  PetscInt, parameter, public :: NPHASEDOF = 2
  PetscInt, parameter, public :: THREENPDOF = 3
  PetscInt, parameter, public :: NFLOWDOF = 4
  PetscInt, parameter, public :: NTRANDOF = 5
  PetscInt, parameter, public :: NGEODOF = 7

  PetscInt, parameter, public :: GLOBAL = 1
  PetscInt, parameter, public :: LOCAL = 2
  PetscInt, parameter, public :: NATURAL = 3

  PetscInt, parameter, public :: NULL_MODE = 0

  ! flow modes
  PetscInt, parameter, public :: MPH_MODE = 1
  PetscInt, parameter, public :: RICHARDS_MODE = 2
  PetscInt, parameter, public :: G_MODE = 3
  PetscInt, parameter, public :: TH_MODE = 4
  PetscInt, parameter, public :: WF_MODE = 5
  PetscInt, parameter, public :: RICHARDS_TS_MODE = 6
  PetscInt, parameter, public :: TH_TS_MODE = 7
  PetscInt, parameter, public :: H_MODE = 8
  PetscInt, parameter, public :: ZFLOW_MODE = 9
  PetscInt, parameter, public :: PNF_MODE = 10

  ! transport modes
  PetscInt, parameter, public :: RT_MODE = 1
  PetscInt, parameter, public :: NWT_MODE = 2
  PetscInt, parameter, public :: EXPLICIT_ADVECTION = 10

  ! geophysics modes
  PetscInt, parameter, public :: ERT_MODE = 1
  PetscInt, parameter, public :: SIP_MODE = 2

  ! condition types
  PetscInt, parameter, public :: NULL_CONDITION = 0
  PetscInt, parameter, public :: DIRICHLET_BC = 1
  PetscInt, parameter, public :: NEUMANN_BC = 2
  PetscInt, parameter, public :: DIRICHLET_ZERO_GRADIENT_BC = 3
  PetscInt, parameter, public :: ZERO_GRADIENT_BC = 4
  PetscInt, parameter, public :: HYDROSTATIC_BC = 5
  PetscInt, parameter, public :: HYDROSTATIC_SEEPAGE_BC = 6
  PetscInt, parameter, public :: MASS_RATE_SS = 7
  PetscInt, parameter, public :: VOLUMETRIC_RATE_SS = 8
  PetscInt, parameter, public :: SCALED_MASS_RATE_SS = 9
  PetscInt, parameter, public :: SCALED_VOLUMETRIC_RATE_SS = 10
  PetscInt, parameter, public :: CONCENTRATION_SS = 11
  PetscInt, parameter, public :: EQUILIBRIUM_SS = 12
  PetscInt, parameter, public :: HYDROSTATIC_CONDUCTANCE_BC = 13
  PetscInt, parameter, public :: UNIT_GRADIENT_BC = 14
  PetscInt, parameter, public :: SATURATION_BC = 15
  PetscInt, parameter, public :: HET_VOL_RATE_SS = 16
  PetscInt, parameter, public :: HET_MASS_RATE_SS = 17
  PetscInt, parameter, public :: HET_DIRICHLET_BC = 18
  PetscInt, parameter, public :: ENERGY_RATE_SS = 19
  PetscInt, parameter, public :: SCALED_ENERGY_RATE_SS = 20
  PetscInt, parameter, public :: HET_ENERGY_RATE_SS = 21
  PetscInt, parameter, public :: HET_SURF_HYDROSTATIC_SEEPAGE_BC = 22
  PetscInt, parameter, public :: SPILLOVER_BC = 23
  PetscInt, parameter, public :: SURFACE_DIRICHLET = 24
  PetscInt, parameter, public :: SURFACE_ZERO_GRADHEIGHT = 25
  PetscInt, parameter, public :: SURFACE_SPILLOVER = 26
  PetscInt, parameter, public :: HET_HYDROSTATIC_SEEPAGE_BC = 27
  PetscInt, parameter, public :: HET_HYDROSTATIC_CONDUCTANCE_BC = 28
  PetscInt, parameter, public :: TOTAL_MASS_RATE_SS = 29
  PetscInt, parameter, public :: DIRICHLET_SEEPAGE_BC = 30
  PetscInt, parameter, public :: DIRICHLET_CONDUCTANCE_BC = 31

  PetscInt, parameter, public :: WELL_SS = 100

  ! source/sink scaling options
  PetscInt, parameter, public :: SCALE_BY_PERM = 1
  PetscInt, parameter, public :: SCALE_BY_NEIGHBOR_PERM = 2
  PetscInt, parameter, public :: SCALE_BY_VOLUME = 3

  ! connection types
  PetscInt, parameter, public :: INTERNAL_CONNECTION_TYPE = 1
  PetscInt, parameter, public :: BOUNDARY_CONNECTION_TYPE = 2
  PetscInt, parameter, public :: INITIAL_CONNECTION_TYPE = 3
  PetscInt, parameter, public :: SRC_SINK_CONNECTION_TYPE = 4

  ! dofs for each mode
  PetscInt, parameter, public :: THC_PRESSURE_DOF = 1
  PetscInt, parameter, public :: THC_TEMPERATURE_DOF = 2
  PetscInt, parameter, public :: THC_CONCENTRATION_DOF = 3
  PetscInt, parameter, public :: THC_MASS_RATE_DOF = 4
  PetscInt, parameter, public :: THC_ENTHALPY_DOF = 5

  PetscInt, parameter, public :: TH_PRESSURE_DOF = 1
  PetscInt, parameter, public :: TH_TEMPERATURE_DOF = 2
  PetscInt, parameter, public :: TH_CONDUCTANCE_DOF = 3

  PetscInt, parameter, public :: MPH_PRESSURE_DOF = 1
  PetscInt, parameter, public :: MPH_TEMPERATURE_DOF = 2
  PetscInt, parameter, public :: MPH_CONCENTRATION_DOF = 3

  PetscInt, parameter, public :: RICHARDS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: RICHARDS_CONDUCTANCE_DOF = 2

  PetscInt, parameter, public :: ZFLOW_PRESSURE_DOF = 1
  PetscInt, parameter, public :: ZFLOW_CONDUCTANCE_DOF = 2

  PetscInt, parameter, public :: MIS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: MIS_CONCENTRATION_DOF = 2

  ! mphase equation of state
  PetscInt, parameter, public :: EOS_SPAN_WAGNER = 1
  PetscInt, parameter, public :: EOS_MRK = 2

  ! phase ids
  PetscInt, parameter, public :: LIQUID_PHASE = 1
  PetscInt, parameter, public :: GAS_PHASE = 2

  PetscInt, parameter, public :: MAX_PHASE = 2

  ! approaches to coupling reactive transport
  PetscInt, parameter, public :: GLOBAL_IMPLICIT = 0
  PetscInt, parameter, public :: OPERATOR_SPLIT = 1

  ! ids of non-petsc arrays
  PetscInt, parameter, public :: MATERIAL_ID_ARRAY = 1
  PetscInt, parameter, public :: CC_ID_ARRAY = 2  ! characteristic curves
  PetscInt, parameter, public :: CCT_ID_ARRAY = 3 ! charact. curves thermal

  ! interpolation methods
  PetscInt, parameter, public :: INTERPOLATION_NULL = 0
  PetscInt, parameter, public :: INTERPOLATION_STEP = 1
  PetscInt, parameter, public :: INTERPOLATION_LINEAR = 2

  ! print secondary continuum variable ids
  PetscInt, parameter, public :: PRINT_SEC_TEMP =           0
  PetscInt, parameter, public :: PRINT_SEC_CONC =           1
  PetscInt, parameter, public :: PRINT_SEC_MIN_VOLFRAC =    2
  PetscInt, parameter, public :: PRINT_SEC_MIN_RATE =       3
  PetscInt, parameter, public :: PRINT_SEC_MIN_SI =         4

  PetscInt, parameter, public :: PROCEED = 0
  PetscInt, parameter, public :: DONE = 1
  PetscInt, parameter, public :: FAIL = 2

  ! Grid type
  PetscInt, parameter, public :: NULL_GRID = 0
  PetscInt, parameter, public :: STRUCTURED_GRID = 1
  PetscInt, parameter, public :: UNSTRUCTURED_GRID = 2
  PetscInt, parameter, public :: IMPLICIT_UNSTRUCTURED_GRID = 3
  PetscInt, parameter, public :: EXPLICIT_UNSTRUCTURED_GRID = 4
  PetscInt, parameter, public :: POLYHEDRA_UNSTRUCTURED_GRID = 5
  PetscInt, parameter, public :: ONE_DIM_GRID = 1
  PetscInt, parameter, public :: TWO_DIM_GRID = 2
  PetscInt, parameter, public :: THREE_DIM_GRID = 3
  PetscInt, parameter, public :: VERTEX_CENTERED_OUTPUT_MESH = 1
  PetscInt, parameter, public :: CELL_CENTERED_OUTPUT_MESH = 2

  ! Geomechanics
  PetscInt, parameter, public :: GEOMECH_DISP_X_DOF = 1
  PetscInt, parameter, public :: GEOMECH_DISP_Y_DOF = 2
  PetscInt, parameter, public :: GEOMECH_DISP_Z_DOF = 3
  PetscInt, parameter, public :: GEOMECH_ONE_WAY_COUPLED = 4
  PetscInt, parameter, public :: GEOMECH_TWO_WAY_COUPLED = 5

  ! Macros that are used as 'vscatter_index' values
  PetscInt, parameter, public :: SUBSURF_TO_GEOMECHANICS = 3
  PetscInt, parameter, public :: GEOMECHANICS_TO_SUBSURF = 4

  ! Ice/water/vapor partitioning model
  PetscInt, parameter, public :: PAINTER_EXPLICIT = 1
  PetscInt, parameter, public :: PAINTER_KARRA_IMPLICIT = 2
  PetscInt, parameter, public :: PAINTER_KARRA_EXPLICIT = 3
  PetscInt, parameter, public :: DALL_AMICO = 4
  PetscInt, parameter, public :: PAINTER_KARRA_EXPLICIT_NOCRYO = 5

  ! Relative permeability averaging
  PetscInt, parameter, public :: UPWIND = 1
  PetscInt, parameter, public :: HARMONIC = 2
  PetscInt, parameter, public :: DYNAMIC_HARMONIC = 3

  ! uninitialized values
  PetscInt, parameter, public :: UNINITIALIZED_INTEGER = -999
  PetscReal, parameter, public :: UNINITIALIZED_DOUBLE = -999.d0
  
  ! maximum values
  PetscReal, parameter, public :: MAX_DOUBLE = 1.d20

  ! global solver convergence criteria
  PetscInt, parameter, public :: CONVERGENCE_OFF = -999
  PetscInt, parameter, public :: CONVERGENCE_CUT_TIMESTEP = -1
  PetscInt, parameter, public :: CONVERGENCE_KEEP_ITERATING = 0
  PetscInt, parameter, public :: CONVERGENCE_FORCE_ITERATION = 1
  PetscInt, parameter, public :: CONVERGENCE_CONVERGED = 2

  ! Dummy value
  PetscReal, parameter, public :: DUMMY_VALUE = UNINITIALIZED_DOUBLE

  interface Uninitialized
    module procedure UninitializedInteger
    module procedure UninitializedDouble
    module procedure UninitializedMatType
  end interface

  interface Initialized
    module procedure InitializedInteger
    module procedure InitializedDouble
    module procedure InitializedMatType
  end interface

  public :: Initialized, &
            Uninitialized, &
            UninitializedMessage, &
            GetVersion

contains

! ************************************************************************** !

function InitializedInteger(value)
  !
  ! Tests whether a variable is initialized based orginally being set to
  ! the value UNINITIALIZED_INTEGER
  !
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  implicit none

  PetscInt :: value
  PetscBool :: InitializedInteger

  InitializedInteger = .not.Uninitialized(value)

end function InitializedInteger


! ************************************************************************** !

function UninitializedInteger(value)
  !
  ! Tests whether a variable is uninitialized based orginally being set to
  ! the value UNINITIALIZED_INTEGER
  !
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  implicit none

  PetscInt :: value
  PetscBool :: UninitializedInteger

  UninitializedInteger = (value == UNINITIALIZED_INTEGER)

end function UninitializedInteger

! ************************************************************************** !

function InitializedDouble(value)
  !
  ! Tests whether a variable is initialized based orginally being set to
  ! the value UNINITIALIZED_INTEGER
  !
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  implicit none

  PetscReal :: value
  PetscBool :: InitializedDouble

  InitializedDouble = .not.Uninitialized(value)

end function InitializedDouble

! ************************************************************************** !

function UninitializedDouble(value)
  !
  ! Tests whether a variable is uninitialized based orginally being set to
  ! the value UNINITIALIZED_INTEGER
  !
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  implicit none

  PetscReal :: value
  PetscBool :: UninitializedDouble

  UninitializedDouble = (dabs(value-UNINITIALIZED_DOUBLE) < 1.d-20)

end function UninitializedDouble

! ************************************************************************** !

function InitializedMatType(value)
  !
  ! Tests whether a variable is initialized based orginally being set to
  ! the value PETSC_NULL_CHARACTER.
  !
#include "petsc/finclude/petscmat.h"
  use petscmat

  implicit none

  MatType :: value
  PetscBool :: InitializedMatType

  InitializedMatType = .not.Uninitialized(value)

end function InitializedMatType

! ************************************************************************** !

function UninitializedMatType(value)
  !
  ! Tests whether a variable is uninitialized based orginally being set to
  ! the value PETSC_NULL_CHARACTER.
  !
#include "petsc/finclude/petscmat.h"
  use petscmat

  implicit none

  MatType :: value
  PetscBool :: UninitializedMatType

  UninitializedMatType = (value == PETSC_NULL_CHARACTER)

end function UninitializedMatType


! ************************************************************************** !

function UninitializedMessage(variable_name,routine_name)
  !
  ! Tests whether a variable is uninitialized based orginally being set to
  ! the value UNINITIALIZED_INTEGER
  !
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  implicit none

  character(len=*) :: variable_name
  character(len=*) :: routine_name

  character(len=MAXSTRINGLENGTH) :: UninitializedMessage

  if (len_trim(routine_name) > 1) then
    UninitializedMessage = trim(variable_name) // &
                           ' uninitialized in ' // &
                           trim(routine_name) // '.'
  else
    UninitializedMessage = trim(variable_name) // &
                           ' uninitialized.'
  endif

end function UninitializedMessage

! ************************************************************************** !

function GetVersion()
  !
  ! Returns the PFLOTRAN version in string format using semantic versioning
  !
  ! Author: Glenn Hammond
  ! Date: 04/23/20
  !
  implicit none

  character(len=MAXWORDLENGTH) :: GetVersion

  character(len=MAXWORDLENGTH) :: word

  if (PFLOTRAN_RELEASE) then
    write(word,*) PFLOTRAN_VERSION_MAJOR
    GetVersion = 'PFLOTRAN v' // trim(adjustl(word))
    write(word,*) PFLOTRAN_VERSION_MINOR
    GetVersion = trim(GetVersion) // '.' // trim(adjustl(word))
    if (PFLOTRAN_VERSION_PATCH > 0) then
      write(word,*) PFLOTRAN_VERSION_PATCH
      GetVersion = trim(GetVersion) // '.' // trim(adjustl(word))
    else if (PFLOTRAN_VERSION_PATCH < -1) then
      GetVersion = trim(GetVersion) // '-alpha'
    else if (PFLOTRAN_VERSION_PATCH < 0) then
      GetVersion = trim(GetVersion) // '-beta'
    endif
  else
    GetVersion = 'PFLOTRAN Development Version'
  endif

end function GetVersion

end module PFLOTRAN_Constants_module
