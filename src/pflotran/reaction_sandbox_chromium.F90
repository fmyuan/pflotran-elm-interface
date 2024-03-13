module Reaction_Sandbox_Chromium_class

! Sandbox reaction for Cr(VI) reduction using bio-reduction with reduced permeability

  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_chromium_type
    character(len=MAXWORDLENGTH) :: name_D_mobile   ! This the mobile food in the reaction
    character(len=MAXWORDLENGTH) :: name_D_immobile ! This the immobile food in the reaction
    character(len=MAXWORDLENGTH) :: name_C     ! This is for Cr(VI) in the reaction
    character(len=MAXWORDLENGTH) :: name_B ! This is for the biomass in the reaction
    character(len=MAXWORDLENGTH) :: name_I ! This is for the alcohol in the reaction
    character(len=MAXWORDLENGTH) :: name_X ! This is for the biocide in the reaction
    character(len=MAXWORDLENGTH) :: name_biomineral ! This is for the dummny bio mineral

    PetscInt :: B_id
    PetscInt :: C_id
    PetscInt :: D_immobile_id  ! Immobile food
    PetscInt :: D_mobile_id    ! Mobile food
    PetscInt :: I_id
    PetscInt :: X_id
    PetscInt :: biomineral_id

    ! Decay and inhibition parameters in our sophisticated model
    PetscReal :: background_concentration_B    ! Minimum background concentration of the biomass
    PetscReal :: rate_B_1
    PetscReal :: rate_B_2
    PetscReal :: rate_C
    PetscReal :: rate_D
    PetscReal :: rate_D_i
    PetscReal :: rate_D_m
    PetscReal :: inhibition_B
    PetscReal :: inhibition_C
    PetscReal :: monod_D
    PetscReal :: inhibition_I
    PetscReal :: mass_action_B
    PetscReal :: mass_action_CD
    PetscReal :: mass_action_X
    PetscReal :: stoichiometric_C
    PetscReal :: stoichiometric_D_1
    PetscReal :: stoichiometric_D_2
    PetscReal :: exponent_B
    PetscReal :: density_B

  contains
    procedure, public :: ReadInput => ChromiumRead
    procedure, public :: Setup => ChromiumSetup
    procedure, public :: Evaluate => ChromiumReact
    procedure, public :: UpdateKineticState => ChromiumKineticState
    procedure, public :: Destroy => ChromiumDestroy
  end type reaction_sandbox_chromium_type

  public :: ChromiumCreate

contains

! ************************************************************************** !

function ChromiumCreate()
  !
  ! Allocates Cr(VI) bio-reduction object
  !
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  !

  implicit none

  class(reaction_sandbox_chromium_type), pointer :: ChromiumCreate

  allocate(ChromiumCreate)
  ChromiumCreate%name_D_mobile = ''
  ChromiumCreate%name_D_immobile = ''
  ChromiumCreate%name_C = ''
  ChromiumCreate%name_B = ''
  ChromiumCreate%name_I = ''
  ChromiumCreate%name_X = ''
  ChromiumCreate%name_biomineral = ''

  ChromiumCreate%D_mobile_id = 0
  ChromiumCreate%D_immobile_id = 0
  ChromiumCreate%C_id = 0
  ChromiumCreate%B_id = 0
  ChromiumCreate%I_id = 0
  ChromiumCreate%X_id = 0
  ChromiumCreate%biomineral_id = 0

  ChromiumCreate%stoichiometric_D_1 = 0.d0
  ChromiumCreate%rate_D = 0.d0
  ChromiumCreate%rate_C = 0.d0
  ChromiumCreate%inhibition_C = 0.d0
  ChromiumCreate%rate_B_2 = 0.d0
  ChromiumCreate%rate_B_1 = 0.d0
  ChromiumCreate%monod_D = 0.d0
  ChromiumCreate%inhibition_B = 0.d0
  ChromiumCreate%background_concentration_B = 0.d0
  ChromiumCreate%mass_action_CD = 0.d0
  ChromiumCreate%stoichiometric_C = 1.d0
  ChromiumCreate%stoichiometric_D_2 = 1.d0
  ChromiumCreate%rate_D_i = 0.d0
  ChromiumCreate%rate_D_m = 0.d0
  ChromiumCreate%exponent_B = 0.d0
  ChromiumCreate%inhibition_I = 0.d0
  ChromiumCreate%mass_action_B = 0.d0
  ChromiumCreate%mass_action_X = 0.d0
  ChromiumCreate%density_B = 0.d0
  nullify(ChromiumCreate%next)

end function ChromiumCreate

! ************************************************************************** !

subroutine ChromiumRead(this,input,option)
  !
  ! Reads input deck for Cr(VI) bio-reduction parameters
  !
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  class(reaction_sandbox_chromium_type) :: this
  type(input_type), pointer  :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
    call StringToUpper(word)

    select case(trim(word))
      case('NAME_B')
        call InputReadWord(input,option,this%name_B,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_B', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('NAME_C')
        call InputReadWord(input,option,this%name_C,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_C', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('NAME_D_MOBILE')
        call InputReadWord(input,option,this%name_D_mobile,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_D_mobile', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('NAME_D_IMMOBILE')
        call InputReadWord(input,option,this%name_D_immobile,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_D_immobile', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('NAME_I')
        call InputReadWord(input,option,this%name_I,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_I', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('NAME_X')
        call InputReadWord(input,option,this%name_X,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_X', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('NAME_BIOMINERAL')
        call InputReadWord(input,option,this%name_biomineral,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_biomineral', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')

      case('EXPONENT_B')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%exponent_B)
        call InputErrorMsg(input,option,'exponent_B', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')

      case('BACKGROUND_CONC_B')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%background_concentration_B )
        call InputErrorMsg(input,option,'background_concentration_B', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')

      case('MASS_ACTION_B')
        call InputReadDouble(input,option,this%mass_action_B)
        call InputErrorMsg(input,option,'mass_action_B', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('MASS_ACTION_CD')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%mass_action_CD)
        call InputErrorMsg(input,option,'mass_action_CD', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('MASS_ACTION_X')
        call InputReadDouble(input,option,this%mass_action_X)
        call InputErrorMsg(input,option,'mass_action_X', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')

      case('INHIBITION_B')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%inhibition_B)
        call InputErrorMsg(input,option,'inhibition_B', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('INHIBITION_C')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%inhibition_C)
        call InputErrorMsg(input,option,'inhibition_C', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('MONOD_D')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%monod_D)
        call InputErrorMsg(input,option,'monod_D', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('INHIBITION_I')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%inhibition_I)
        call InputErrorMsg(input,option,'inhibition_I', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')

      case('RATE_B_1')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_B_1)
        call InputErrorMsg(input,option,'rate_B_1', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('RATE_B_2')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_B_2)
        call InputErrorMsg(input,option,'rate_B_2', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('RATE_C')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_C)
        call InputErrorMsg(input,option,'rate_C', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('RATE_D')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_D)
        call InputErrorMsg(input,option,'rate_D', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('RATE_D_IMMOB')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%rate_D_i)
        call InputErrorMsg(input,option,'rate_D_i', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('RATE_D_MOBIL')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%rate_D_m)
        call InputErrorMsg(input,option,'rate_D_m', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')

      case('DENSITY_B')
        call InputReadDouble(input,option,this%density_B)
        call InputErrorMsg(input,option,'density_B', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')

      case('STOICHIOMETRIC_C')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%stoichiometric_C )
        call InputErrorMsg(input,option,'stoichiometric_C', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('STOICHIOMETRIC_D_1')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%stoichiometric_D_1)
        call InputErrorMsg(input,option,'stoichiometric_D_1', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')
      case('STOICHIOMETRIC_D_2')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%stoichiometric_D_2)
        call InputErrorMsg(input,option,'stoichiometric_D_2', &
                           'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION')

      case default
        call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,CHROMIUM_REDUCTION',option)
    end select
  enddo

end subroutine ChromiumRead

! ************************************************************************** !

subroutine ChromiumSetup(this,reaction,option)
  !
  ! Sets up the Cr(VI) bio-reduction reaction
  !
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  !

  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module
  use Reaction_Mineral_Aux_module
  use Option_module

  implicit none

  class(reaction_sandbox_chromium_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  this%D_mobile_id = &
    ReactionAuxGetPriSpecIDFromName(this%name_D_mobile, &
                                    reaction,option)

  this%C_id = &
    ReactionAuxGetPriSpecIDFromName(this%name_C, &
                                    reaction,option)

  this%I_id = &
    ReactionAuxGetPriSpecIDFromName(this%name_I, &
                                    reaction,option)

  this%X_id = &
    ReactionAuxGetPriSpecIDFromName(this%name_X, &
                                    reaction,option)

  this%B_id = &
    ReactionImGetSpeciesIDFromName(this%name_B, &
                                   reaction%immobile,option)

  this%D_immobile_id = &
    ReactionImGetSpeciesIDFromName(this%name_D_immobile, &
                                   reaction%immobile,option)
  this%biomineral_id = &
    ReactionMnrlGetMnrlIDFromName(this%name_biomineral, &
                                  reaction%mineral,option)

end subroutine ChromiumSetup

! ************************************************************************** !

subroutine ChromiumReact(this,Residual,Jacobian,compute_derivative, &
                       rt_auxvar,global_auxvar,material_auxvar, &
                       reaction,option)
  !
  ! Evaluates reaction storing residual and/or Jacobian
  !
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_chromium_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water
  PetscReal :: mu_B, mu_CD
  PetscReal :: sum_food
  PetscInt :: idof_food_mobile, idof_food_immobile, idof_biomass, idof_Cr
  PetscInt :: idof_alcohol, idof_biocide
  PetscReal :: immobile_to_water_vol
  PetscReal :: immobile_mole_fraction, mobile_mole_fraction
  PetscReal :: biomass_residual_delta

  ! Description of subroutine arguments:

  ! Residual - 1D array storing residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entires in units kg water/sec
  !
  !  Jacobian [kg water/sec] * dc [mol/kg water] = -Res [mol/sec]
  !
  ! compute_derivative - Flag indicating whether analtical derivative should
  !   be calculated.  The user must provide either the analytical derivatives
  !   or a numerical approximation unless always running with
  !   NUMERICAL_JACOBIAN_RXN defined in input deck.  If the use of
  !   NUMERICAL_JACOBIAN_RXN is assumed, the user should provide an error
  !   message when compute_derivative is true.  E.g.
  !
  !   option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
  !                      'due to assumptions in ChromeAlcohol'
  !   call printErrMsg(option)
  !
  ! rt_auxvar - Object holding chemistry information (e.g. concentrations,
  !   activity coefficients, mineral volume fractions, etc.).  See
  !   reactive_transport_aux.F90.
  !
  !   Useful variables:
  !     rt_auxvar%total(:,iphase) - total component concentrations
  !                                 [mol/L water] for phase
  !     rt_auxvar%pri_molal(:) - free ion concentrations [mol/kg water]
  !     rt_auxvar%pri_act_coef(:) - activity coefficients for primary species
  !     rt_auxvar%aqueous%dtotal(:,iphase) - derivative of total component
  !                 concentration with respect to free ion [kg water/L water]
  !
  ! global_auxvar - Object holding information on flow (e.g. saturation,
  !   density, viscosity, temperature, etc)
  !
  !   Useful variables:
  !     global_auxvar%den(iphase) - fluid density [mol/m^3]
  !     global_auxvar%den_kg(iphase) - fluid density [kg/m^3]
  !     global_auxvar%sat(iphase) - saturation
  !     global_auxvar%temp - temperature [C]
  !
  ! porosity - effective porosity of grid cell [m^3 pore/m^3 bulk]
  ! volume - volume of grid cell [m^3]
  ! reaction - Provides access to variable describing chemistry.  E.g.
  !   reaction%ncomp - # chemical degrees of freedom (mobile and immobile)
  !   reaction%naqcomp - # chemical degrees of freedom on water
  !   reaction%primary_species_names(:) - names of primary species
  !
  ! option - Provides handle for controlling simulation, catching and
  !          reporting errors.

  ! Unit of the residual must be in moles/second
  ! global_auxvar%sat(iphase) = saturation of cell
  ! 1.d3 converts m^3 water -> L water
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  ! always subtract contribution from residual

  idof_food_mobile = this%D_mobile_id
  idof_Cr = this%C_id
  idof_alcohol = this%I_id
  idof_biocide = this%X_id
  idof_biomass = reaction%offset_immobile + this%B_id
  idof_food_immobile = reaction%offset_immobile + this%D_immobile_id

  immobile_to_water_vol = &
     material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0                   ! L water/ m3 bulk

  sum_food = rt_auxvar%total(idof_food_mobile,iphase) + &
             rt_auxvar%immobile(this%D_immobile_id)/ &
             immobile_to_water_vol                                                ! in mol/L water; Note that food_immobile is divided by porosity*saturation

  mu_B = this%rate_B_1*rt_auxvar%immobile(this%B_id)* &      ! mol/m3 bulk/s
        ! F monod term, unitless
        (sum_food/(sum_food + this%monod_D))* &
        ! B monod inhibition term, unitless
        (this%inhibition_B/ &
        (rt_auxvar%immobile(this%B_id) + &
         this%inhibition_B))**this%exponent_B* &
        ! I Monod inhibition term, unitless
        (this%inhibition_I/ &
        (this%inhibition_I + &
         rt_auxvar%total(idof_alcohol,iphase)))

  mu_CD = this%mass_action_CD*sum_food*rt_auxvar%total(idof_Cr,iphase)    ! mol/L/s

  Residual(idof_Cr) =      Residual(idof_Cr) + &
                           ! Biological reaction, mol/s
                           this%rate_C* &                              ! /s
                           rt_auxvar%immobile(this%B_id)* &                 ! mol/m3 bulk
                           rt_auxvar%total(idof_Cr,iphase)/ &                     ! mol/L
                           (this%inhibition_C + &
                           rt_auxvar%total(idof_Cr,iphase))* &                    ! mol/L
                           material_auxvar%volume + &                             ! m3 bulk
                           ! Direct reaction, mol/s
                           this%stoichiometric_C*mu_CD* &                 ! mol/L water/s
                           material_auxvar%volume* &                              ! m3 bulk
                           immobile_to_water_vol                                  ! L water/m3 bulk

  biomass_residual_delta = &                                                      ! Growth usage, mol/s
                           - mu_B*material_auxvar%volume + &                      ! mol/m3 bulk/s * m3 bulk
                           ! Natural decay, mol/s
                           this%rate_B_2* &                         ! 1/s
                           (rt_auxvar%immobile(this%B_id) - &
                            this%background_concentration_B)* &                                  ! mol/m3 bulk
                           material_auxvar%volume + &                             ! m3 bulk
                           ! Biocide reaction, mol/s
                           this%mass_action_B* &                 ! L/mol/s
                           (rt_auxvar%immobile(this%B_id) - &
                           this%background_concentration_B)* &                                   ! mol/m3 bulk
                           rt_auxvar%total(idof_biocide,iphase)* &                ! mol/L
                           material_auxvar%volume                                 ! m3 bulk

  Residual(idof_biomass) = Residual(idof_biomass) + biomass_residual_delta        ! mol/s

  mobile_mole_fraction = rt_auxvar%total(idof_food_mobile,iphase)/sum_food
  immobile_mole_fraction = 1 - mobile_mole_fraction


  Residual(idof_food_mobile) = Residual(idof_food_mobile) + &
                           ! Growth usage, mol/s
                           this%stoichiometric_D_1* &
                           mobile_mole_fraction* &                                ! dimensionless
                           mu_B*material_auxvar%volume + &                        ! mol/m3 bulk/s * m3 bulk
                           ! Direct usage, mol/s
                           this%rate_D* &                          ! 1/s
                           mobile_mole_fraction* &                                ! dimensionless
                           rt_auxvar%immobile(this%B_id)* &                 ! mol/m3 bulk
                           material_auxvar%volume + &                             ! m3 bulk
                           ! Direct reaction, mol/s
                           this%stoichiometric_D_2* &
                           mobile_mole_fraction* &                                ! dimensionless
                           mu_CD*material_auxvar%volume* &                         ! mol/L/s * m3 bulk
                           immobile_to_water_vol + &                              ! L water/m3 bulk
                           ! immobilization, mol/s
                           this%rate_D_i* &                   ! 1/s
                           rt_auxvar%total(idof_food_mobile,iphase)* &            ! mol/L
                           material_auxvar%volume* &                              ! m3 bulk
                           immobile_to_water_vol - &                              ! L water/m3 bulk
                           ! remobilization, mol/s
                           this%rate_D_m* &                   ! 1/s
                           rt_auxvar%immobile(this%D_immobile_id)* &           ! mol/m3 bulk
                           material_auxvar%volume                                 ! m3 bulk


  Residual(idof_food_immobile) = Residual(idof_food_immobile) + &
                           ! Growth usage, mol/s
                           this%stoichiometric_D_1* &                          ! unitless
                           immobile_mole_fraction* &                              ! dimensionless
                           mu_B*material_auxvar%volume + &                        ! mol/m3 bulk/s * m3 bulk
                           ! Direct usage, mol/s
                           this%rate_D* &                          ! 1/s
                           immobile_mole_fraction* &                              ! dimensionless
                           rt_auxvar%immobile(this%B_id)* &                 ! mol/m3 bulk
                           material_auxvar%volume + &                             ! m3 bulk
                           ! Direct reaction, mol/s
                           this%stoichiometric_D_2* &
                           immobile_mole_fraction* &                              ! dimensionless
                           mu_CD*material_auxvar%volume* &                         ! mol/L/s * m3 bulk
                           immobile_to_water_vol - &                              ! L water/m3 bulk
                           ! immobilization, mol/s
                           this%rate_D_i* &                   ! 1/s
                           rt_auxvar%total(idof_food_mobile,iphase)* &            ! mol/L
                           material_auxvar%volume* &                              ! m3 bulk
                           immobile_to_water_vol + &                              ! L water/m3 bulk
                           ! remobilization, mol/s
                           this%rate_D_m* &                   ! 1/s
                           rt_auxvar%immobile(this%D_immobile_id)* &           ! mol/m3 bulk
                           material_auxvar%volume                                 ! m3 bulk

  Residual(idof_biocide) = Residual(idof_biocide) + &
                          this%mass_action_X* &                           ! L/mol/s
                          rt_auxvar%total(idof_biocide,iphase)* &                  ! mol/L
                          rt_auxvar%immobile(this%B_id)* &                   ! mol/m3 bulk
                          material_auxvar%volume                                   ! m3 bulk

  if (compute_derivative) then
    option%io_buffer = 'Analytical Jacobian calculation for this is a pain.' // &
      ' Use Numerical Jacobian only!!!!'
    call printErrMsg(option)
  endif

end subroutine ChromiumReact

! ************************************************************************** !

subroutine ChromiumKineticState(this,rt_auxvar,global_auxvar, &
                                               material_auxvar, &
                                               reaction,option)
  !
  ! Updates the kinetic state for the sandbox
  !
  ! Author: Satish Karra, Scott Hansen, Sachin Pandey, LANL
  ! Date: 09/15/2016
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module, only: material_auxvar_type

  implicit none

  class(reaction_sandbox_chromium_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  ! the following arrays must be declared after reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: biomass_residual_delta, delta_volfrac

  PetscReal :: mu_B
  PetscReal :: sum_food
  PetscInt :: idof_food_mobile, idof_food_immobile, idof_biomass, idof_Cr
  PetscInt :: idof_alcohol, idof_biocide
  PetscReal :: immobile_to_water_vol

  idof_food_mobile = this%D_mobile_id
  idof_Cr = this%C_id
  idof_alcohol = this%I_id
  idof_biocide = this%X_id
  idof_biomass = reaction%offset_immobile + this%B_id
  idof_food_immobile = reaction%offset_immobile + this%D_immobile_id

  immobile_to_water_vol = &
            material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0            ! L water/ m3 bulk

  sum_food = rt_auxvar%total(idof_food_mobile,iphase) + &
            rt_auxvar%immobile(this%D_immobile_id)/ &
            immobile_to_water_vol                                                 ! in mol/L water; Note that food_immobile is divided by porosity*saturation

  mu_B = this%rate_B_1*rt_auxvar%immobile(this%B_id)* &      ! mol/m3 bulk/s
            ! F monod term, unitless
            (sum_food/(sum_food + this%monod_D))* &
            ! B monod inhibition term, unitless
            (this%inhibition_B/ (rt_auxvar%immobile(this%B_id) + this%inhibition_B))**this%exponent_B * &
            ! I inhibition term, unitless
            (this%inhibition_I/ (rt_auxvar%total(idof_alcohol,iphase)+this%inhibition_I))

  biomass_residual_delta = &                                       ! Growth usage, mol/s
            - mu_B*material_auxvar%volume + &                      ! mol/m3 bulk/s * m3 bulk
            ! Natural decay, mol/s
            this%rate_B_2* &                         ! 1/s
            (rt_auxvar%immobile(this%B_id) - &
            this%background_concentration_B)* &                                   ! mol/m3 bulk
            material_auxvar%volume + &                             ! m3 bulk
            ! Biocide reaction, mol/s
            this%mass_action_B* &                 ! L/mol/s
            (rt_auxvar%immobile(this%B_id) - &
            this%background_concentration_B)* &                                   ! mol/m3 bulk
            rt_auxvar%total(idof_biocide,iphase)* &                ! mol/L
            material_auxvar%volume                                 ! m3 bulk

  delta_volfrac = &
            - biomass_residual_delta / &                           ! mol/s
            (this%density_B*1000.d0) / &                           ! mol/L * L/m3
            material_auxvar%volume * &                             ! m3 bulk
            option%tran_dt                                         ! s

  rt_auxvar%mnrl_volfrac(this%biomineral_id) = rt_auxvar%mnrl_volfrac(this%biomineral_id) + &
                                      delta_volfrac

end subroutine ChromiumKineticState

! ************************************************************************** !

subroutine ChromiumDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Satish Karra and Scott Hansen, LANL
  ! Date: 08/19/2015
  !

  implicit none

  class(reaction_sandbox_chromium_type) :: this

! 12. Add code to deallocate contents of the Cr_bio object

end subroutine ChromiumDestroy

end module Reaction_Sandbox_Chromium_class
