module Auxiliary_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Global_Aux_module
  use TH_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module
  use NW_Transport_Aux_module
  use Mphase_Aux_module
  use General_Aux_module
  use Hydrate_Aux_module
  use WIPP_Flow_Aux_module
  use Material_Aux_module
  use ERT_Aux_module
  use ZFlow_Aux_module
  use PNF_Aux_module
  use Secondary_Continuum_Aux_module
  use InlineSurface_Aux_module
  use Inversion_TS_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: auxiliary_type
    type(global_type), pointer :: Global
    type(reactive_transport_type), pointer :: RT
    type(nw_transport_type), pointer :: NWT
    type(th_type), pointer :: TH
    type(richards_type), pointer :: Richards
    type(zflow_type), pointer :: ZFlow
    type(pnf_type), pointer :: PNF
    type(mphase_type), pointer :: Mphase
    type(general_type), pointer :: General
    type(hydrate_type), pointer :: Hydrate
    type(wippflo_type), pointer :: WIPPFlo
    type(material_type), pointer :: Material
    type(ert_type), pointer :: ERT
    type(sc_heat_type), pointer :: SC_heat
    type(sc_rt_type), pointer :: SC_RT
    type(inlinesurface_type), pointer :: InlineSurface
    type(inversion_forward_aux_type), pointer :: inversion_forward_aux
  end type auxiliary_type

  public :: AuxInit, &
            AuxDestroy

contains

! ************************************************************************** !

subroutine AuxInit(aux)
  !
  ! Nullifies pointers in auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 04/09/08
  !

  implicit none

  type(auxiliary_type) :: aux

  nullify(aux%Global)
  nullify(aux%RT)
  nullify(aux%NWT)
  nullify(aux%TH)
  nullify(aux%Richards)
  nullify(aux%ZFlow)
  nullify(aux%PNF)
  nullify(aux%ERT)

  nullify(aux%Mphase)
  nullify(aux%General)
  nullify(aux%Hydrate)
  nullify(aux%WIPPFlo)
  nullify(aux%Material)
  nullify(aux%SC_heat)
  nullify(aux%SC_RT)
  nullify(aux%InlineSurface)
  nullify(aux%inversion_forward_aux)

end subroutine AuxInit

! ************************************************************************** !

subroutine AuxDestroy(aux)
  !
  ! Deallocates any allocated pointers in auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 04/09/08
  !

  implicit none

  type(auxiliary_type) :: aux

  call GlobalAuxDestroy(aux%Global)
  call RTAuxDestroy(aux%RT)
  call NWTAuxDestroy(aux%NWT)
  call THAuxDestroy(aux%TH)
  call RichardsAuxDestroy(aux%Richards)
  call ZFlowAuxDestroy(aux%ZFlow)
  call PNFAuxDestroy(aux%PNF)
  call MphaseAuxDestroy(aux%Mphase)
  call GeneralAuxDestroy(aux%General)
  call HydrateAuxDestroy(aux%Hydrate)
  call WIPPFloAuxDestroy(aux%WIPPFlo)
  call MaterialAuxDestroy(aux%Material)
  call ERTAuxDestroy(aux%ERT)
  call SecondaryAuxHeatDestroy(aux%SC_heat)
  call SecondaryAuxRTDestroy(aux%SC_RT)
  call InlineSurfaceAuxDestroy(aux%InlineSurface)
  ! DO NOT destroy aux%inversion_forward_aux; it is destroyed elsewhere

  nullify(aux%Global)
  nullify(aux%RT)
  nullify(aux%NWT)
  nullify(aux%Richards)
  nullify(aux%ZFlow)
  nullify(aux%PNF)
  nullify(aux%Mphase)
  nullify(aux%General)
  nullify(aux%Hydrate)
  nullify(aux%WIPPFlo)
  nullify(aux%Material)
  nullify(aux%ERT)
  nullify(aux%SC_Heat)
  nullify(aux%SC_RT)
  nullify(aux%InlineSurface)
  nullify(aux%inversion_forward_aux)

end subroutine AuxDestroy

end module Auxiliary_module
