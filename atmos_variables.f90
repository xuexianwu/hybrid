!======================================================================!
module atmos_variables
!----------------------------------------------------------------------!
! Definitions for all atmosphere-related variables.
!----------------------------------------------------------------------!
use control_parameters
use constants
!----------------------------------------------------------------------!
implicit none
save
!----------------------------------------------------------------------!
! Global climate fields.
!---------------------------------------------------------------------!
! Precipitation (mm/day).
!---------------------------------------------------------------------!
real*8, allocatable :: raini_ij    (:,:)
!---------------------------------------------------------------------!
! Min. 24-hr temperature (degC).
!---------------------------------------------------------------------!
real*8, allocatable :: tmin_ij     (:,:)
!---------------------------------------------------------------------!
! Max. 24-hr temperature (degC).
!---------------------------------------------------------------------!
real*8, allocatable :: tmax_ij     (:,:)
!---------------------------------------------------------------------!
! Atmos. vapour pressure (hPa or g/kg > kg/kg (isi)).
!---------------------------------------------------------------------!
real*8, allocatable :: vap_ij      (:,:)
!---------------------------------------------------------------------!
! Solar SW radiation (W/m^2 > MJ/m^2/d).
!---------------------------------------------------------------------!
real*8, allocatable :: radi_ij     (:,:)
!---------------------------------------------------------------------!
! Maximum solar radiation (MJ/m^2/d).
!---------------------------------------------------------------------!
real*8, allocatable :: rc_ij       (:,:)
!---------------------------------------------------------------------!
! Wind speed (m/s).
!---------------------------------------------------------------------!
real*8, allocatable :: windsp_ij   (:,:)
!---------------------------------------------------------------------!
! Surface atmos. pressure (MPa or kPa (isi)).
!---------------------------------------------------------------------!
real*8, allocatable :: pressure_ij (:,:)
!---------------------------------------------------------------------!
! Ld radiation (W/m^2).
!---------------------------------------------------------------------!
real*8, allocatable :: ldi_ij     (:,:)
!----------------------------------------------------------------------!
! Total precipitation mass in timepoint dtpoint (kg/m^2).
!----------------------------------------------------------------------!
real*8, allocatable :: prec        (:,:)
!----------------------------------------------------------------------!
! Atmosphere surface level pressure at timepoint dtpoint (hPa).
!----------------------------------------------------------------------!
real*8, allocatable :: pedn        (:,:)
!----------------------------------------------------------------------!
! Atmosphere surface level temperature at timepoint dtpoint (K).
!----------------------------------------------------------------------!
real*8, allocatable :: t           (:,:)
!----------------------------------------------------------------------!
! Surface SW radiation at timstep dtpoint (W/m^2).
!----------------------------------------------------------------------!
real*8, allocatable :: fsf         (:,:)
!----------------------------------------------------------------------!
! Mass of air in surface level (kg/m^2).
!----------------------------------------------------------------------!
real*8, allocatable :: am          (:,:)
!----------------------------------------------------------------------!
! Atmosphere surface level specific humidity at timestep (kg/kg).
!----------------------------------------------------------------------!
real*8, allocatable :: q_surf      (:,:)
!----------------------------------------------------------------------!
! Surface LWd radiation at timstep dtpoint (W/m^2).
!----------------------------------------------------------------------!
real*8, allocatable :: trhr        (:,:)
!----------------------------------------------------------------------!
! Sine of solar elevation (ratio).
!----------------------------------------------------------------------!
real*8, allocatable :: sbeta_can   (:,:)
!----------------------------------------------------------------------!
! Downwelling diffuse PAR (umol/m^2/s).
!----------------------------------------------------------------------!
real*8, allocatable :: ipar_df_can   (:,:)
!----------------------------------------------------------------------!
! Downwelling direct PAR (umol/m^2/s).
!----------------------------------------------------------------------!
real*8, allocatable :: ipar_dr_can   (:,:)
!----------------------------------------------------------------------!
! Local sine of solar elevation (ratio).
!----------------------------------------------------------------------!
real*8 :: sbeta_l
!----------------------------------------------------------------------!
! Local downwelling diffuse PAR (umol/m^2/s).
!----------------------------------------------------------------------!
real*8 :: ipar_df_l
!----------------------------------------------------------------------!
! Local downwelling direct PAR (umol/m^2/s).
!----------------------------------------------------------------------!
real*8 :: ipar_dr_l
!----------------------------------------------------------------------!
! Local precipitation rate (m/s == 10^3 kg/m2/s).
!----------------------------------------------------------------------!
real*8 :: pr
!----------------------------------------------------------------------!
! Heat in precipitation (W/m^2).
!----------------------------------------------------------------------!
real*8 :: htpr
!----------------------------------------------------------------------!
! Local fraction (by area) of precipitation (fraction).
!----------------------------------------------------------------------!
real*8 :: prfr
!----------------------------------------------------------------------!
! Local atmospheric surface layer pressure (hPa).
!----------------------------------------------------------------------!
real*8 :: ps
!----------------------------------------------------------------------!
! Local incident SW (W/m^2).
!----------------------------------------------------------------------!
real*8 :: srheat
real*8 :: srht
!----------------------------------------------------------------------!
! Local incident LW (W/m^2).
!----------------------------------------------------------------------!
real*8 :: trheat
real*8 :: trht
!----------------------------------------------------------------------!
! Local atmosphere surface level temperature (K).
!----------------------------------------------------------------------!
real*8 :: ts
!----------------------------------------------------------------------!
! Local atmosphere surface level sp. humidity (kg/kg).
!----------------------------------------------------------------------!
real*8 :: qs
!----------------------------------------------------------------------!
! Local surface atmosphere level density (kg/m^3).
!----------------------------------------------------------------------!
real*8 :: rho
!----------------------------------------------------------------------!
! Air density at surface height (10^3 kg/m^3).
!----------------------------------------------------------------------!
real*8 :: rho3
!----------------------------------------------------------------------!
! Mass of water in first atmosphere layer (kg/m^2).
!----------------------------------------------------------------------!
real*8 :: qm1
!----------------------------------------------------------------------!
! Heat drag coefficient (dimensionless).
!----------------------------------------------------------------------!
real*8 :: ch
!----------------------------------------------------------------------!
! Wind speed at surface height (m/s).
!----------------------------------------------------------------------!
real*8 :: vsm
!----------------------------------------------------------------------!
! Mean global atmospheric co2 mixing ratio equivalent (ppm).
!----------------------------------------------------------------------!
real*8 :: co2eq_ppm
!----------------------------------------------------------------------!
! Mean global atmospheric co2 mixing ratio Kyoto equivalent (ppm).
!----------------------------------------------------------------------!
real*8 :: kyoto_co2eq_ppm
!----------------------------------------------------------------------!
! Mean global atmospheric co2 mixing ratio (ppm).
!----------------------------------------------------------------------!
real*8 :: co2_ppm
!----------------------------------------------------------------------!
! Used in evap_limits.
!----------------------------------------------------------------------!
!real*8, allocatable :: evap_max_ij (:,:)
!real*8, allocatable :: fr_sat_ij (:,:)
!----------------------------------------------------------------------!
end module atmos_variables
!======================================================================!
