!======================================================================!
MODULE GEN_ATMOS_WG
!----------------------------------------------------------------------!
use control_parameters
!----------------------------------------------------------------------!
IMPLICIT NONE
save
!----------------------------------------------------------------------!
!  Daily generated atmosphere properties.
!----------------------------------------------------------------------!
REAL, DIMENSION(366) :: rain     ! Precipitation (mm/day)
REAL, DIMENSION(366) :: dtmin    ! Daily minimum temperature (oF)
REAL, DIMENSION(366) :: dtmax    ! Daily maximum temperature (oF)
REAL, DIMENSION(366) :: rad      ! Daily radiation (Langleys/day)
REAL, DIMENSION(366) :: vap      ! Mean daily vapour pressure (MPa)
REAL, DIMENSION(366) :: windsp   ! Mean daily windspeed (m/s)
real, dimension(366) :: pressure ! Mean daily surface pressure (Pa)
real*8, allocatable :: rain_ij (:,:,:) ! Precipitation    (mm/s)
real*8, allocatable :: pedn_ij (:,:,:) ! Surface pressure  (hPa)
real*8, allocatable :: tak_ij  (:,:,:) ! Surface temperautre  (K)
real*8, allocatable :: qsrf_ij (:,:,:) ! Surface sp. hum  (kg/kg)
real*8, allocatable :: rad_ij  (:,:,:) ! SW radiation  (W/m^2)
real*8, allocatable :: ld_ij   (:,:,:) ! LW radiation  (W/m^2)
!---------------------------------------------------------------------!
! sin of solar elevation (ratio).
!---------------------------------------------------------------------!
real*8, allocatable :: sbeta_ij   (:,:,:)
!---------------------------------------------------------------------!
! Downwelling diffuse irradiance (umol(PAR)/m^2/s).
!---------------------------------------------------------------------!
real*8, allocatable :: ipar_df_ij (:,:,:)
!---------------------------------------------------------------------!
! Downwelling direct irradiance (umol(PAR)/m^2/s).
!---------------------------------------------------------------------!
real*8, allocatable :: ipar_dr_ij (:,:,:)
!----------------------------------------------------------------------!
!  Sub-daily generated atmosphere properties.
!----------------------------------------------------------------------!
real*8, allocatable :: rg (:)    ! Solar flux at earth's sfc    (W/m2)
real*8, allocatable :: ld (:)    ! Downwelling longwave rad     (W/m2)
real*8, allocatable :: ta (:)    ! Air temperature                (oC)
real*8, allocatable :: qa  (:)   ! Atmos. water vap. mix. rat. (kg/kg)
real*8, allocatable :: ppt (:)   ! Precipitation                 (m/s)
real*8, allocatable :: p   (:)   ! Atmospheric pressure           (Pa)
real*8, allocatable :: tak (:)   ! Air temperature                 (K)
real*8, allocatable :: u (:)       ! Wind speed                    (m/s)
real*8, allocatable :: sbeta (:)   ! Sine of solar elevation     (ratio)
real*8, allocatable :: ipar_dr (:) ! Downwelling direct PAR
                                   ! above canopy            (umol/m2/s)
real*8, allocatable :: ipar_df (:) ! Downwelling diffuse PAR
                                   ! above canopy            (umol/m2/s)
!----------------------------------------------------------------------!
END MODULE GEN_ATMOS_WG
!======================================================================!
