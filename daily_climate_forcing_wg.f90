!======================================================================!
SUBROUTINE DAILY_CLIMATE_FORCING_WG
!----------------------------------------------------------------------!
! Generates sub-daily weather at tsteph times from daily values.
!----------------------------------------------------------------------!
use constants
USE GEN_ATMOS_WG
USE WGEN_PARAMS_WG
USE CONTROL_DEFINITIONS_WG
USE PHYSICAL_PARAMETERS_WG
use atmos_variables
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
!  Internal parameters.
!----------------------------------------------------------------------!
REAL                 :: td    ! Time (day)
REAL                 :: th    ! Solar time (hr)
REAL                 :: tmin  ! 24-hr min. temperature (degC)
REAL                 :: tmax  ! 24-hr max. temperature (degC)
REAL                 :: radMJ ! Daily radiation (MJ/m2/day)
REAL                 :: ratio ! Daily actual/potential radiation (ratio)
REAL                 :: delta ! Solar declination (rad)
REAL, DIMENSION(300) :: so    ! Solar flux with no atmosphere (W/m2)
REAL                 :: amp   ! Amplitude of daily radiation (W/m2)
REAL                 :: slope ! Slope of regression between radiation
                              ! and temperature (oC/(W/m2))
REAL                 :: intc  ! Intercept of regression between
                              ! radiation and temperature (oC)
REAL                 :: ipar  ! Downwelling total PAR (umol/m2/s)
REAL                 :: rf    ! Factor for separating dr/df radiation
REAL                 :: kl    ! Factor for separating dr/df radiation
REAL                 :: fdf   ! Fraction of diffuse PAR (fraction)
REAL*8               :: qsat  ! Saturation mixing ratio (kg/kg)
REAL*8               :: QSAT_WG  ! Saturation mixing ratio (kg/kg)
REAL                 :: cloud ! Cloud fraction (fraction)
real :: ld_mean ! Daily mean simulated Ld (W/m^2)
REAL :: Rchl,Rfrc
INTEGER, PARAMETER :: t1_doy = 305
REAL, PARAMETER :: Ccrit = 57.40
REAL, PARAMETER :: Fcrit = 262.53
REAL, PARAMETER :: TPmin = -3.8
REAL, PARAMETER :: TPopt =  0.41
REAL, PARAMETER :: TPmax = 12.43
REAL, PARAMETER :: TPb   =  0.0
REAL, PARAMETER :: aPH   =  26.49
REAL, PARAMETER :: bPH   =  -0.19
REAL, PARAMETER :: cPH   =  -17.07
real*8 :: emm
!----------------------------------------------------------------------!
!  Convert daily minimum temperature from oF to oC.
!----------------------------------------------------------------------!
tmin = (dtmin (jday) - 32.0) * 5.0 / 9.0
!----------------------------------------------------------------------!
!  Convert daily maximum temperature from oF to oC.
!----------------------------------------------------------------------!
tmax = (dtmax (jday) - 32.0) * 5.0 / 9.0
!----------------------------------------------------------------------!
!  Limit as found necessary for some simulations. May need to be
!  careful with this in future!
!----------------------------------------------------------------------!
tmax = MIN (tmax, 70.0)
!----------------------------------------------------------------------!
!  Ratio of radiation to potential for day (fraction).
!----------------------------------------------------------------------!
!IF (rad (kday) .GT. 500.0) rad (kday) = 0.8 * rad (kday)
!rad (kday) = 0.8 * rad (kday)
IF (rc (jday) .GT. 0.0) THEN
  ratio = rad (jday) / rc (jday)
ELSE
  ratio = 1.0
ENDIF
IF (ratio > 1.0) ratio = 1.0
!ratio=0.8 !****
!----------------------------------------------------------------------!
!  Convert daily radiation from Langleys/day to MJ/m2/day.
!----------------------------------------------------------------------!
radMJ = 0.04189 * rad (jday)
!----------------------------------------------------------------------!
!  Time (day).
!----------------------------------------------------------------------!
td = FLOAT (jday)
!----------------------------------------------------------------------!
!  Solar declination (rad). Taken from 'powerofsun'.
!----------------------------------------------------------------------!
delta = ASIN (0.39795 * COS (rfac * 0.98563 * (td - 173.0)))
!----------------------------------------------------------------------!
!  Initialise solar time for day (hr).
!----------------------------------------------------------------------!
th = 0.0
!----------------------------------------------------------------------!
!  Loop over tpoints in day.
!----------------------------------------------------------------------!
DO ktpoint = 1, ntpoints
!----------------------------------------------------------------------!
!  Sub-daily climate.
!----------------------------------------------------------------------!
!  Sine of solar elevation (ratio).
!----------------------------------------------------------------------!
  sbeta (ktpoint) = SIN (latr) * SIN (delta) + COS (latr) * COS (delta)&
  &                 * COS (rfac * 15.0 * (th - 12.0))
!----------------------------------------------------------------------!
!  Solar flux with no atmosphere (W/m2).
!----------------------------------------------------------------------!
  so (ktpoint) = scs * (1.0 + 0.034 *                                 &
  &              COS (rfac * 360.0 * td / 365.25)) * sbeta (ktpoint)
!----------------------------------------------------------------------!
!  Solar flux on earth's surface (W/m2).
!----------------------------------------------------------------------!
  rg (ktpoint) = MAX (ratio * so (ktpoint), 0.0)
  !--------------------------------------------------------------------!
  !  Add solar radiation to energy balance diagnostic (J).
  !--------------------------------------------------------------------!
  !ebal (:) = ebal (:) + tstep * rg (ktpoint)
!----------------------------------------------------------------------!
!  Downwelling total PAR (umol/m2/s).
!----------------------------------------------------------------------!
  ipar = 2.3 * rg (ktpoint)
!----------------------------------------------------------------------!
!  Separate direct and diffuse fractions.
!----------------------------------------------------------------------!
  rf = 0.847 - 1.61 * sbeta (ktpoint) + 1.04 * (sbeta (ktpoint) **2)
!----------------------------------------------------------------------!
  kl = (1.47 - rf) / 1.66
!----------------------------------------------------------------------!
  IF (ratio .LE. 0.22) THEN
    fdf = 1.0
  ELSEIF (ratio .LE. 0.35) THEN
    fdf = 1.0 - 6.4 * ((ratio - 0.22) **2)
  ELSEIF (ratio .LE. kl) THEN
    fdf = 1.47 - 1.66 * ratio
  ELSE
    fdf = rf
  ENDIF
!----------------------------------------------------------------------!
  ipar_dr (ktpoint) = (1.0 - fdf) * ipar
  ipar_df (ktpoint) = fdf * ipar
!----------------------------------------------------------------------!
!  Atmospheric pressure (Pa).
!----------------------------------------------------------------------!
  p (ktpoint) = pressure (jday)
!----------------------------------------------------------------------!
!  Wind speed (m/s). Fixed for now because varies wildly otherwise.
!----------------------------------------------------------------------!
  u (ktpoint) = windsp (jday)
!----------------------------------------------------------------------!
!  Update solar time for next time point (ht).
!----------------------------------------------------------------------!
  th = th + tsteph
!----------------------------------------------------------------------!
!  Next time point in day (ktpoint = 1, ntpoints).
!----------------------------------------------------------------------!
END DO
!----------------------------------------------------------------------!
!  Amplitude of daily radiation (W/m2).
!----------------------------------------------------------------------!
amp = rg ((ntpoints / 2) + 1) - rg (1)
!----------------------------------------------------------------------!
!  Slope of regression between radiation and temperature (oC/(W/m2)).
!----------------------------------------------------------------------!
IF (amp .GT. 0.0) THEN
  slope = (tmax - tmin) / amp
ELSE
  slope = 0.0
ENDIF
!----------------------------------------------------------------------!
!  Intercept of regression between radiation and temperature (oC).
!----------------------------------------------------------------------!
intc = tmin - slope * rg (1)
!----------------------------------------------------------------------!
!  Initialise solar time for day (hr).
!----------------------------------------------------------------------!
th = 0.0
!--------------------------------------------------------------------!
! Initialise Ld daily mean (W/m^2).
!--------------------------------------------------------------------!
ld_mean = 0.0d0
!----------------------------------------------------------------------!
!  Loop over time points in day to calculate tempereatures and Ld.
!----------------------------------------------------------------------!
DO ktpoint = 1, ntpoints
!----------------------------------------------------------------------!
!  Air temperature (oC).
!----------------------------------------------------------------------!
  ta (ktpoint) = intc + slope * rg (ktpoint)
!----------------------------------------------------------------------!
!  Air temperature (K).
!----------------------------------------------------------------------!
  tak (ktpoint) = ta (ktpoint) + tfrz
!----------------------------------------------------------------------!
!  Precipitation spread through day.
!----------------------------------------------------------------------!
  ppt (ktpoint) = 1.0e-3 * rain (jday) / (60.0 * 60.0 * tsteph * &
  &               float (ntpoints))
  !ppt (ktpoint) = 0.0d0
  !if (ktpoint == 12) then
  !  ppt (ktpoint) = 1.0d-3 * rain (jday) / (3600.0d0 * tsteph)
  !endif
!----------------------------------------------------------------------!
!  Generated VAP for each day is in kPa. Convert to
!  mixing ratio (kg/kg) if necessary.
!----------------------------------------------------------------------!
  if (isi) then
    qa (ktpoint) = vap (jday)
  else
    qa (ktpoint) = 0.622e3 * vap (jday) / (p (ktpoint) - &
    &              1.0e3 * vap (jday))
  endif
  !--------------------------------------------------------------------!
  ! Calculate qsat (kg/kg) from t (degC) and p (mb).
  !--------------------------------------------------------------------!
  qsat = QSAT_WG (ta (ktpoint), p (ktpoint))
  !--------------------------------------------------------------------!
  IF (qa (ktpoint) .GT. qsat) qa (ktpoint) = qsat
  !--------------------------------------------------------------------!
  !  Cloud cover (fraction).
  !--------------------------------------------------------------------!
  cloud = 1.0d0 - ratio
  cloud = min (1.0d0,cloud)
  cloud = max (0.0d0,cloud)
  !--------------------------------------------------------------------!
  !  Downwelling longwave radiation (W/m2). Uses formulation from
  !  Pirazzini with eclr = 0.69.
  !--------------------------------------------------------------------!
  !ld (ktpoint) = (0.69 * (1.0 - (cloud ** 6)) + 0.979 * (cloud ** 4))&
  !&    * stblz * (tak (ktpoint) ** 4)
  emm = 0.787812 * (1.0d0 - (cloud ** 6)) + 0.971812 * (cloud ** 4)
  ld (ktpoint) = emm * stblz * (tak (ktpoint) ** 4)
  !--------------------------------------------------------------------!
  ! Downwelling longwave radiation from Dilley & O'Brien (1998) (W/m^2).
  !--------------------------------------------------------------------!
  !ld (ktpoint) = 59.38d0 + (113.7d0 * ((ta (ktpoint)+ tf)/tf) ** 6)  &
  !&  + 96.96d0 * ((465.0d1 * vap (jday) / (ta (ktpoint) + tf)) /     &
  !& 25.0d0) ** 0.5
  !--------------------------------------------------------------------!
  ! Sum for Ld daily mean (W/m^2).
  !--------------------------------------------------------------------!
  ld_mean = ld_mean + ld (ktpoint)
  !--------------------------------------------------------------------!
  !  Update solar time for next time point (hr)
  !--------------------------------------------------------------------!
  th = th + tsteph
!----------------------------------------------------------------------!
ENDDO
!----------------------------------------------------------------------!
! Adjust sub-daily Ld values to give desired daily mean (W/m^2).
!----------------------------------------------------------------------!
ld_mean = ld_mean / float (ntpoints)
if (ld_mean > 0.0d0) then
  ld (:) = ld (:) * ldi_ij (i0,j0) / ld_mean
else
  ld (:) = ldi_ij (i0,j0)
endif
!----------------------------------------------------------------------!
RETURN
END SUBROUTINE DAILY_CLIMATE_FORCING_WG
!======================================================================!
FUNCTION QSAT_WG (tc,p)
!----------------------------------------------------------------------!
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
!@ver  1.0
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
!@var A,B,C   expansion coefficients for QSAT
REAL, PARAMETER :: A = 3.797915e0    !3.797915d0
REAL, PARAMETER :: B = 7.93252e-6    !7.93252d-6
REAL, PARAMETER :: C = 2.166847e-3   !2.166847d-3
REAL*8 :: tc
REAL*8 :: p
REAL*8 :: QL
REAL*8 :: TM
REAL*8 :: PR
REAL*8 :: QSAT_WG
!----------------------------------------------------------------------!
!**** Note that if QL is considered to be a function of temperature, the
!**** correct argument in QSAT is the average QL from t=0 (C) to TM, ie.
!**** QL = 0.5*(QL(0)+QL(t))
!      REAL*8, INTENT(IN) :: TM  !@var TM   potential temperature (K)
!      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
!      REAL*8, INTENT(IN) :: QL  !@var QL   lat. heat of vap./sub. (J/kg)
!      REAL*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
QL = 2500800.0 - 2360.0 * tc
TM = tc + 273.16
PR = p / 100.0
QSAT_WG = A * EXP (QL * (B - C / TM)) / PR
!----------------------------------------------------------------------!
RETURN
END
!======================================================================!
