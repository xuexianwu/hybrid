!======================================================================!
SUBROUTINE NUDGE_WG
!----------------------------------------------------------------------!
! Forces generated monthly climate to equal input climate by adjusting
! daily values of temperature, radiation, vapour pressure, and
! precipitation proportionally.
!----------------------------------------------------------------------!
USE PHYSICAL_PARAMETERS_WG
USE GEN_ATMOS_WG
USE WGEN_params_WG
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
INTEGER :: ndm(nmo),kdayi,rf(nmo)
REAL :: rain_mo(nmo),rd_mo(nmo),raint(nmo),adj
REAL :: tmx_mo(nmo),tmn_mo(nmo),rg_mo(nmo),vap_mo(nmo)
REAL :: td,th,delta
REAL :: so(300)
!----------------------------------------------------------------------!
ndm = (/31,28,31,30,31,30,31,31,30,31,30,31/)
!----------------------------------------------------------------------!
!  Adjust rainfall so that has same number of wet days and the same
!  total rainfall as in climatology.
!----------------------------------------------------------------------!
kdayi = 1
raint(:) = 0.0
rd_mo (:) = 0.0
tmx_mo (:) = 0.0
tmn_mo (:) = 0.0
rg_mo (:) = 0.0
vap_mo (:) = 0.0
DO kmo = 1, nmo
  raint (kmo) = raint (kmo)+ pwet_mm (kmo) * rfrac_mm (kmo) * ndm (kmo)
  DO kday = 1, ndm (kmo)
    IF (rain (kdayi) > 0.0) rd_mo (kmo) = rd_mo (kmo) + 1.0
    tmx_mo (kmo) = tmx_mo (kmo) + dtmax (kdayi)
    tmn_mo (kmo) = tmn_mo (kmo) + dtmin (kdayi)
    rg_mo (kmo) = rg_mo (kmo) + rad (kdayi)
    vap_mo (kmo) = vap_mo (kmo) + vap (kdayi)
    kdayi = kdayi + 1
  ENDDO
ENDDO
!----------------------------------------------------------------------!
DO kmo = 1, nmo
  rf (kmo) = NINT (rd_mo (kmo) - rfrac_mm (kmo) * FLOAT (ndm (kmo)))
  tmx_mo (kmo) = tmx_mo (kmo) / FLOAT (ndm (kmo))
  tmn_mo (kmo) = tmn_mo (kmo) / FLOAT (ndm (kmo))
  rg_mo (kmo) = 0.04189 * rg_mo (kmo) / FLOAT (ndm (kmo))
  vap_mo (kmo) = vap_mo (kmo) / FLOAT (ndm (kmo))
ENDDO
!----------------------------------------------------------------------!
kdayi = 1
rain_mo (:) = 0.0
DO kmo = 1, nmo
  DO kday = 1, ndm (kmo)
    IF ((rf (kmo) < 0) .AND. (rain (kdayi) == 0.0)) THEN
      rain (kdayi) = pwet_mm (kmo)
      rf (kmo) = rf (kmo) + 1
    ENDIF
    IF ((rf (kmo) > 0) .AND. (rain (kdayi) > 0.0)) THEN
      rain (kdayi) = 0.0
      rf (kmo) = rf (kmo) - 1
    ENDIF
    rain_mo (kmo) = rain_mo (kmo) + rain (kdayi)
    dtmax (kdayi) = dtmax (kdayi) + (tmax_mm (kmo) * 9.0 / 5.0 + 32.0 &
    &               - tmx_mo (kmo))
    dtmin (kdayi) = dtmin (kdayi) + (tmin_mm (kmo) * 9.0 / 5.0 + 32.0 &
    &               - tmn_mo (kmo))
    IF (rg_mo (kmo) > 0.0) &
    & rad (kdayi) = rad (kdayi) * solrad_mm (kmo) / rg_mo (kmo)
    IF (vap_mo (kmo) > 0.0) &
    & vap (kdayi) = vap (kdayi) * dq_mm (kmo) / vap_mo (kmo)
    kdayi = kdayi + 1
  ENDDO
ENDDO
!----------------------------------------------------------------------!
kdayi = 1
DO kmo = 1, nmo
  IF (rain_mo (kmo) > 0.0) THEN
    adj = raint (kmo) / rain_mo (kmo)
  ELSE
    adj = 1.0
  ENDIF
  DO kday = 1, ndm (kmo)
    rain (kdayi) = adj * rain (kdayi)
    !------------------------------------------------------------------!
    kdayi = kdayi + 1
  ENDDO
ENDDO
!----------------------------------------------------------------------!
END SUBROUTINE NUDGE_WG
!======================================================================!
