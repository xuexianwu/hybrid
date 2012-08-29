subroutine pbl (ch_out,tg_in)
! Based on Ch of Hansen et al. (1983).
use atmos_variables
use soil_variables
use constants
implicit none
real*8 :: ch_out,tg_in
real*8, parameter :: zs = 10.0d0 ! Height of surface layer (m)
! Holder.
real*8, parameter :: z0 =  1.0d0  ! Roughness length (m)
real*8, parameter :: k  =  0.35d0
real*8, parameter :: a_ris  = 10.40d0
real*8, parameter :: b_ris  =  8.45d-1
real*8, parameter :: c_ris  =  1.68d0
real*8, parameter :: d_ris  =  0.81d0
real*8, parameter :: f_ris  =  0.14d0
real*8 :: ris,dm,cd,cdn
! Drag coefficient for neutral stability.
cdn = k * k / (log (zs / z0) * log (zs / z0))
! Bulk Richardson number.
ris = zs * grav * (ts - tg_in) / (tg_in * vsm * vsm)
if (ris < 0.0d0) then
  dm = (((1.0d0 - a_ris * ris) * (1.0d0 - b_ris * ris)) / &
  &    (1.0d0 - c_ris * ris)) ** 0.5d0
  cd = cdn * dm
  ch_out = cd * 1.35 * ((1.0d0 - d_ris * ris) / &
  &    (1.0d0 - f_ris * ris)) ** 0.5d0
else
  dm = 1.0d0 / (1.0d0 + (11.2d0 + 90.0d0 * ris) * ris)
  cd = cdn * dm
  ch_out = cd * 1.35 / (1.0d0 + 1.93 * ris)
endif
!write (*,*) 'pbl',ts-tg_in
end subroutine pbl
