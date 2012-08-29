!======================================================================!
subroutine simveg
!----------------------------------------------------------------------!
use veg_variables
use constants
use control_parameters
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
 ra         = cland (1) * 0.5d0 / (14.0d0 * sday)
 cgrow      = cland (1) * 0.5d0 / (14.0d0 * sday)
 clit       = cland (2)         / (13.0d0 * 365.0d0 * sday)
 clit_trans = cland (3)         / ( 0.5d0 * 365.0d0 * sday)
 rh_clit    = cland (3)         / ( 0.5d0 * 365.0d0 * sday)
 rh_som     = cland (4)         / ( 9.0d0 * 365.0d0 * sday)
!----------------------------------------------------------------------!
! Diagnostics.
!----------------------------------------------------------------------!
npp = gpp - ra
nee = rh_clit + rh_som - npp
!----------------------------------------------------------------------!
end subroutine simveg
!======================================================================!
subroutine veg_conductance
!----------------------------------------------------------------------!
! Local canopy conductance.
!----------------------------------------------------------------------!
use veg_parameters
use veg_variables
use soil_variables
use constants
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
if (physiol) then
  !--------------------------------------------------------------------!
  call canopy
  !--------------------------------------------------------------------!
else
  !--------------------------------------------------------------------!
  call cond
  !--------------------------------------------------------------------!
endif
!----------------------------------------------------------------------!
end subroutine veg_conductance
!======================================================================!
subroutine cond
!----------------------------------------------------------------------!
! Local canopy conductance.
! Outputs are:
! cnc
! gpp
! veg_ci
!----------------------------------------------------------------------!
use constants
use soil_variables
use veg_variables
use atmos_variables
use gen_atmos_wg
use control_definitions_wg
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
real*8 :: dqs
!----------------------------------------------------------------------!
rt = gasc * (tp (0,2) + tf)
!----------------------------------------------------------------------!
! Total canopy transmittance of shortwave (fraction).
! Simple for now.
!----------------------------------------------------------------------!
trans_sw = exp (-0.5d0 * alai)
fpart (:) = 1.0d0 - exp (-0.5d0 * alait (:))
!----------------------------------------------------------------------!
! Compute humidity deficit across canopy surface and conductance factor.
!----------------------------------------------------------------------!
dqs = qv - qf
if (dqs < zero) dqs = zero
f_dq = 2.8d0 ** (-80.0d0 * dqs)
!----------------------------------------------------------------------!
! Compute cnc (m/s) and gpp (umol/m2/s).
!----------------------------------------------------------------------!
call zbrent_simple (veg_ci)
!----------------------------------------------------------------------!
! Convert internal CO2 from Pa to mol/m^3.
!----------------------------------------------------------------------!
veg_ci = veg_ci / rt
!----------------------------------------------------------------------!
end subroutine cond
!======================================================================!
subroutine zbrent_simple (b)
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
integer, parameter :: itmax = 100 ! ZBRENT loop max. (n)
real*8, parameter :: tol = 0.1 ! ZBRENT solution accuracy (umol/m^2/s)
real*8, parameter :: x1 = 0.001  ! Min. ci (Pa)
real*8, parameter :: x2 = 1.0e20 ! Max. ci (Pa)
real*8, parameter :: eps = 1.0e-12
real*8 :: afunc_simple
integer :: iter
real*8 :: a,b,c,d,e,fa,fb,fc,p,q,ri,s,tol1,xm
!----------------------------------------------------------------------!
a=x1
b=x2
fa=afunc_simple(a)
fb=afunc_simple(b)
if((fa > 0.0 .and. fb > 0.0).or.(fa < 0.0.and.fb < 0.0)) then
  write (*,*)a,b
  write (*,*)fa,fb
  write (*,*) 'root must be bracketed for zbrent_simple'
  write (99,*)a,b
  write (99,*)fa,fb
  write (99,*) 'root must be bracketed for zbrent_simple'
  stop
endif
c=b
fc=fb
do iter = 1, itmax
  if((fb>0.0.and.fc>0.0).or.(fb<0.0.and.fc<0.0))then
    c=a
    fc=fa
    d=b-a
    e=d
  endif
  if(abs(fc)<abs(fb))then
    a=b
    b=c
    c=a
    fa=fb
    fb=fc
    fc=fa
  endif
  tol1=2.0*EPS*abs(b)+0.5*tol
  xm=0.5*(c-b)
  if(abs(xm)<=tol1 .or. fb==0.0) return
  if(abs(e)>=tol1 .and. abs(fa)>abs(fb)) then
    s=fb/fa
    if(a==c)then
      p=2.0*xm*s
      q=1.0-s
    else
      q=fa/fc
      ri=fb/fc
      p=s*(2.0*xm*q*(q-ri)-(b-a)*(ri-1.0))
      q=(q-1.0)*(ri-1.0)*(s-1.0)
    endif
    if(p>0.0) q=-q
    p=abs(p)
    if(2.0*p < min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
      e=d
      d=p/q
    else
      d=xm
      e=d
    endif
  else
    d=xm
    e=d
  endif
  a=b
  fa=fb
  if(abs(d)>tol1) then
    b=b+d
  else
    b=b+sign(tol1,xm)
  endif
  fb=afunc_simple(b)
enddo
write (*,*) 'zbrent_simple exceeding maximum iterations'
write (99,*) 'zbrent_simple exceeding maximum iterations'
stop
end subroutine zbrent_simple
!======================================================================!
function afunc_simple (ci)
!----------------------------------------------------------------------!
use constants
use soil_variables
use veg_variables
use atmos_variables
use gen_atmos_wg
use control_definitions_wg
use control_parameters
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
real*8, parameter :: eps = 1.0d-12
real*8, parameter :: c1 = 90.0d0
real*8 :: afunc_simple,ci,acan_sat,srht0,rcan,acan_bio
real*8 :: f_ci,rtot,acan_gcan
 !--------------------------------------------------------------------!
 ! CO2 partial pressure at surface level (Pa).
 !--------------------------------------------------------------------!
 capa = co2_ppm * ps / 1.0d4
 !--------------------------------------------------------------------!
 cnc = f_swp * alaie / rs
 srht0 = max (srheat,zero)
 cnc = cnc * (srht0 / c1) / (1.0d0 + srht0 / c1)
 cnc = cnc / (1.0d0 + ((tp (0,2) + tf - 296.d0) / 15.0d0) ** 4)
 acan_sat = 1.8d0 * cnc / 650.0d-6 ! really need to calib. to canopy
 rcan = 0.1d0 * acan_sat ! improve
 acan_bio = acan_sat * ci / (23.0d0 + ci) - rcan
 f_ci = (ci + 9.6) / (ci + eps)
 cnc = 650.0e-6 * acan_sat * f_ci * f_dq * f_swp
 !  if (mod (itime+nday/2, nday) == 0) then
 !write (*,*) jday
 !write (*,*) betad,alai,1.0e3*cnc
 !endif
 cnc = max (cnc,alai*0.0001d0) ! need to improve
 cnc = min (cnc,alai*0.005d0) ! need to improve
 !  if (mod (itime+nday/2, nday) == 0) then
 !write (*,*) betad,alai,1.0e3*cnc
 !endif
 rtot = 1.0 / (cnc + eps) + 1.0 / (ch * vsm + eps)
 acan_gcan = (capa - ci) * 1.0d6 / (1.6 * rt * rtot + eps)
 afunc_simple = acan_bio - acan_gcan
!----------------------------------------------------------------------!
! kgC/m^2/s
!----------------------------------------------------------------------!
gpp = 0.012d0 * (acan_bio + rcan) / 1.0d6
!----------------------------------------------------------------------!
end function afunc_simple
!======================================================================!
subroutine update_veg_locals
!----------------------------------------------------------------------!
use atmos_variables
use constants
use veg_variables
use soil_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
real*8 :: evap_tot2
!----------------------------------------------------------------------!
 cna = ch * vsm
!----------------------------------------------------------------------!
! New water vapour mixing ratio for next timestep (kg/kg).
!----------------------------------------------------------------------!
qf = (fd * evapvd + fw * evapvw) / (rho3 * cna + 1.0d-12) + qs
!----------------------------------------------------------------------!
end subroutine update_veg_locals
!======================================================================!
