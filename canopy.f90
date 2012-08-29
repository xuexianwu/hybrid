!======================================================================!
subroutine canopy
!----------------------------------------------------------------------!
! Computes:
! cnc
! gpp
! veg_ci
!----------------------------------------------------------------------!
use control_parameters
use veg_parameters
use veg_variables
use soil_variables
use atmos_variables
use gen_atmos_wg
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
integer :: nlayer
real*8 :: ftmp,ntop,lcu,dqs
!----------------------------------------------------------------------!
rt = gasc * (tp (0,2) + tf)
!----------------------------------------------------------------------!
! Compute humidity deficit across canopy surface and conductance factor.
!----------------------------------------------------------------------!
dqs = qv - qf
if (dqs < zero) dqs = zero
f_dq = 2.8d0 ** (-80.0d0 * dqs)
!----------------------------------------------------------------------!
! CO2 partial pressure at surface level (Pa).
!----------------------------------------------------------------------!
  capa = co2_ppm * ps / 1.0d4
!----------------------------------------------------------------------!
n_layers = nint (alai / laid)
n_layers = max (1, n_layers)
n_layers = min (ncanopy_layers_max,n_layers)
!----------------------------------------------------------------------!
 c3 = .true.
 !c3 = .false.
!----------------------------------------------------------------------!
if (c3) then
  !--------------------------------------------------------------------!
  ! C3 canopy mitochondrial respiration (umol/m^2/s).
  !--------------------------------------------------------------------!
  rpcrn = 0.2e3 * 14.0e-6 * nfp * ncan * exp (18.72 - 46.39e3 / rt)
  !--------------------------------------------------------------------!
else
  !--------------------------------------------------------------------!
  ! C4 canopy mitochondrial respiration (umol/m^2/s).
  !--------------------------------------------------------------------!
  ftmp = 2.0 ** ((tp (0,2) - 25.0) / 10.0)
  rpcrn = 4.0e-3 * nfp * ncan * ftmp / &
  &       (1.0 + exp (1.3 * (tp (0,2) - 55.0)))
  !--------------------------------------------------------------------!
endif
!----------------------------------------------------------------------!
ntop = kn * ncan / (1.0 - exp (-kn * alai) + eps_sing)
lcu = laid / 2.0
do layer = 1, n_layers
  nl (layer) = exp (-kn * lcu) * ntop
  if (c3) then
    n3 (layer) = 6.0 - 3.6 * exp (-kn3 * lcu)
  else
    n3 (layer) = 0.35 * (6.0 - 3.6 * exp (-kn3 * lcu))
  endif
  fabsb (layer) = exp (-ka * n3 (layer) * nl (layer))
  lcu = lcu + laid
enddo
!----------------------------------------------------------------------!
pcp = 1.0e-4 * exp (19.02 - 37.83e3 / rt ) * ps ! Pa
kc  = 1.0e-4 * exp (38.05 - 79.43e3 / rt ) * ps ! Pa
ko  = 1.0e-4 * exp (20.30 - 36.38e3 / rt ) * ps ! kPa
km = kc * (1.0 + 20.9 / ko)
n1 = nfp
n2 = nfp
if (c3) then
  !--------------------------------------------------------------------!
  n1 = 0.12 * 1.3 * n1
  n1 = f_swp * exp (-((tp (0,2) - topt_c3) / omega) ** 2) * n1
  !--------------------------------------------------------------------!
  n2 = 0.23 * f_swp * exp (26.35 - 65.33e3 / rt) * n2
  !--------------------------------------------------------------------!
else
  !--------------------------------------------------------------------!
  n1 = 0.23 * n1 * f_swp * exp  (-((tp (0,2) - topt_c4) / omega) ** 2)
  !--------------------------------------------------------------------!
  n2 = 0.10 * f_swp * exp (22.55 - 55.9e3 / rt) * n2
  !--------------------------------------------------------------------!
endif
!----------------------------------------------------------------------!
call zbrent (veg_ci)
!----------------------------------------------------------------------!
! kgC/m^2/s
!----------------------------------------------------------------------!
gpp = 0.012d0 * (acann + rpcrn) / 1.0d6
!----------------------------------------------------------------------!
! Convert internal CO2 from Pa to mol/m^3.
!----------------------------------------------------------------------!
veg_ci = veg_ci / rt
!----------------------------------------------------------------------!
end subroutine canopy
!======================================================================!
subroutine zbrent (b)
!----------------------------------------------------------------------!
use control_parameters
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
integer, parameter :: itmax = 100 ! ZBRENT loop max. (n)
real*8, parameter :: tol = 0.1 ! ZBRENT solution accuracy (umol/m^2/s)
real*8, parameter :: x1 = 0.001  ! Min. ci (Pa)
real*8, parameter :: x2 = 1.0e20 ! Max. ci (Pa)
real*8 :: afunc,a,b,c,d,e,fa,fb,fc,p,q,ri,s,tol1,xm
integer :: iter
!----------------------------------------------------------------------!
a=x1
b=x2
fa=afunc(a)
fb=afunc(b)
if((fa>0.0.and.fb>0.0).or.(fb<0.0.and.fb.lt.0.0)) then
  write (*,*)a,b
  write (*,*)fa,fb
  write (*,*) 'root must be bracketed for zbrent'
  write (99,*)a,b
  write (99,*)fa,fb
  write (99,*) 'root must be bracketed for zbrent'
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
  tol1=2.0*eps_sing*abs(b)+0.5*tol
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
  fb=afunc(b)
enddo
write (99,*) 'zbrent exceeding maximum iterations'
write (99,*) 'zbrent exceeding maximum iterations'
stop
!----------------------------------------------------------------------!
end subroutine zbrent
!======================================================================!
function afunc (co2i)
!----------------------------------------------------------------------!
use control_parameters
use veg_parameters
use gen_atmos_wg
use soil_variables
use veg_variables
use atmos_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
integer :: ice
real*8 :: afunc
real*8 :: co2i,co2
real*8 :: temp1,temp2,radf1,radf2,radf3,radf4,radf3b
real*8 :: acang,lcu,acanni,rtot,f_ci
real*8 :: isha,fsl,isla,idfa,idra,idrdra,ash,asl,acangl,ansat,acann_save
real*8 :: rhopar ! PAR canopy reflection (fraction)
real*8 :: kblp   ! Ext. coeff for PAR over black leaves (unitless)
!----------------------------------------------------------------------!
do ice = 1, 2
  if (ice == 1) co2 = co2i
  if (ice == 2) co2 = 1.0e6
  if (sbeta_l > 0.0) then
    !------------------------------------------------------------------!
    radf3b = 1.0 - sigma
    temp1 = sqrt (radf3b)
    temp2 = 2.0 * (1.0 - temp1) / (1.0 + temp1)
    rhopar = temp2 / (1.0 + 1.6 * sbeta_l)
    kblp   = 0.5 * kdf / (0.8 * temp1 * sbeta_l)
    radf1 = (1.0 - rhopar) * ipar_df_l * kdf
    radf2 = ipar_dr_l * kblp
    radf3 = radf3b * radf2
    radf2 = (1.0 - rhopar) * temp1 * radf2
    radf4 = -temp1 * kblp
    !------------------------------------------------------------------!
    acang = 0.0
    lcu = laid / 2.0
    !------------------------------------------------------------------!
    do layer = 1, n_layers
      if (radf2 == 0.0) then ! No direct light
        isha = radf1 * exp (-kdf * lcu)
        fsl = 0.0
        isla = 0.0
      else ! Is direct light
        fsl = exp (-kblp * lcu)
        idfa = radf1 * exp (-kdf * lcu)
        idra = radf2 * exp (radf4 * lcu)
        idrdra = radf3 * exp (radf4 * lcu)
        isha = idfa + (idra - idrdra)
        isla = isha + radf3
      endif
      if (c3) then
        call c3phot (co2,isha,nl(layer),n3(layer),fabsb(layer),ash)
      else
        call c4phot (co2,isha,tp(0,2),nl(layer),n3(layer), &
      &             fabsb(layer),ash)
      endif
      if (radf2 > 0.0) then
        if (c3) then
          call c3phot (co2,isla,nl(layer),n3(layer),fabsb(layer),asl)
        else
          call c4phot (co2,isla,tp(0,2),nl(layer),n3(layer), &
          &           fabsb(layer),asl)
        endif
      endif
      acangl = fsl * asl + (1.0 - fsl) * ash
      acang = acang + acangl
      lcu = lcu + laid
    enddo
    !------------------------------------------------------------------!
  else
    acang = 0.0
  endif
  acang = laid * acang
  acann = acang - rpcrn
  if (ice == 1) acann_save = acann
  if (ice == 2) ansat      = acann
enddo
acann = acann_save
!----------------------------------------------------------------------!
 f_ci = (co2i + 9.6247) / (co2i + eps_sing)
 f_ci = max (0.0, f_ci)
 f_ci = min (1.0, f_ci)
 cnc = f_ci * f_dq * f_swp * 650e-6 * ansat
 cnc = max (cnc,alai*0.0001d0) ! need to improve
 cnc = min (cnc,alai*0.005d0) ! need to improve
 rtot = 1.0 / (cnc + eps_sing) + 1.0 / (ch * vsm + eps_sing)
 acanni = (capa - co2i) * 1.0e6 / (1.6 * rt * rtot + eps_sing)
!----------------------------------------------------------------------!
afunc = acann - acanni
!----------------------------------------------------------------------!
end function afunc
!======================================================================!
subroutine c3phot (co2,light,nl_l,n3_l,fabsb_l,ag)
!----------------------------------------------------------------------!
use veg_parameters
use veg_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
real*8 :: co2,light,n3_l,ag,m1,m2,msat,agsf1,agsf2,nco,nlim,nsat,nl_l
real*8 :: fabs,fabsb_l
!----------------------------------------------------------------------!
m1 = (co2 + eps_sing) / (co2 + 2.0 * pcp + eps_sing)
m2 = (co2 + eps_sing) / (co2 + km        + eps_sing)
msat = min (m1 * n1, m2 * n2)
agsf1 = 1.0 - pcp / (co2 + eps_sing)
agsf2 = agsf1 * alpha_c3 * m1
agsf1 = agsf1 * msat
nco = msat / (alpha_c3 * ka * m1)
nlim = -dlog (eps_sing + nco / (light * n3_l + eps_sing)) / &
&      (ka * n3_l)
if (nlim < 1.0) then
  nsat = 0.0
elseif (nlim > nl_l) then
  nsat = nl_l
else
  nsat = nlim
endif
fabs = exp (-ka * n3_l * nsat) - fabsb_l
ag = agsf1 * nsat + agsf2 * light * fabs
!----------------------------------------------------------------------!
end subroutine c3phot
!======================================================================!
subroutine c4phot (co2,light,tcan,nl_l,n3_l,fabsb_l,ag)
!----------------------------------------------------------------------!
use veg_parameters
use veg_variables
use atmos_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
real*8 :: co2,light,tcan,nl_l,n3_l,ag,m1,m2,m3,msat,n3i,kti,agsf1,agsf2
real*8 :: nco,nlim,nsat,fabs,fabsb_l
real*8, parameter :: alpha_c4 = 0.05
real*8, parameter :: alphai = 0.08 ! Intrinsic quantum yield
real*8, parameter :: alphar = 0.11 ! RuBP quantum requirement (mol/mol)
real*8, parameter :: f      = 0.6 ! Fractional RuBP quantum requirement.
!----------------------------------------------------------------------!
n3i = 0.7 / 220.0
kti = 2 ** ((tcan - 25.0) / 10.0)
m1 = 7800.0 / (7800.0 + 2.0 * pcp + eps_sing)
m2 = 7800.0 / (7800.0 + km        + eps_sing)
m3 = nfp * 1.0e4 * (kti * co2 / ps) / 1.6
msat = min (m1 * n1, m2 * n2, m3 * n3i)
agsf1 = (1.0 - pcp / 7800.0)
agsf2 = agsf1 * f * alphai * m1
agsf1 = agsf1 * msat
nco = msat / (f * alphai * ka * m1)
nlim = -dlog (nco / (light * n3_l + eps_sing)) / &
&      (ka * n3_l)
if (nlim < 0.0) then
  nsat = 0.0
elseif (nlim > nl_l) then
  nsat = nl_l
else
  nsat = nlim
endif
fabs = exp (-ka * n3_l * nsat) - fabsb_l
ag = agsf1 * nsat + agsf2 * light * fabs
!----------------------------------------------------------------------!
end subroutine c4phot
!======================================================================!
