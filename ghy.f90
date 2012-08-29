!======================================================================!
subroutine reth
!----------------------------------------------------------------------!
! Calculate relative volumetric soil moisture (theta) from water depths
! in soil layers; snow depth from snow in each snow layer; fraction of
! wet canopy; fraction of snow that is exposed (not masked).
!----------------------------------------------------------------------!
use control_parameters
use constants
use soil_variables
use veg_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  do  nsoil_layer = 1, nsoil_layers
    !------------------------------------------------------------------!
    ! Absolute volumetric soil moisture (m^3/m^3).
    !------------------------------------------------------------------!
    theta (nsoil_layer,ibv) = w (nsoil_layer,ibv) / dz (nsoil_layer)
    !------------------------------------------------------------------!
    ! Fractional saturation (fraction).
    !------------------------------------------------------------------!
    if (ws (nsoil_layer,ibv) > 0.0d0) then
      soilsat (nsoil_layer,ibv) = w  (nsoil_layer,ibv) / &
      &                           ws (nsoil_layer,ibv)
    else
      soilsat (nsoil_layer,ibv) = 0.0d0
    endif
    !------------------------------------------------------------------!
  end do
end do
!----------------------------------------------------------------------!
! Fraction of canopy covered by water (fraction).
!----------------------------------------------------------------------!
if (process_vege .and. ws (0,2) > 0.0d0 ) then
  theta (0,2) = (w (0,2) / ws (0,2)) ** (2.0d0 / 3.0d0)
else
  theta (0,2) = 0.0d0
endif
theta (0,2) = min (theta (0,2), one)
!----------------------------------------------------------------------!
! Set up snow depth variables.
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  snowd (ibv) = 0.0d0
  do lsn = 1, nsnow_layers (ibv)
    ! Compute snowd as if all snow was distributed uniformly
    ! over the cell (i.e. snowd = wsn * sn_frac / 1.d0 )
    snowd (ibv) = snowd (ibv) + wsn (lsn,ibv) * fr_snow (ibv)
  enddo
enddo
!----------------------------------------------------------------------!
! Fraction of wet and dry canopy (fraction).
!----------------------------------------------------------------------!
fw = theta (0,2)
!----------------------------------------------------------------------!
! Fraction of snow that is exposed (fraction).
!----------------------------------------------------------------------!
fm = 1.0d0 - exp (-snowd (2) / (snowm + 1.0d-12))
! Correct fraction of wet canopy by snow fraction.
fd = 1.0d0 - fw
!----------------------------------------------------------------------!
! Actual fraction of canopy that is dry (fraction).
!----------------------------------------------------------------------!
fd0 = fd
!----------------------------------------------------------------------!
end subroutine reth
!======================================================================!
subroutine hydra
!----------------------------------------------------------------------!
! Computes water potentials of layers (matric and gravitational).
!----------------------------------------------------------------------!
! Input:
! theta : relative volumetric water in layers (m)
!----------------------------------------------------------------------!
! Output
! h : water potential in layers (m);
!----------------------------------------------------------------------!
use constants
use soil_variables
use veg_variables
use control_parameters
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
real*8 :: hl,thr,thr0,thr1,thr2
real*8 :: temp,d1,d2,dl,xku1,xku2
real*8 :: xklu,xk1,xkl,xkud
integer :: i,j,j1,j2,jc
integer :: ith
!----------------------------------------------------------------------!
xkud = 2.78d-5
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  xk (nsoil_layers+1,ibv) = 0.0d0
  xku (0,ibv) = 0.0d0
  do nsoil_layer = 1, nsoil_layers
    !------------------------------------------------------------------!
    ! Calculate water potential (m).
    !------------------------------------------------------------------!
    j1 = 0
    j2 = nth
    thr1 = thets (nsoil_layer,ibv)
    thr2 = thetm (nsoil_layer,ibv)
    thr0 = theta (nsoil_layer,ibv)
    thr0 = min (thr1,thr0)
    thr0 = max (thr2,thr0)
    do jc = 1, jcm
      j = (j1 + j2) / 2
      thr = 0.0d0
      do i = 1, ntextures - 1
        thr = thr + thm (j,i) * q (i,nsoil_layer)
      enddo
      if ((thr - thr0) < 0.0d0) then
!     here thr is too small, bisect on low j end
        j2   = j
        thr2 = thr
      elseif ((thr - thr0) > 0.0d0) then
!     here thr is too large, bisect on high j end
        j1   = j
        thr1 = thr
      else                ! i.e. .eq.
!     here thr is equal to thr0
        hl   = hlm (j)
        j1   = j
        thr1 = thr0
!     the strange value for thr2 below is only for calculating temp
        thr2 = -10.0d0
        goto 500
      endif
    enddo                ! jc
!     here theta is between two adjacent thr''s. interpolate.
    hl = (hlm (j1) * (thr0 - thr2) + hlm (j2) * (thr1 - thr0)) &
    &    / (thr1 - thr2)
 500      continue
    h (nsoil_layer,ibv) = hl
    !------------------------------------------------------------------!
    ! Calculate diffusivity (kg/m/s).
    !------------------------------------------------------------------!
    ith = j1
    temp = (thr1 - thr0) / (thr1 - thr2)
    d1 = 0.0d0
    d2 = 0.0d0
    xku1 = 0.0d0
    xku2 = 0.0d0
    xkus (nsoil_layer,ibv) = 0.0d0
    do i = 1, ntextures - 1
      d1 = d1 + q (i,nsoil_layer) * dlm (ith  ,i)
      d2 = d2 + q (i,nsoil_layer) * dlm (ith+1,i)
      xku1 = xku1 + q (i,nsoil_layer) * xklm (ith  ,i)
      xku2 = xku2 + q (i,nsoil_layer) * xklm (ith+1,i)
      xkus(nsoil_layer,ibv) = xkus (nsoil_layer,ibv) + &
      &                       q (i,nsoil_layer) * xklm (0,i)
    enddo ! i
    dl = (1.0d0 - temp) * d1 + temp * d2
    dl = (1.0d0 - fice (nsoil_layer,ibv)) * dl
    d (nsoil_layer,ibv) = dl
    !------------------------------------------------------------------!
    ! Calculate conductivity (m/s).
    !------------------------------------------------------------------!
    xklu = (1.0d0 - temp) * xku1 + temp * xku2
    xklu = (1.0d0 - fice (nsoil_layer,ibv)) * xklu
    xku (nsoil_layer,ibv) = xklu
    if (nsoil_layer == 1) then
      xk1 = 0.0d0
      do i = 1, ntextures - 1
        xk1 = xk1 + qk (i,1) * xklm (0,i)
      enddo
      xkl = xk1
      xkl = xkl / (1.0d0 + xkl / (-zc (1) * xkud))
      xkl = (1.0d0 - fice (1,ibv) * theta (1,ibv) / thets (1,ibv)) * xkl
      xkl = max (zero, xkl)
      xk (1,ibv) = xkl
    else
      xk (nsoil_layer,ibv) = sqrt (xku (nsoil_layer - 1,ibv) &
      & * xku (nsoil_layer,ibv))
    endif
    !------------------------------------------------------------------!
  enddo ! nsoil_layer
enddo ! ibv
!----------------------------------------------------------------------!
! Add gravitational potential to h.
!----------------------------------------------------------------------!
do nsoil_layer = 1, nsoil_layers
  do ibv = i_bare, i_vege
    h (nsoil_layer,ibv) = h (nsoil_layer,ibv) + zc (nsoil_layer)
  enddo
enddo
!----------------------------------------------------------------------!
end subroutine hydra
!======================================================================!
subroutine hl0
!----------------------------------------------------------------------!
! Set up table of volumetric moisture contents at specific matric
! potentials. h, matric potential, is tabulated in a geometric series
! from 0 to hmin, with a first step of delh1. Theta (moisture content
! relative to saturation value for that texture)
! depends not only on the matric potential, but also on the soil
! texture. Cubic equation solved to determine theta as a function of h.
! Also set up conductivity and diffusivity tables.
!----------------------------------------------------------------------!
use constants
use control_parameters
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
integer :: i ! Texture loop count                                    (n)
integer :: j ! Table loop count                                      (n)
integer :: m ! Loop count for solution to geometric function         (n)
integer :: k ! Loop count for conductivity coefficients              (n)
real*8 :: testh ! Solution to geometric function                 (m3/m3)
real*8 :: func,dfunc  ! Results of geometric functions (m?)
real*8 :: diff ! Ratio of functions
real*8 :: hs,a1,a2,a3 ! Matric potential function parameters (?)
real*8 :: delh1 ! First step for matric potential increment          (m)
real*8 :: delhn ! Matric potential increment                         (m)
real*8 :: hmin  ! Minimum matric potential                           (m)
real*8 :: s,alph0,alph0o,alpls1 ! Step calculators
real*8 :: arg(ntextures-1) ! Arguments to calculate conductivities   (?)
                           ! and diffusivities                       (?)
!----------------------------------------------------------------------!
! Number used to index maximum element of thm (n).
!----------------------------------------------------------------------!
nth = 2 ** nexp
!----------------------------------------------------------------------!
! Set up array of matric potentials, hlm                             (m)
!----------------------------------------------------------------------!
hlm (0) = 0.0d0 ! First matric potential (m).
delh1 = -0.00625d0 ! First step.
hmin = -1000.d0 ! Minimum matric potential (m).
delhn = delh1 ! Initialise step.
s = hmin / delh1
alph0 = 1.0d0 / 8.0d0
 10   alph0o = alph0
alph0 = (s * alph0 + 1.0d0) ** (1.0d0 / nth) - 1.0d0
if (abs (alph0o - alph0) >= 1.0d-8) go to 10
alpls1 = 1.0d0 + alph0 ! Step factor.
do j = 1, nth ! Loop over matric potentials.
  hlm (j) = hlm (j-1) + delhn ! Set matric potential (m).
  delhn = alpls1 * delhn ! New hlm increment (m).
enddo
!----------------------------------------------------------------------!
! Set up table of thetas at various matric potentials for each texture.
!----------------------------------------------------------------------!
do i = 1, ntextures - 1
  thm (0,i) = 1.00d0 ! Initial estimate for table values.
  do j = 1, nth
    ! Calculate geometric function parameters.
    hs= -exp (c * (a (1,i) + a (2,i) + a (3,i) + a (4,i)))
    a1 = a (3,i) / a (4,i)
    a2 = (a (2,i) - (log (-hlm (j) - hs)) / c) / a (4,i)
    a3 = a (1,i) / a (4,i)
    testh = thm (j-1,i)
    do m = 1, mmax
      func = (testh ** 3) + (a1 * (testh ** 2)) + (a2 * testh) + a3
      dfunc= (3 * testh ** 2) + (2 * a1 * testh) + a2
      diff = func / dfunc
      testh = testh - diff
      if (abs (diff) < xtol) goto 140
    enddo
    write (98,*) 'max # iterations for thm',mmax
140 thm (j,i) = testh
  enddo
enddo
!----------------------------------------------------------------------!
! Set up table of soil moisture conductivities.
!----------------------------------------------------------------------!
do j = 0, nth
  do i = 1, ntextures - 1
    xklm (j,i) = 0.0d0 ! Initial hydraulic conductivity (kg/m^2/s)
    arg (i) = 0.0d0
    do k = -1, 2
      arg (i) = arg (i) + b (k+2,i) * thm (j,i) ** k
    enddo
    arg (i) = min (arg (i), sxtn)
    arg (i) = max (arg (i), -sxtn)
    xklm (j,i) = exp (c * arg(i))
  enddo
  !--------------------------------------------------------------------!
  ! Set up table of soil moisture diffusivities.
  !--------------------------------------------------------------------!
  do i = 1, ntextures - 1
    dlm (j,i) = 0.0d0 ! Initial hydraulic diffusivity (kg/m/s).
    arg (i) = 0.0d0
    do k = -1, 2
      arg (i) = arg (i) + pq (k+2,i) * thm (j,i) ** k
    enddo
    arg (i) = min (arg(i),  sxtn)
    arg (i) = max (arg(i), -sxtn)
    dlm (j,i) = exp (c * arg (i))
  enddo
  !--------------------------------------------------------------------!
enddo
!----------------------------------------------------------------------!
! Convert volumetric moisture contents to relative values (fraction).
!----------------------------------------------------------------------!
do j = 0, nth
  do i = 1, ntextures - 1
    thm (j,i) = thm (j,i) * sat (i)
  enddo
enddo
!----------------------------------------------------------------------!
end subroutine hl0
!======================================================================!
subroutine evap_limits (compute_evap,evap_max_out,fr_sat_out)
!----------------------------------------------------------------------!
! Compute transpiration from each soil layer (evapdl).
!----------------------------------------------------------------------!
use soil_variables
use veg_variables
use constants
use atmos_variables
use control_parameters
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
logical, intent(in) :: compute_evap
real*8, intent(out) :: evap_max_out,fr_sat_out
!----------------------------------------------------------------------!
! Wilting point (m).
!----------------------------------------------------------------------!
real*8, parameter :: hw = -100.0d0
!----------------------------------------------------------------------!
! Potential evaporation for vegetated soil (m/s).
!----------------------------------------------------------------------!
real*8 :: epvg
!----------------------------------------------------------------------!
! Specific humidity of bare soil (kg/kg).
!----------------------------------------------------------------------!
real*8 :: qb
real*8 :: qbs
!----------------------------------------------------------------------!
! Specific humidity of vegetated ground? (kg/kg).
!----------------------------------------------------------------------!
real*8 :: qvg
real*8 :: qvs
!----------------------------------------------------------------------!
! Canopy evap limited by soil surface water, precipitation, and total
! soil water (m/s).
!----------------------------------------------------------------------!
real*8 :: evap_max_vegsoil
real*8 :: evap_max_nsat
real*8 :: evap_max_snow(2)
real*8 :: evap_max_sat
real*8 :: fr_sat
!----------------------------------------------------------------------!
! Condensation limit due to atmosphere water (m/s).
!----------------------------------------------------------------------!
real*8 :: qm1dt
!----------------------------------------------------------------------!
! Surface atmosphere level saturated vapour mixing ratio (kg/kg).
!----------------------------------------------------------------------!
real*8 :: qsat_modele
!----------------------------------------------------------------------!
! Dry canopy evapoporation limited by stomatal conductance (m/s).
!----------------------------------------------------------------------!
real*8 :: pot_evap_can
!----------------------------------------------------------------------!
! Total canopy beta (fraction).
!----------------------------------------------------------------------!
real*8 :: betat
!----------------------------------------------------------------------!
! Soil water potential in layer (MPa).
!----------------------------------------------------------------------!
real*8 :: swpl
!----------------------------------------------------------------------!
! Place holder.
!----------------------------------------------------------------------!
tsn1 (:) = tp (1,1)
!----------------------------------------------------------------------!
! cna is the conductance of the atmosphere (m/s).
!----------------------------------------------------------------------!
cna = ch * vsm
!----------------------------------------------------------------------!
! Divide by rho_water to get flux in m/s.
!----------------------------------------------------------------------!
rho3 = rho / rhow
!----------------------------------------------------------------------!
! Make sure important variables are initialised (need for ibv hack).
!----------------------------------------------------------------------!
evap_max      (:) = 0.0d0
evap_max_snow (:) = 0.0d0
evap_max_wet  (:) = 0.0d0
evap_max_dry  (:) = 0.0d0
betadl (:) = 0.0d0 ! Transp. efficiency of each soil layer    (fraction)
betad      = 0.0d0 ! Overall transpiration efficiency         (fraction)
f_swpl (:) = 0.0d0 ! Phys. impact of each soil layer          (fraction)
f_swp      = 0.0d0 ! Overall soil moist. impact on physiology (fraction)
!----------------------------------------------------------------------!
! Soil moisture evaporation limitation.
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  evap_max (ibv) = 0.0d0
  do nsoil_layer = 1, nsoil_layers
    evap_max (ibv) = evap_max (ibv) + &
    &                (w (nsoil_layer,ibv) - dz (nsoil_layer) * &
    &                thetm (nsoil_layer,ibv)) / dt
  enddo
enddo
! Maximal evaporation from the snow fraction may be too restrictive
! with resp. to pr. but will leave for now.
do ibv = i_bare, i_vege
  evap_max_snow (ibv) = pr
  do lsn = 1, nsnow_layers (ibv)
    evap_max_snow (ibv) = evap_max_snow (ibv) + &
    &                     wsn (lsn,ibv) / dt
  enddo
enddo
!----------------------------------------------------------------------!
! Evaporation from bare soil>
!----------------------------------------------------------------------!
if (process_bare) then
  !--------------------------------------------------------------------!
  ibv = 1
  !--------------------------------------------------------------------!
  ! "No support for saturated soil yet, setting just in case...'
  !--------------------------------------------------------------------!
  evap_max_wet (ibv) = evap_max (ibv) + pr
  !--------------------------------------------------------------------!
  ! Evap limited by diffusion and precipitation (m/s).
  !--------------------------------------------------------------------!
  evap_max_dry (ibv) = min (evap_max (ibv), &
  &  2.467d0 * d (1,1) * (theta (1,1) - thetm (1,1)) / dz (1) + pr)
  !--------------------------------------------------------------------!
endif
!----------------------------------------------------------------------!
! Evaporation from the canopy.
!----------------------------------------------------------------------!
if (process_vege) then
  !--------------------------------------------------------------------!
  ibv = 2
  !--------------------------------------------------------------------!
  ! "pr does not work for snow"
  !--------------------------------------------------------------------!
  evap_max_wet (ibv) = w (0,2) / dt
  !--------------------------------------------------------------------!
  ! Evaporation limit for soil under the canopy (m/s).
  !--------------------------------------------------------------------!
  evap_max_vegsoil = min (evap_max (ibv), &
  &  2.467d0 * d (1,2) * (theta (1,2) - thetm (1,2)) / dz (1) + pr)
  !--------------------------------------------------------------------!
  ! Root f_swp for physiology (fraction).
  !--------------------------------------------------------------------!
  betad = 0.0d0
  f_swp = 0.0d0
  do nsoil_layer = 1, nsoil_layers
    betadl (nsoil_layer) = (1.0d0 - fice (nsoil_layer,2)) * &
    &                      fr (nsoil_layer) * max ((hw -    &
    &                      h (nsoil_layer,2)) / hw, zero)
    swpl = rhow * grav * h (nsoil_layer,2) / 1.0d6
    f_swpl (nsoil_layer) = 1.0d0 - swpl / (-1.5d0)
    f_swpl (nsoil_layer) = min (1.0d0,f_swpl(nsoil_layer))
    f_swpl (nsoil_layer) = max (0.0d0,f_swpl(nsoil_layer))
    f_swpl (nsoil_layer) = fr (nsoil_layer) * f_swpl(nsoil_layer)
    betad = betad + betadl (nsoil_layer)
    f_swp  = f_swp  + f_swpl  (nsoil_layer)
  enddo
  betadl (:) = f_swpl (:)
  betad = f_swp
  if (betad < 1.0d-12) betad = 0.0d0 ! to avoid 0/0 divisions
  if (f_swp  < 1.0d-12) f_swp  = 0.0d0 ! to avoid 0/0 divisions
  !--------------------------------------------------------------------!
  ! Canopy water vapour mixing ratio (kg/kg).
  !--------------------------------------------------------------------!
  qv = qsat_modele (tp(0,2)+tf,lhe,ps)
  !--------------------------------------------------------------------!
  ! Get canopy conductance.
  !--------------------------------------------------------------------!
  call veg_conductance
  !--------------------------------------------------------------------!
  ! Total beta due to canopy conductance (fraction).
  !--------------------------------------------------------------------!
  betat = cnc / (cnc + cna + 1.0d-12)
  !--------------------------------------------------------------------!
  ! Dry canopy evaporation limited by stomatal conductance (m/s).
  !--------------------------------------------------------------------!
  pot_evap_can = betat * rho3 * cna * (qv - qs)
  !--------------------------------------------------------------------!
  evap_max_dry (ibv) = 0.0d0
  if (betad > 0.0d0) then
    do nsoil_layer = 1, nsoil_layers
      evap_max_dry (ibv) = evap_max_dry (ibv) + min (pot_evap_can * &
      &                    betadl (nsoil_layer) / betad, &
      &                    (w (nsoil_layer,ibv) - dz (nsoil_layer) * &
      &                    thetm (nsoil_layer,ibv)) / dt )
    enddo
  endif
  !--------------------------------------------------------------------!
endif
! Now we have to add the fluxes according to fractions.
! Bare soil.
evap_max_sat = fb * fr_snow (1) * evap_max_snow (1)
evap_max_nsat = fb * (1.0d0 - fr_snow (1)) * evap_max_dry (1)
! Canopy
evap_max_sat = evap_max_sat + fv * (fr_snow (2) * fm * &
&              evap_max_snow (2) + (1.0d0 - fr_snow (2) * fm) * &
&              theta (0,2) * evap_max_wet (2))
evap_max_nsat = evap_max_nsat + fv * (( 1.0d0 - fr_snow (2) * fm) &
&               * (1.0d0 - theta (0,2)) * evap_max_dry (2))
fr_sat = fb * fr_snow (1) + fv * (fr_snow (2) * fm + &
&        (1.0d0 - fr_snow (2) * fm) * theta (0,2))
! Set variables for output.
evap_max_out = evap_max_nsat
fr_sat_out = fr_sat
! Have to pass these flux limits.
if (.not. (evap_max_out > 0.0d0 .or. evap_max_out <= 0.0d0)) then
  write (99,*) 'evap_max_out bounds',i0,j0
  write (99,*) 'cnc cna',cnc,cna
  write (99,*) 'evap_max_out,evap_max(1),evap_max(2)',&
    &           evap_max_out,evap_max(1),evap_max(2)
  write (99,*) 'evap_max_dry(1),evap_max_dry(2)', &
  &             evap_max_dry(1),evap_max_dry(2)
  write (99,*) 'pot_evap_can',pot_evap_can
  write (99,*) 'betat',betat
  write (99,*) 'betad',betad
  write (99,*) 'betadl',betadl(:)
  write (99,*) 'w(1)',w(:,1)
  write (99,*) 'w(2)',w(:,2)
  write (99,*) 'thetm(1)',thetm(:,1)
  write (99,*) 'thetm(2)',thetm(:,2)
  write (99,*) 'evapvw evapvd',evapvw,evapvd
  stop
endif
!----------------------------------------------------------------------!
! qm1 has mass of water vapor in first atmosphere layer (g/m^2).
qm1dt = 0.001d0 * qm1 / dt
!----------------------------------------------------------------------!
! Bare soil water vapour mixing ratio (kg/kg).
!----------------------------------------------------------------------!
qb = qsat_modele (tp(1,1)+tf,lhe,ps)
qbs = qsat_modele (tsn1(1)+tf,lhe,ps)
!----------------------------------------------------------------------!
! Canopy water vapour mixing ratio (kg/kg).
!----------------------------------------------------------------------!
qv = qsat_modele (tp(0,2)+tf,lhe,ps)
qvs = qsat_modele (tsn1(2)+tf,lhe,ps)
!----------------------------------------------------------------------!
! Vegetated ground? water vapour mixing ratio (kg/kg).
!----------------------------------------------------------------------!
qvg = qsat_modele (tp(1,2)+tf,lhe,ps)
!----------------------------------------------------------------------!
! Potential evaporation for bare soil (m/s).
!----------------------------------------------------------------------!
epb = rho3 * cna * (qb - qs)
epbs = rho3 * cna * (qbs - qs)
!----------------------------------------------------------------------!
! Potential evaporation for canopy (m/s).
!----------------------------------------------------------------------!
epv = rho3 * cna * (qv - qs)
epvs = rho3 * cna * (qvs - qs)
!----------------------------------------------------------------------!
! Potential evaporation for vegetated soil (m/s).
! Stated in modelE as not actually correct (possibly because qvg is
! not relevant qv for total evaporation here). Later, make this correct
! and see what effect is.
!----------------------------------------------------------------------!
epvg = rho3 * cna * (qvg - qs)
!----------------------------------------------------------------------!
! Bare soil evaporation.
!----------------------------------------------------------------------!
if (process_bare) then
  !--------------------------------------------------------------------!
  evapb = min (epb, evap_max_dry (1))
  !--------------------------------------------------------------------!
  ! Condensation limit due to atmosphere water (m/s).
  !--------------------------------------------------------------------!
  evapb = max (evapb, -qm1dt)
  evapbs = min (epbs, evap_max_snow (1))
  evapbs = max (evapbs, -qm1dt)
else
  evapb = 0.0d0
  evapbs = 0.0d0
  !--------------------------------------------------------------------!
endif
!----------------------------------------------------------------------!
! Vegetated soil evaporation.
!----------------------------------------------------------------------!
if (process_vege) then
  ! In case.
  evapvg = 0.0d0
  !--------------------------------------------------------------------!
  ! Evaporation from wet canopy (m/s).
  !--------------------------------------------------------------------!
  evapvw = min (epv, evap_max_wet (2))
  evapvw = max (evapvw,-qm1dt)
  if (.not. (evapvw > 0.0d0 .or. evapvw <= 0.0d0)) then
    write (99,*) 'evapvw',evapvw
    write (99,*) 'epv evap_max_wet(2)',epv,evap_max_wet(2)
    write (99,*) '-qm1dt',-qm1dt
  endif
  !--------------------------------------------------------------------!
  ! Dry evaporation (transpiration) from canopy (m/s).
  !--------------------------------------------------------------------!
  evapvd = min (epv, evap_max_dry (2))
  evapvd = max (evapvd, 0.0d0) ! No dew allowed?
  evapvs = min (epvs, evap_max_snow (2))
  evapvs = max (evapvs, -qm1dt)
  !--------------------------------------------------------------------!
  ! Evaporation from soil under canopy (m/s).
  !--------------------------------------------------------------------!
  evapvg = min (epvg, evap_max_vegsoil)
  evapvg = min (evapvg, epv-evapvd*fd-evapvw*fw )
  evapvg = max (evapvg, 0.0d0)
  !--------------------------------------------------------------------!
  ! See if have dew, if so let fall on entire canopy.
  !--------------------------------------------------------------------!
  if (evapvw < 0.0d0) then
    fw = 1.0d0
    fd = 0.0d0
  endif
  !--------------------------------------------------------------------!
else
  evapvw = 0.0d0
  evapvd = 0.0d0
  evapvs = 0.0d0
  evapvg = 0.0d0
endif
!----------------------------------------------------------------------!
! Compute transpiration for separate layers (=0 for bare soil).
!----------------------------------------------------------------------!
evapdl (1:nsoil_layers,1:2) = 0.0d0
if (betad > 0.0d0) then
  do nsoil_layer = 1, nsoil_layers
    evapdl (nsoil_layer,2) = evapvd * betadl (nsoil_layer) / betad
  enddo
endif
!----------------------------------------------------------------------!
end subroutine evap_limits
!======================================================================!
subroutine sensible_heat
!----------------------------------------------------------------------!
! Compute sensible heat of bare/vegetated soil and snow.
!----------------------------------------------------------------------!
use soil_variables
use veg_variables
use atmos_variables
use constants
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
cna = ch * vsm
!----------------------------------------------------------------------!
! Sensible heat flux from bare soil to air, no snow (J/s).
!----------------------------------------------------------------------!
snsh (1) = sha * rho * cna * (tp (1,1) - ts + tf)
!----------------------------------------------------------------------!
! Canopy (J/s).
!----------------------------------------------------------------------!
snsh (2) = sha * rho * cna * (tp (0,2) - ts + tf)
!----------------------------------------------------------------------!
! Bare soil snow (J/s).
!----------------------------------------------------------------------!
snshs (1) = sha * rho * cna * (tsn1 (1) - ts + tf)
!----------------------------------------------------------------------!
! Canopy snow (J/s).
!----------------------------------------------------------------------!
snshs (2) = sha * rho * cna * (tsn1 (2) - ts + tf)
!----------------------------------------------------------------------!
! Derivative is the same for all above. used?
!----------------------------------------------------------------------!
!dsnsh_dt = sha * rho * cna (2)
!----------------------------------------------------------------------!
end subroutine sensible_heat
!======================================================================!
subroutine drip_from_canopy
!----------------------------------------------------------------------!
! Computes the flux of drip water (and its heat cont.) from canopy
! The drip water is split into liquid water dripw() and snow drips()
! For bare soil fraction the canopy is transparent
!----------------------------------------------------------------------!
use control_parameters
use constants
use atmos_variables
use soil_variables
use veg_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
! Place holders.
!----------------------------------------------------------------------!
real*8, parameter :: snow_cover_same_as_rad = 0
!----------------------------------------------------------------------!
real*8 ptmp,ptmps,pfac,pm,pmax
real*8 snowf,dr,drs
real*8 melt_w, melt_h
!----------------------------------------------------------------------!
! Calculate snow fall.  snowf is snow fall, m s-1 of water depth.
!----------------------------------------------------------------------!
snowf = 0.0d0
!***fix for now if (htpr < 0.0d0) snowf = min(-htpr/fsn, pr)
!----------------------------------------------------------------------!
if (process_vege) then
  ptmps = -evapvw * fw
  ptmp = pr - snowf
  !---------------------------------------------------------------------!
  ! Use effects of subgrid-scale precipitation to calculate drip.
  !---------------------------------------------------------------------!
  pm = 1.0d-6
  pmax = fd0 * pm
  drs = max (ptmps-pmax, zero)
  dr = drs
  if (ptmp > 0.0d0) then
    pfac = (pmax - ptmps) * prfr / ptmp
    if (pfac >= 0.0d0) then
      if (pfac < 30.0d0) dr = ptmp * exp (-pfac)
    else
      dr = ptmp + ptmps - pmax
    endif
  endif
  !---------------------------------------------------------------------!
  ! Make sure "dr" makes sense.
  !---------------------------------------------------------------------!
  dr = min (dr, pr-snowf-evapvw*fw )
  dr = max (dr, pr-snowf-evapvw*fw - (ws(0,2)-w(0,2))/dts)
  dr = max (dr, 0.0d0 )    ! just in case (probably don''t need it)
  dripw   (2) = dr
  ! Don't allow it to freeze.
  htdripw (2) = shw * dr * max (tp(0,2),0.0d0) 
  !---------------------------------------------------------------------!
  ! Snow falls through the canopy.
  !---------------------------------------------------------------------!
  drips  (2) = snowf
  htdrips(2) = min (htpr, 0.d0 ) ! liquid H20 is 0 C, so no heat
else  ! no vegetated fraction.
  dripw   (2) = 0.0d0
  htdripw (2) = 0.0d0
  drips   (2) = 0.0d0
  htdrips (2) = 0.0d0
endif
!-----------------------------------------------------------------------!
! For bare soil drip is precipitation.
!-----------------------------------------------------------------------!
drips   (1) = snowf
htdrips (1) = min (htpr, 0.0d0 )
dripw   (1) = pr   - drips   (1)
htdripw (1) = htpr - htdrips (1)
if (snow_cover_same_as_rad .ne. 0 ) then 
  do ibv=1, 2
    if (tp (1,ibv) > 0.0d0 ) then
      melt_w = (1.0d0 - fr_snow (ibv)) * drips   (ibv)
      melt_h = (1.0d0 - fr_snow (ibv)) * htdrips (ibv)
      drips   (ibv) = drips   (ibv) - melt_w
      htdrips (ibv) = htdrips (ibv) - melt_h
      dripw   (ibv) = dripw   (ibv) + melt_w
      htdripw (ibv) = htdripw (ibv) + melt_h
    endif
  enddo
endif
!----------------------------------------------------------------------!
end subroutine drip_from_canopy
!======================================================================!
subroutine flg
!----------------------------------------------------------------------!
! Computes the water fluxes at the surface.
!----------------------------------------------------------------------!
use veg_variables
use soil_variables
use atmos_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
! Place holders.
!----------------------------------------------------------------------!
real*8, parameter :: flmlt (nsurf_max) = 0.0d0
real*8, parameter :: flmlt_scale (nsurf_max) = 0.0d0
!----------------------------------------------------------------------!
! Bare soil.
!----------------------------------------------------------------------!
if (process_bare) then
  f (1,1) = -flmlt (1) * fr_snow (1) - flmlt_scale (1) &
  &         - (dripw (1) - evapb) * (1.0d0 - fr_snow (1))
endif
!----------------------------------------------------------------------!
! Vegetated soil.
!----------------------------------------------------------------------!
if (process_vege) then
  !--------------------------------------------------------------------!
  ! Upward flux from wet canopy.
  !--------------------------------------------------------------------!
  fc (0) = -pr + evapvw * fw * (1.0d0 - fm * fr_snow (2))
  !--------------------------------------------------------------------!
  ! Snow masking of pr is ignored since it is not included into drip.
  !--------------------------------------------------------------------!
  fc (1) = -dripw (2) - drips (2)
  !--------------------------------------------------------------------!
  ! f(1,2) is a flux up from the soil.
  !--------------------------------------------------------------------!
  f (1,2) = -flmlt (2) * fr_snow (2) - flmlt_scale (2) &
  &         - (dripw (2) - evapvg) * (1.0d0 - fr_snow (2))
  !--------------------------------------------------------------------!
endif
!----------------------------------------------------------------------!
! Compute evap_tot for accumulators.
!----------------------------------------------------------------------!
evap_tot (1) = evapb * (1.0d0 - fr_snow (1)) + evapbs * fr_snow (1)
evap_tot (2) = (evapvw * fw + evapvd * fd) * (1.0d0 - fr_snow (2) * fm)&
&              + evapvs * fr_snow (2) * fm + evapvg * (1.0d0 - &
&              fr_snow (2))
!----------------------------------------------------------------------!
end subroutine flg
!======================================================================!
subroutine flhg
!----------------------------------------------------------------------!
! Calculate the ground heat fluxes (to the surface).
!----------------------------------------------------------------------!
use constants
use soil_variables
use control_parameters
use veg_variables
use atmos_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
real*8 :: thrm_can,thrm_soil(2)
!----------------------------------------------------------------------!
! Place holders.
!----------------------------------------------------------------------!
fhsng (:) = 0.0d0
fhsng_scale (:) = 0.0d0
if (htpr < 0.0d0) fhsng_scale (:) = htpr ! fix for now
thrmsn (:) = 0.0d0
!----------------------------------------------------------------------!
thrm_can = stbo * (tp (0,2) + tf) ** 4
thrm_soil (1:2) = stbo * (tp(1,1:2) + tf) ** 4
!thrm_can = 0.0d0 ! holder
!thrm_soil (:) = 0.0d0 ! holder
!----------------------------------------------------------------------!
if (process_bare) then
  fh (1,1) = -fhsng (1) * fr_snow (1) - fhsng_scale (1) &
  &          + (-htdripw (1) + &
  &          evapb * elh + snsh (1) + thrm_soil (1) - srht - trht) &
  &          * (1.0d0 - fr_snow (1))
endif
!----------------------------------------------------------------------!
! Vegetated soil.
!----------------------------------------------------------------------!
if (process_vege) then
  fh (1,2) = -fhsng (2) * fr_snow (2) - fhsng_scale (2) &
  &          + (-htdripw (2) &
  &          + evapvg * elh + thrm_soil (2) - thrm_can &
  &          * (1.0d0 - trans_sw) &
  &          - trans_sw * (srht + trht)) * (1.0d0 - fr_snow (2))
  !--------------------------------------------------------------------!
  ! Canopy
  !--------------------------------------------------------------------!
  fch (0) = -htpr + &
  &         (evapvw * elh * fw + snsh (2) + thrm_can &
  &         * (1.0d0 - trans_sw) &
  &         - (1.0d0 - trans_sw) * (srht + trht) &
  &         + evapvd * elh * fd) * (1.0d0 - fm * fr_snow (2))
  fch (1) = -(thrm_can - thrm_soil (2)) * (1.0d0 - trans_sw) &
  &         * (1.0d0 - fr_snow (2)) & !rad soil
  &         -(thrm_can - thrmsn (2)) * fr_snow (2) * &
  &         (1.d0 - fm) &   
  !rad snow
  &        * (1.0d0 - trans_sw) &
  &        - htdripw (2) - htdrips (2)                  !heat of precip
  !--------------------------------------------------------------------!
endif
!----------------------------------------------------------------------!
!XXXXXXXXXXXXXXXX don't forget to add trans_sw stuff to snow !
!----------------------------------------------------------------------!
! Compute total LW heat (J/s).
!----------------------------------------------------------------------!
thrm_tot (1) = thrm_soil (1) * (1.0d0 - fr_snow (1)) + &
&              thrmsn (1) * fr_snow (1)
thrm_tot (2) = thrm_can * (1.0d0 - fr_snow (2) * fm) * &
&              (1.d0 - trans_sw) + &
&              thrmsn (2) * fr_snow (2) * fm + &
&              thrmsn (2) * trans_sw * fr_snow (2) * (1.0d0 - fm) + &
&              trans_sw * thrm_soil (2) * (1.0d0 - fr_snow (2))
!----------------------------------------------------------------------!
! Compute total sensible heat (J/s).
!----------------------------------------------------------------------!
snsh_tot (1) = snsh (1) * (1.0d0 - fr_snow (1)) + snshs (1) * &
&              fr_snow (1)
snsh_tot (2) = snsh (2) * (1.0d0 - fr_snow (2) * fm) + snshs (2) * &
&              fr_snow (2) * fm
!----------------------------------------------------------------------!
end subroutine flhg
!======================================================================!
subroutine runoff
!----------------------------------------------------------------------!
! Calculates surface and underground runoffs.
!----------------------------------------------------------------------!
! input:
!   xinfc  - infiltration capacity, m s-1
!   prfr   - fraction of precipitation
!   xk     - conductivity, m s-1
!   dz     - layer thicknesses, m
!   sl     - slope
!   sdstnc - interstream distance, m
!----------------------------------------------------------------------!
! output:
!   rnf  - surface runoff
!   rnff - underground runoff, m s-1
!----------------------------------------------------------------------!
! use effects of subgrid scale rain
! use precipitation that includes smow melt
!----------------------------------------------------------------------!
use soil_variables
use veg_variables
use atmos_variables
use constants
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
! Looks like these are coefficients for Topmodel. Leave for now in case
! use this at a later stage.
!----------------------------------------------------------------------!
real*8 f_k0_exp_k
!----------------------------------------------------------------------!
! water_down: flux of water at the soil surface.
! satfrac: fraction of saturated soil.
! prec_fr: soil fraction in which precipitation is falling.
!----------------------------------------------------------------------!
real*8 water_down, satfrac, prec_fr
!----------------------------------------------------------------------!
! sdstnc: interstream distance (m).
!----------------------------------------------------------------------!
real*8, parameter :: sdstnc = 100.d0
!----------------------------------------------------------------------!
! rosmp: used to compute saturated fraction ((w/ws)**rosmp).
!----------------------------------------------------------------------!
real*8, parameter :: rosmp = 8.0d0
!----------------------------------------------------------------------!
! Zero out underground runoff for all surfaces and layers (m/s).
!----------------------------------------------------------------------!
rnff (:,:) = 0.d0
!----------------------------------------------------------------------!
! Surface runoff hack to conserve water for ibv !=0,1 (m/s).
! Should be set to 0 after testing.
!----------------------------------------------------------------------!
rnf (:) = pr
!----------------------------------------------------------------------!
! Compute surface runoff for each surface (m/s).
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  water_down = -f (1,ibv)
  water_down = max (water_down, zero) ! to make sure rnf > 0.
  !--------------------------------------------------------------------!
  ! Everything that falls on saturated fraction goes to runoff.
  !--------------------------------------------------------------------!
  satfrac = (w (1,ibv) / ws (1,ibv)) ** rosmp
  rnf (ibv) = satfrac * water_down
  water_down = (1.0d0 - satfrac) * water_down
  !--------------------------------------------------------------------!
  ! if we introduce large scale precipitation it should be
  ! applied here
  !!! the following line is a hack. in a more precise approach
  !!! one should treat snow-free fraction separately
  !--------------------------------------------------------------------!
  prec_fr = max (prfr, fr_snow (ibv))
  !--------------------------------------------------------------------!
  if (water_down * 30.0d0 > xinfc (ibv) * prec_fr) then
    rnf (ibv) = rnf (ibv) + water_down * exp (-xinfc (ibv) * &
    &           prec_fr / water_down)
  endif
enddo
!----------------------------------------------------------------------!
! Compute underground runoff for each surface and layer (m/s).
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  do nsoil_layer = 1, nsoil_layers
    rnff (nsoil_layer,ibv) = xku (nsoil_layer,ibv) * sl * &
    &                        dz (nsoil_layer) / sdstnc
  end do
end do
!----------------------------------------------------------------------!
end subroutine runoff
!======================================================================!
subroutine fllmt
!----------------------------------------------------------------------!
use control_parameters
use veg_variables
use soil_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
real*8 :: trunc
!real*8 :: trunc = 1.0d-12 ! Fix for truncation
real*8 :: wn    ! New w to test saturation
real*8 :: dflux ! Water that can be taken from layer to comp. -rnf
real*8 :: drnf
trunc = 0.0d0 ! Fix for truncation.
!----------------------------------------------------------------------!
! Prevent over/underestimation of layers 2-n.
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  do nsoil_layer = nsoil_layers, 2, -1
    wn = w (nsoil_layer,ibv) + (f (nsoil_layer+1,ibv) - &
    &    f (nsoil_layer,ibv) - rnff (nsoil_layer,ibv) - &
    &    fd * (1.0d0 - fr_snow (2) * fm) * &
    &    evapdl (nsoil_layer,ibv)) * dts
    !-------------------------------------------------------------------!
    ! Compensate oversaturation by increasing flux up
    !-------------------------------------------------------------------!
    if (wn - ws (nsoil_layer,ibv) > trunc) &
    & f (nsoil_layer,ibv) = f(nsoil_layer,ibv) + (wn - &
    & ws (nsoil_layer,ibv) + trunc) / dts
    !-------------------------------------------------------------------!
    ! Compensate undersaturation by decreasing runoff
    !-------------------------------------------------------------------!
    if (wn - dz (nsoil_layer) * thetm (nsoil_layer,ibv) < trunc) then
      rnff (nsoil_layer,ibv) = rnff (nsoil_layer,ibv) + &
      & (wn - dz (nsoil_layer) * thetm (nsoil_layer,ibv) - trunc) / dts
      !-----------------------------------------------------------------!
      ! See if have to compensate with f
      !-----------------------------------------------------------------!
      if (rnff (nsoil_layer,ibv) < 0.0d0 ) then
        f (nsoil_layer,ibv) = f (nsoil_layer,ibv) + &
        &  rnff (nsoil_layer,ibv)
        rnff (nsoil_layer,ibv) = 0.0d0
      endif
      !-----------------------------------------------------------------!
    endif
  enddo
enddo
!----------------------------------------------------------------------!
! Prevent over/undersaturation of first layer.
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  wn = w (1,ibv) + (f (2,ibv) - f (1,ibv) &
  & - rnf (ibv) - rnff (1,ibv) &
  & - fd * (1.0d0 - fr_snow (2) * fm) * evapdl (1,ibv)) * dts
  if (wn - ws (1,ibv) > trunc) &
  & rnf (ibv) = rnf (ibv) + (wn - ws (1,ibv) + trunc) / dts
  if (wn - dz (1) * thetm (1,ibv) < trunc) &
  & rnf (ibv) = rnf (ibv) + &
  & (wn - dz (1) * thetm (1,ibv) - trunc) / dts
enddo
!----------------------------------------------------------------------!
! Now trying to remove negative runoff.
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  nsoil_layer = 1
  do while ((rnf (ibv) < 0.0d0) .and. (nsoil_layer <= nsoil_layers))
    if (nsoil_layer > 1) then
      !----------------------------------------------------------------!
      ! This is how much water we can take from layer nsoil_layer.
      !----------------------------------------------------------------!
      dflux = f (nsoil_layer+1,ibv) + &
      & (w (nsoil_layer,ibv) - dz (nsoil_layer) * &
      & thetm (nsoil_layer,ibv)) / dts &
      & - f (nsoil_layer,ibv) - rnff (nsoil_layer,ibv) &
      & - fd * (1.0d0 - fr_snow (2) * fm) * evapdl (nsoil_layer,ibv)
      !----------------------------------------------------------------!
      f (nsoil_layer,ibv) = f (nsoil_layer,ibv) - rnf (ibv)
      rnf (ibv) = rnf (ibv) + min (-rnf (ibv), dflux)
    endif
    !-------------------------------------------------------------------!
    ! rnff always >= 0, use it also to compensate rnf < 0.
    !-------------------------------------------------------------------!
    if (rnff (nsoil_layer,ibv) < 0.0d0) then
      write (98,*) 'fllmt: negative underground runoff' 
      stop
    endif
    drnf = min (-rnf (ibv), rnff (nsoil_layer,ibv))
    rnf (ibv) = rnf (ibv) + drnf
    rnff (nsoil_layer,ibv) = rnff (nsoil_layer,ibv) - drnf
    nsoil_layer = nsoil_layer + 1
  enddo
  !---------------------------------------------------------------------!
  ! Check if rnf == 0 up to machine accuracy.
  !---------------------------------------------------------------------!
  if (rnf (ibv) < -1.0d-12 ) then
    write (98,*) 'fllmt: rnf<0, ibv=',ibv,rnf(ibv)
    write (98,*) 'fllmt: negative runoff'
    stop
  endif
  !---------------------------------------------------------------------!
  ! If -1d-12 < rnf < 0. put it to 0 to avoid possible problems
  ! actually for ground hydrology it is not necessary.
  !---------------------------------------------------------------------!
  rnf (ibv) = max (rnf(ibv), 0.0d0)
enddo
!----------------------------------------------------------------------!
end subroutine fllmt
!======================================================================!
subroutine xklh
!----------------------------------------------------------------------!
! Evaluates the heat conductivity between layers using method of
! Devries.
!----------------------------------------------------------------------!
use control_parameters
use soil_variables
use veg_variables
use constants
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
integer :: i,j,k
real*8 :: gaa,gabc(3),gca,xa,xb,xden,xi,xnum,xs,xw
real*8 :: ba,hcwt(ntextures-1),hcwta,hcwtb,hcwti,hcwtw
real*8 :: xsha(nsoil_layers_max,2),xsh(nsoil_layers_max,2)
!----------------------------------------------------------------------!
! The alam''s are the heat conductivities.
!----------------------------------------------------------------------!
real*8, parameter :: alamw = 0.573345d0 &
&     ,alami = 2.1762d0 &
&     ,alama = 0.025d0 &
&     ,alambr= 2.9d0 &
&     ,alams (ntextures-1) = (/ 8.8d0, 2.9d0, 2.9d0, .25d0 /)
!----------------------------------------------------------------------!
hcwtw = 1.0d0
hcwti=0.d0
hcwtb=1.d0
!----------------------------------------------------------------------!
do i=1,ntextures-1
  hcwt(i)=0.d0
end do
do j=1,3
  hcwti=hcwti+1.d0/(1.d0+(alami/alamw-1.d0)*gabc(j))
  do i=1,ntextures-1
    hcwt(i)=hcwt(i)+1.d0/(1.d0+(alams(i)/alamw-1.d0)*gabc(j))
  end do
end do
hcwti=hcwti/3.d0
do i=1,ntextures-1
  hcwt(i)=hcwt(i)/3.d0
end do
do ibv=1,2        ! i_bare,i_vege
  do k=1,nsoil_layers
    xsha(k,ibv)=0.d0
    xsh(k,ibv)=0.d0
    do i=1,ntextures-1
      xs=(1.d0-thm(0,i))*q(i,k)
      xsha(k,ibv)=xsha(k,ibv)+xs*hcwt(i)*alams(i)
      xsh(k,ibv)=xsh(k,ibv)+xs*hcwt(i)
    end do
  end do
end do
ba=alama/alamw-1.d0
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  do nsoil_layer = 1, nsoil_layers
    gaa = 0.298d0 * theta (nsoil_layer,ibv) / &
    &     (thets (nsoil_layer,ibv) + 1.0d-6) + 0.035d0
    gca = 1.0d0 - 2.0d0 * gaa
    hcwta = (2.0d0 / (1.0d0 + ba * gaa) + &
    &       1.0d0 / (1.0d0 + ba * gca)) / 3.0d0
    ! xw,xi,xa are the volume fracs. don''t count snow in soil lyr 1.
    xw = w (nsoil_layer,ibv) * (1.0d0 - fice (nsoil_layer,ibv)) / &
    &    dz (nsoil_layer)
    xi = w (nsoil_layer,ibv) * fice (nsoil_layer,ibv) / dz (nsoil_layer)
    xa = (thets (nsoil_layer,ibv) - theta (nsoil_layer,ibv))
    xb = q (ntextures,nsoil_layer)
    xnum = xw * hcwtw * alamw + xi * hcwti * alami + &
    &      xa * hcwta * alama + xsha (nsoil_layer,ibv) &
    &      + xb * hcwtb * alambr
    xden = xw * hcwtw + xi * hcwti + xa * hcwta + &
    &      xsh (nsoil_layer,ibv) + xb * hcwtb
    xkh (nsoil_layer,ibv) = xnum / xden
    if (xkh (nsoil_layer,ibv) < 0.0d0 ) then
      write (99,*) 'xklh: heat conductivity<0'
      stop
    endif
  enddo
enddo
!----------------------------------------------------------------------!
! Get the average conductivity between layers
!----------------------------------------------------------------------!
do ibv =i_bare, i_vege
  do nsoil_layer = 2, nsoil_layers
    xkhm (nsoil_layer,ibv) = ((zb (nsoil_layer)- zc (nsoil_layer-1)) &
    &  * xkh (nsoil_layer,ibv) + (zc (nsoil_layer) - &
    &  zb (nsoil_layer)) * xkh (nsoil_layer-1,ibv)) / &
    &  (zc (nsoil_layer) - zc (nsoil_layer - 1))
  enddo
enddo
!----------------------------------------------------------------------!
end subroutine xklh
!======================================================================!
subroutine fl
!----------------------------------------------------------------------!
! Computes water fluxes between soil layers.
! Input:
!   h     ; water potential         (m)
!   xk    ; conductivity          (m/s)
!   zc    ; layer centres           (m)
! Output:
!   f     ; flux upwards          (m/s)
!   xinfc ; infiltration capacity (m/s)
!----------------------------------------------------------------------!
use control_parameters
use soil_variables
use veg_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
! Set no flux across bottom of profile.
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  f (nsoil_layers+1,ibv) = 0.0d0
enddo
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  do nsoil_layer = 2, nsoil_layers
    f (nsoil_layer,ibv) = -xk (nsoil_layer,ibv) *                      &
    &                     (h (nsoil_layer-1,ibv) - h (nsoil_layer,ibv))&
    &                     / (zc (nsoil_layer-1) - zc (nsoil_layer))
  enddo
enddo
!----------------------------------------------------------------------!
! Put infiltration maximum into xinfc.
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  xinfc (ibv) = xk (1,ibv) * h (1,ibv) / zc (1)
enddo
!----------------------------------------------------------------------!
end subroutine fl
!======================================================================!
subroutine flh
!----------------------------------------------------------------------!
use control_parameters
use veg_variables
use soil_variables
use constants
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
do ibv= i_bare, i_vege
  fh (nsoil_layers+1,ibv) = 0.0d0
  ! Total heat flux is heat carried by water flow plus heat conduction.
  do nsoil_layer = 2, nsoil_layers
    fh (nsoil_layer,ibv) = -xkhm (nsoil_layer,ibv) * &
    &  (tp (nsoil_layer-1,ibv) - tp (nsoil_layer,ibv)) / &
    &  (zc (nsoil_layer-1) - zc (nsoil_layer))
    if (f (nsoil_layer,ibv) > 0.0d0) then
      fh (nsoil_layer,ibv) = fh (nsoil_layer,ibv) + &
      &  f (nsoil_layer,ibv) * tp (nsoil_layer,ibv) * shw
    else
      fh (nsoil_layer,ibv) = fh (nsoil_layer,ibv) + &
      &  f (nsoil_layer,ibv) * tp (nsoil_layer-1,ibv) * shw
    endif
  enddo
enddo
!----------------------------------------------------------------------!
end subroutine flh
!======================================================================!
subroutine retp
!----------------------------------------------------------------------!
! Evaluates temperatures in soil layers based on the heat values and
! executes snowmelt.
!----------------------------------------------------------------------!
! Inputs:
!----------------------------------------------------------------------!
! w  : water in soil layers                                          (m)
! ht : heat in soil layers                                       (J/m^2)
!----------------------------------------------------------------------!
! Output:
!----------------------------------------------------------------------!
! tp : temperatures of layers                                    (deg C)
!----------------------------------------------------------------------!
use constants
use control_parameters
use soil_variables
use veg_variables
use atmos_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
integer :: kk          ! Index of first layer in surface type        (n)
!----------------------------------------------------------------------!
! Place holder for heat in snow layers.
!----------------------------------------------------------------------!
hsn (:,:) = 0.0d0
!----------------------------------------------------------------------!
! Initialise soil temperatures and ice fractions.
!----------------------------------------------------------------------!
tp   (:,:) = 0.0d0
fice (:,:) = 0.0d0
!----------------------------------------------------------------------!
! Loop over surfaces.
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  ! Set index of first layer in surface type
  ! (to get canopy if vegetated).
  kk = 2 - ibv
  do nsoil_layer = kk, nsoil_layers
    ! If all frozen.
    if (fsn * w (nsoil_layer,ibv) + ht (nsoil_layer,ibv) < 0.d0 ) then
      tp (nsoil_layer,ibv) = (ht (nsoil_layer,ibv) + &
      & w (nsoil_layer,ibv) * fsn) / (shc (nsoil_layer,ibv) &
      & + w (nsoil_layer,ibv) * shi)
      fice (nsoil_layer,ibv) = 1.0d0
    elseif (ht (nsoil_layer,ibv) > 0.0d0) then ! all melted
      tp (nsoil_layer,ibv) = ht (nsoil_layer,ibv) / &
      & (shc (nsoil_layer,ibv) + w (nsoil_layer,ibv) * shw)
    elseif (w (nsoil_layer,ibv) >= 1.0d-12) then  ! part frozen
      fice (nsoil_layer,ibv) = -ht (nsoil_layer,ibv) / &
      & (fsn * w (nsoil_layer,ibv))
    endif
  enddo
enddo
!----------------------------------------------------------------------!
! Fix for now. tsn1 setting should be moved elsewhere.
!----------------------------------------------------------------------!
tsn1 (:) = 0.0d0
do ibv = i_bare, i_vege
  if (wsn (1,ibv) > 1.0d-6 .and. &
  &   hsn (1,ibv) + wsn (1,ibv) * fsn < 0.0d0) then
      tsn1 (ibv) = (hsn (1,ibv) + wsn (1,ibv) * fsn) / &
  &   (wsn (1,ibv) * shi)
  endif
enddo
if (tp (1,1) > 100.0d0 .or. tp (0,2) > 100.0d0) then
  write (*,*) 'retp tp bounds error',tp(1,1),tp(0,2)
  write (*,*) 'i0 j0 Ch',i0,j0,Ch
  write (*,*) 't ',t(i0,j0)
  write (99,*) 'retp tp bounds error',tp(1,1),tp(0,2)
  write (99,*) 'i0 j0 Ch',i0,j0,Ch
  write (99,*) 't ',t(i0,j0)
  stop
endif
!----------------------------------------------------------------------!
end subroutine retp
!======================================================================!
subroutine apply_fluxes
!----------------------------------------------------------------------!
! Apply computed fluxes to advance w.
!----------------------------------------------------------------------!
use control_parameters
use veg_variables
use soil_variables
use constants
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
! The canopy.
!----------------------------------------------------------------------!
w  (0,2) = w  (0,2) + (fc  (1) - fc  (0)) * dts
ht (0,2) = ht (0,2) + (fch (1) - fch (0)) * dts
!----------------------------------------------------------------------!
! Carbon.
!----------------------------------------------------------------------!
 cland (1) = cland (1) + (gpp - ra - cgrow)            * dts
 cland (2) = cland (2) + (cgrow - clit)                * dts
 cland (3) = cland (3) + (clit - rh_clit - clit_trans) * dts
 cland (4) = cland (4) + (clit_trans - rh_som)         * dts
!----------------------------------------------------------------------!
! The soil water and heat.
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  !--------------------------------------------------------------------!
  ! Surface runoff.
  !--------------------------------------------------------------------!
  w  (1,ibv) = w (1,ibv) - rnf (ibv) * dts
  ht (1,ibv) = ht (1,ibv) - shw * &
  &            max (tp(1,ibv),0.0d0) * rnf (ibv) * dts
  !--------------------------------------------------------------------!
  ! Rest of fluxes.
  !--------------------------------------------------------------------!
  do nsoil_layer = 1, nsoil_layers
    w (nsoil_layer,ibv) = w (nsoil_layer,ibv) +        &
    & (f (nsoil_layer+1,ibv) - f (nsoil_layer,ibv) - &
    & rnff (nsoil_layer,ibv)                           &
    & - fd * (1.0d0 - fr_snow (2) * fm) * evapdl (nsoil_layer,ibv) &
    & ) * dts
    ht (nsoil_layer,ibv) = ht (nsoil_layer,ibv) + &
    &  (fh (nsoil_layer+1,ibv) - fh (nsoil_layer,ibv) - &
    &   shw * max (tp(nsoil_layer,ibv),0.0d0) * &
    &  (rnff (nsoil_layer,ibv))) * dts
  enddo
enddo
!----------------------------------------------------------------------!
! Check for under/over-saturation.
!----------------------------------------------------------------------!
do ibv = i_bare, i_vege
  do nsoil_layer = 1, nsoil_layers
    if (w (nsoil_layer,ibv) < dz (nsoil_layer) * &
    & thetm (nsoil_layer,ibv) - 1.0d-14) then
      write (98,*) 'ghy:',nsoil_layer,ibv,w(nsoil_layer,ibv), &
      &            dz(nsoil_layer),thetm(nsoil_layer,ibv)
      write (98,*) 'ghy: w < dz*thetm'
      stop
    endif
    if (w (nsoil_layer,ibv) > ws (nsoil_layer,ibv) + 1.0d-14) then
    !if (w (nsoil_layer,ibv) > ws (nsoil_layer,ibv) + 1.0d-12) then
      write (98,*) 'ghy: w > ws'
      stop
    endif
    w (nsoil_layer,ibv) = max (w (nsoil_layer,ibv), &
    & dz (nsoil_layer) * thetm (nsoil_layer,ibv))
    w (nsoil_layer,ibv) = min (w (nsoil_layer,ibv), &
    & ws (nsoil_layer,ibv))
  enddo
enddo
!----------------------------------------------------------------------!
end subroutine apply_fluxes
!======================================================================!
subroutine advnc
!----------------------------------------------------------------------!
! Advance quantities by one timestep.
!----------------------------------------------------------------------!
use veg_variables
use soil_variables
use atmos_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
real*8 :: dum1,dum2
!----------------------------------------------------------------------!
! Make sure that dts is always initialised.
!----------------------------------------------------------------------!
nit = 0
dtr = dt
dts = dt
!----------------------------------------------------------------------!
! Reset main water fluxes so that they are always initialsed.
!----------------------------------------------------------------------!
f   (:,:)  = 0.0d0 ! Water flux between layers (>0 is up) (m/s)
fh  (:,:)  = 0.0d0 ! Heat flux between layers (>0 is up) (W)
fc  (:)    = 0.0d0 ! Water flux at canopy boundaries (>0 is up) (m/s)
fch (:)    = 0.0d0 ! Heat flux at canopy boundaries (>0 is up) (W)
!----------------------------------------------------------------------!
! Normal case (both present).
i_bare = 1
i_vege = 2
process_bare = .true.
process_vege = .true.
if (fb == 0.0d0) then ! bare fraction is missing
  i_bare = 2
  process_bare = .false.
endif
if (fv == 0.0d0) then ! vege fraction is missing
  i_vege = 1
  process_vege = .false.
endif
call reth ! theta; snowd; fw,fd,fm
call retp ! tp; fice
!----------------------------------------------------------------------!
trnf  (:) = 0.0d0 ! Total surface runoff (m).
trnff (:) = 0.0d0 ! Total underground runoff (m).
trnfe (:) = 0.0d0 ! Total surface runoff energy (J).
trnffe(:) = 0.0d0 ! Total underground runoff energy (J).
aevapb    = 0.0d0 ! Bare soil evaporation (m/s).
aevapw    = 0.0d0 ! Wet canopy evaporation (m/s).
aevapd    = 0.0d0 ! Total dry canopy evaporation (m).
aevapvg   = 0.0d0 ! Soil under canopy evaporation (m/s).
asnsh_tot (1:2) = 0.0d0 ! Total sensible heat flux (J).
athrm_tot (1:2) = 0.0d0 ! Total LWu heat flux (J).
agpp      = 0.0d0 ! Total gross primary production (kgC/m^2).
anpp      = 0.0d0 ! Total net primary production (kgC/m^2).
anee      = 0.0d0 ! Total net exchange with atmos. (+ve up) (kgC/m^2).
!----------------------------------------------------------------------!
do while (dtr > 0.0d0)
  nit = nit + 1
  if (nit > limit) then
    write (*,*) 'Max. no. iterations exceeded in advnc'
    write (98,*) 'Max. no. iterations exceeded in advnc'
    stop
  endif
  call hydra
  ! Compute water potentials of soil layers
  call evap_limits (.true.,dum1,dum2) ! Computes evapdl
  call sensible_heat
  call xklh
  call gdtm
  if (dtm >= dtr) then
    dts = dtr
    dtr = 0.0d0
  else
    dts = min (dtm, dtr * 0.5d0)
    dtr = dtr - dts
  endif
!write (88,*)jday,dts
  call drip_from_canopy ! Compute drip of water from canopy.
  call fl               ! Compute water fluxes between soil layers
  call flg              ! Compute water fluxes at surface
  call runoff           ! Compute surface and underground runoffs
  call simveg           ! Compute vegetation fluxes
  call fllmt            ! Place limits on soil water fluxes
  call flh              ! Evaluate heat flux between layers.
  call flhg             ! Calculate ground heat fluxes (to surface)
  call apply_fluxes     ! Advance w and ht
  !--------------------------------------------------------------------!
  ! Diagnostics to check water balance.
  !--------------------------------------------------------------------!
  aevapb = aevapb + evapb * dts       ! Bare soil evaporation (m/s)
  aevapw = aevapw + evapvw * fw * dts ! Wet canopy evaporation (m/s)
  aevapd = aevapd + evapvd * fd * dts ! Total dry canopy evaporation (m)
  aevapvg = aevapvg + evapvg * dts ! Soil under canopy evap.   (m/s)
  !--------------------------------------------------------------------!
  ! Diagnostics to check energy balance.
  !--------------------------------------------------------------------!
  asnsh_tot (1) = asnsh_tot (1) + snsh_tot (1) * dts
  asnsh_tot (2) = asnsh_tot (2) + snsh_tot (2) * dts
  athrm_tot (1) = athrm_tot (1) + thrm_tot (1) * dts
  athrm_tot (2) = athrm_tot (2) + thrm_tot (2) * dts
  !--------------------------------------------------------------------!
  ! Carbon diagnostics.
  !--------------------------------------------------------------------!
  agpp = agpp + fv * gpp * dts
  anpp = anpp + fv * npp * dts
  anee = anee + fv * nee * dts
  !--------------------------------------------------------------------!
  !call reth             ! Compute new thetas
  !call retp             ! Compute new temperatures
  !call update_veg_locals
  do ibv = i_bare, i_vege
    ! Total surface runoff (m).
    trnf  (ibv) = trnf  (ibv) + rnf  (ibv) * dts
    ! Total surface runoff energy (J).
    trnfe (ibv) = trnfe (ibv) + shw * max (tp(1,ibv),0.0d0) * &
    &             rnf  (ibv) * dts
  enddo
  !--------------------------------------------------------------------!
  ! Total underground runoff (m).
  !--------------------------------------------------------------------!
  do ibv = i_bare, i_vege
    do nsoil_layer = 1, nsoil_layers
      trnff (ibv) = trnff (ibv) + rnff (nsoil_layer,ibv) * dts
      trnffe (ibv) = trnffe (ibv) + shw * &
      &              max (tp (nsoil_layer,ibv),0.0d0) * &
      &              rnff  (nsoil_layer,ibv) * dts
    enddo
  enddo
  !--------------------------------------------------------------------!
  call reth ! Compute new thetas
  call retp ! Compute new temperatures
  call update_veg_locals
  !--------------------------------------------------------------------!
  ! Skin temperatures (degC).
  !--------------------------------------------------------------------!
  tearth (i0,j0) = fb * tp (1,1) + fv * &
  & ((trans_sw * tp (1,2)) + (1.0d0 - trans_sw) * tp (0,2))
  !--------------------------------------------------------------------!
enddo
call hydra
!----------------------------------------------------------------------!
end subroutine advnc

      subroutine gdtm
!**** calculates the maximum time step allowed by stability
!**** considerations.
!cc   include 'soils45.com'
!**** soils28   common block     9/25/90
      use atmos_variables
      use veg_variables
      use soil_variables
      implicit none
      real*8 ak1,ak2(2),betas(2),xk2(2),dldz2,dqdt,sgmm
      real*8 dtm1,dtm2,dtm3,dtm4,t450,xk1
      real*8 :: dqsatdt,qsat_modele
      integer k,n
      n = nsoil_layers
      t450=450.d0
      dqdt=dqsatdt(ts,lhe)*qsat_modele(ts,lhe,ps)
!****
!**** first calculate timestep for water movement in soil.
      sgmm=1.0d0
      dldz2=0.d0
      do ibv=i_bare,i_vege
        do k=1,n
          dldz2=max(dldz2,d(k,ibv)/dz(k)**2)
        end do
      end do
      dtm=sgmm/(dldz2+1d-12)
      ! Following commented out 9/7/12 to try and speed things up.
      if(q(4,1).gt.0.d0)dtm=min(dtm,t450) !**********???
      dtm1=dtm
      if ( dtm .lt. 0.d0 ) then
        write (*,*) 'gdtm: dt1_ghy<0'
        write (99,*) 'gdtm: dt1_ghy<0'
        write (99,*) 'i0,j0',i0,j0
        stop
      endif
!****
!**** next calculate timestep for heat movement in soil.
      do ibv=i_bare,i_vege
        do k=1,n
          xk1=xkh(k,ibv)
          ak1=(shc(k,ibv)+((1.d0-fice(k,ibv))*shw+fice(k,ibv)*shi) &
     &         *w(k,ibv))/dz(k)
          dtm=min(dtm,.5d0*ak1*dz(k)**2/(xk1+1d-12))
        end do
      end do
      dtm2=dtm
      if ( dtm .lt. 0.d0 ) then
        write (*,*) 'gdtm: dt2_ghy<0'
        write (99,*) 'gdtm: dt2_ghy<0'
        write (99,*) 'i0,j0',i0,j0
        stop
      endif
!****
!**** finally, calculate max time step for top layer bare soil
!**** and canopy interaction with surface layer.
!**** use timestep based on coefficient of drag
      cna=ch*vsm
      !***rho3=.001d0*rho
      betas(1:2) = 1.d0 ! it''s an overkill but it makes the things
                        ! simpler.
      if(epb.le.0.d0)then
       betas(1)=1.0d0
      else
       betas(1)=evapb/epb
      endif
      if(epv.le.0.d0)then
       betas(2)=1.0d0
      else
       betas(2)=(evapvw*fw+evapvd*(1.d0-fw))/epv
      endif
      do ibv=i_bare,i_vege
        k=2-ibv
        xk2(ibv)=sha*rho*cna &
     &       + betas(ibv)*rho3*cna*elh*dqdt &
     &       + 8.d0*stbo*(tp(k,ibv)+tf)**3
        ak2(ibv)=shc(k,ibv)+((1.d0-fice(k,ibv))*shw+fice(k,ibv)*shi) &
     &       *w(k,ibv)
        dtm=min(dtm,0.5*ak2(ibv)/(xk2(ibv)+1d-12))
        if(ibv.eq.1)dtm3=dtm
        if(ibv.eq.2)dtm4=dtm
!
! prevent oscillation of top snow layer
!     if(isn(ibv).ne.0.or.snowd(ibv).ne.0.d0)then
!      ak3(ibv)=.05d0*shi*spgsn
!      dtm=min(dtm,ak3(ibv)/(xk2(ibv)+1.d-12))
!      if(ibv.eq.1)dtm5=dtm
!      if(ibv.eq.2)dtm6=dtm
!     endif
      end do
      if(dtm.lt.5.d0)then
       write(99,*) '*********** gdtm: ijdebug,fb,fv',fb,fv
       write(99,*) 'i0,j0',i0,j0
       write(99,*)'dtm',dtm1,dtm2,dtm3,dtm4
       write(99,*)'xk2',xk2
       write(99,*)'ak2',ak2
       write(99,*)'snsh',snsh
       write(99,*)'xlth',elh*evap_tot(1:2)
       write(99,*)'dqdt',dqdt
       write(99,*)'ts,tf',ts,tf
       write(99,*)'dlt',tp(1,1)-ts+tf,tp(0,2)-ts+tf
       write (*,*) "gdtm: time step < 5 s"
       write (99,*) "gdtm: time step < 5 s"
       stop
      endif
!****
      end subroutine gdtm
!======================================================================!
      FUNCTION QSAT_modele (TM,LH,PR)
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
!@ver  1.0
      USE constants, only : mrat,rvap,tf
      IMPLICIT NONE
!@var A,B,C   expansion coefficients for QSAT
      REAL*8, PARAMETER :: A=6.108d0*MRAT    !3.797915d0
      REAL*8, PARAMETER :: B= 1./(RVAP*TF)   !7.93252d-6
      REAL*8, PARAMETER :: C= 1./RVAP        !2.166847d-3
!**** Note that if LH is considered to be a function of temperature, the
!**** correct argument in QSAT is the average LH from t=0 (C) to TM, ie.
!**** LH = 0.5*(LH(0)+LH(t))
      REAL*8, INTENT(IN) :: TM  !@var TM   temperature (K)
      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
      REAL*8, INTENT(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
      REAL*8 :: QSAT_modele     !@var QSAT sat. vapour mixing ratio
      QSAT_modele = A*EXP(LH*(B-C/max(130.d0,TM)))/PR
      RETURN
      END
!======================================================================!

      FUNCTION DQSATDT (TM,LH)
!@sum  DQSATDT calculates change of sat. vapour mixing ratio with temp.
!@auth Gary Russell
!@ver  1.0
!**** Note that d(qsat)/dt = qsat * lh * c / T*T
!**** Only the factor of qsat is given here
      USE CONSTANTs, only : rvap
      IMPLICIT NONE
!@var C coefficient for QSAT
      REAL*8, PARAMETER :: C = 1./RVAP        !2.166847d-3
!**** Note that if LH is considered to be a function of temperature, the
!**** correct argument in DQSATDT is the actual LH at TM i.e. LH=LH(TM)
      REAL*8, INTENT(IN) :: TM  !@var TM   temperature (K)
      REAL*8, INTENT(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
      REAL*8 :: DQSATDT         !@var DQSATDT d(qsat)/dT factor only.
      DQSATDT = LH*C/(TM*TM)    ! * QSAT(TM,LH,PR)
      RETURN
      END
