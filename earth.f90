!======================================================================!
subroutine earth
!----------------------------------------------------------------------!
! Make surface calculations for each timestep.
!----------------------------------------------------------------------!
use control_parameters
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
if (sing) then
  !--------------------------------------------------------------------!
  do j0 = olat, olat
    do i0 = olon, olon
      call proc_surf
    enddo
  enddo
  !--------------------------------------------------------------------!
else
  !--------------------------------------------------------------------!
  ! Loop over global grid to compute fluxes.
  !--------------------------------------------------------------------!
  do j0 = 1, jm ! Loop over longitudes.
  !--------------------------------------------------------------------!
    do i0 = 1, im ! Loop over latitudes.
    !------------------------------------------------------------------!
      call proc_surf
    !------------------------------------------------------------------!
    enddo ! End of latitude loop (do i0 = 1, im).
  !--------------------------------------------------------------------!
  enddo ! End of longitude loop (do j0 = 1, jm).
  !--------------------------------------------------------------------!
endif
!----------------------------------------------------------------------!
end subroutine earth
!======================================================================!
subroutine proc_surf
!----------------------------------------------------------------------!
! Surface calculations for each timestep.
!----------------------------------------------------------------------!
use global_variables
use soil_variables
use atmos_variables
use constants
use veg_variables
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
real*8 :: pearth     ! Fraction of grid box with soil                (n)
real*8 :: q1         ! Water vap mix rat in first atmos l        (kg/kg)
real*8 :: tsv        ! Virtual potential temperature of the surface  (K)
real*8 :: qsrf       ! Atmosphere surface level sp. humidity     (kg/kg)
real*8 :: rhosrf     ! Air density of surface layer             (kg/m^3)
real*8 :: ma1        ! Mass of first atmosphere layer           (kg^m-2)
real*8 :: eprec      ! Energy in precipitation                   (J/m^2)
real*8 :: tprcp      ! Temperature of precipiation                   (K)
real*8    prcp       ! Amount of precipitation                  (kg/m^2)
real*8 :: waterb_old ! To test for cons. of water, bare fraction     (m)
real*8 :: waterv_old ! To test for cons. of water, vege fraction     (m)
real*8 :: waterb_new ! To test for cons. of water, bare fraction     (m)
real*8 :: waterv_new ! To test for cons. of water, vege fraction     (m)
real*8 :: wdiffb     ! To test for cons. of water, bare fraction     (m)
real*8 :: wdiffv     ! To test for cons. of water, vege fraction     (m)
real*8 :: energyb_old ! To test for cons. of energy, bare fraction   (J)
real*8 :: energyv_old ! To test for cons. of energy, vege fraction   (J)
real*8 :: energyb_new ! To test for cons. of energy, bare fraction   (J)
real*8 :: energyv_new ! To test for cons. of energy, vege fraction   (J)
real*8 :: ediffb      ! To test for cons. of energy, bare fraction   (J)
real*8 :: ediffv      ! To test for cons. of energy, vege fraction   (J)
!----------------------------------------------------------------------!
! Fraction of grid box with soil (n).
!----------------------------------------------------------------------!
pearth = fearth (i0,j0)
!----------------------------------------------------------------------!
! Atmospheric pressure at surface level (hPa).
!----------------------------------------------------------------------!
ps = pedn (i0,j0)
!----------------------------------------------------------------------!
!  Water vapour mixing ratio in first atmosphere layer (kg/kg).
!----------------------------------------------------------------------!
q1 = q_surf (i0,j0)
!----------------------------------------------------------------------!
! Mass of air at surface level (kg/m^2).
!----------------------------------------------------------------------!
ma1 = am (i0,j0)
!----------------------------------------------------------------------!
! Mass of water in air at surface level (kg/m^2).
!----------------------------------------------------------------------!
qm1 = q1 * ma1
!----------------------------------------------------------------------!
if (pearth > 0.0d0) then
  !--------------------------------------------------------------------!
  ! Local sine of solar elevation (ratio).
  !--------------------------------------------------------------------!
  sbeta_l   = sbeta_can (i0,j0)
  !--------------------------------------------------------------------!
  ! Local downwelling diffuse PAR (umol/m^2/s).
  !--------------------------------------------------------------------!
  ipar_df_l = ipar_df_can (i0,j0)
  !--------------------------------------------------------------------!
  ! Local downwelling direct PAR (umol/m^2/s).
  !--------------------------------------------------------------------!
  ipar_dr_l = ipar_dr_can (i0,j0)
  !--------------------------------------------------------------------!
  ! Local precipitation rate (m/s == 10^3 kg/m^2/s).
  !--------------------------------------------------------------------!
  prcp = prec (i0,j0) ! kg/m^2
  pr = prcp / (dtsrc * rhow) ! (kg/m^2/s)/1000
  tprcp = t (i0,j0) - tf             ! deg C
  ! Energy of precipitation (J/m^2).
  if (tprcp > 0.0) then
    eprec = prcp * tprcp * shw_kg
  else
    eprec = prcp * tprcp * shi_kg
  endif
  !--------------------------------------------------------------------!
  ! Heat in precipitation (W/m^2).
  !--------------------------------------------------------------------!
  htpr = eprec / dtsrc
  !--------------------------------------------------------------------!
  ! Bare and vegetated fractions (fraction)
  !--------------------------------------------------------------------!
  fb = afb (i0,j0)
  fv = 1.d0 - fb
  !--------------------------------------------------------------------!
  ! Normal case (both bare and vegetated surfaces present).
  !--------------------------------------------------------------------!
  i_bare = 1
  i_vege = 2
  process_bare = .true.
  process_vege = .true.
  if (fb == 0.0d0) then ! bare fraction is missing
    i_bare = 2
    process_bare = .false.
  endif
  if (fv == 0.0d0) then ! vegetated fraction is missing
    i_vege = 1
    process_vege = .false.
  endif
  !--------------------------------------------------------------------!
  if (process_bare) then
    waterb_old = sum (wbare (1:nsoil_layers_max,i0,j0)) + &
    &            prec (i0,j0) / (rhow * float (nisurf))
    energyb_old = sum (htbare (1:nsoil_layers_max,i0,j0)) + &
    &             htpr * dtsrc / float(nisurf) +        &
    &             fsf (i0,j0) * dtsrc / float(nisurf) +        &
    &             trhr (i0,j0) * dtsrc / float(nisurf)
  endif
  if (process_vege) then
    waterv_old = sum (wvege (0:nsoil_layers_max,i0,j0)) + &
    &            prec(i0,j0) / (rhow * float (nisurf))
    energyv_old = sum (htvege (0:nsoil_layers_max,i0,j0)) + &
    &             htpr * dtsrc / float(nisurf) +        &
    &             fsf (i0,j0) * dtsrc / float(nisurf) +        &
    &             trhr (i0,j0) * dtsrc / float(nisurf)
  endif
  !--------------------------------------------------------------------!
  ! Local fraction (by area) of precipitation (fraction).
  ! Seems fixed to 10% in GCM!
  !--------------------------------------------------------------------!
  prfr = 0.1d0
  !--------------------------------------------------------------------!
  ! Local soil water depths (m).
  !--------------------------------------------------------------------!
  w (1:nsoil_layers_max,1) = wbare (1:nsoil_layers_max,i0,j0)
  w (0:nsoil_layers_max,2) = wvege (0:nsoil_layers_max,i0,j0)
  !--------------------------------------------------------------------!
  ! Local soil heat values (J/m^2).
  !--------------------------------------------------------------------!
  ht (0:nsoil_layers_max,1) = htbare (0:nsoil_layers_max,i0,j0)
  ht (0:nsoil_layers_max,2) = htvege (0:nsoil_layers_max,i0,j0)
  !--------------------------------------------------------------------!
  ! Snow depths on each surface (m water equivalent).
  !--------------------------------------------------------------------!
  snowd (1:2) = snowbv (1:2,i0,j0)
  !--------------------------------------------------------------------!
  ! Number of snow layers on each surface (m water equivalent).
  !--------------------------------------------------------------------!
  nsnow_layers (1:2) = nsnow_layers_ij (1:2,i0,j0)
  !--------------------------------------------------------------------!
  ! Snow depths in each layer (m water equivalent).
  !--------------------------------------------------------------------!
  wsn (1:nsnow_layers_max,1:2) = wsn_ij (1:nsnow_layers_max,1:2,i0,j0)
  !--------------------------------------------------------------------!
  ! Fraction of snow on each surface (fraction).
  !--------------------------------------------------------------------!
  fr_snow (1:2) = fr_snow_ij (1:2,i0,j0)
  !--------------------------------------------------------------------!
  ! Incident SW (W/m^2). N.B. Already cosine corrected.
  !--------------------------------------------------------------------!
  srheat = fsf (i0,j0)
  !--------------------------------------------------------------------!
  ! Surface atmosphere layer specific humidity in modelE given by
  ! PBL. Here taken from weather input data (kg/kg).
  !--------------------------------------------------------------------!
  qsrf = q_surf (i0,j0)
  !--------------------------------------------------------------------!
  ! Surface atmosphere level specific humidity (kg/kg).
  !--------------------------------------------------------------------!
  qs = qsrf
  !--------------------------------------------------------------------!
  ! Virtual potential temperature of the surface (K).
  ! Just do inverse of what modelE does, assuming t here is ts there.
  !--------------------------------------------------------------------!
  ts = t (i0,j0) ! K
  tsv = ts * (1.0d0 + qs * deltx)
  !--------------------------------------------------------------------!
  ! Simply PBL scheme to compute atmospheric conductance.
  !--------------------------------------------------------------------!
  vsm = 2.0d0
  call pbl (ch,tearth(i0,j0)+tf)
  !--------------------------------------------------------------------!
  ! Surface level air density (kg/m^3).
  !--------------------------------------------------------------------!
  rhosrf = 100.0d0 * ps / (rgas * tsv)
  !--------------------------------------------------------------------!
  ! Surface level air density to be passed (kg/m^3).
  !--------------------------------------------------------------------!
  rho = rhosrf
  !--------------------------------------------------------------------!
  ! Incident SW (W/m^2).
  !--------------------------------------------------------------------!
  srht = srheat
  !--------------------------------------------------------------------!
  ! Incident LWd (W/m^2).
  !--------------------------------------------------------------------!
  trheat = trhr (i0,j0)
  trht = trheat
  !--------------------------------------------------------------------!
  ! Calculate soil properties.
  !--------------------------------------------------------------------!
  call ghinij
  !--------------------------------------------------------------------!
  ! Initialise local vegetation properties.
  !--------------------------------------------------------------------!
  call veg_set_cell
  !--------------------------------------------------------------------!
  ! Advance soil state variables by one timstep.
  !--------------------------------------------------------------------!
  call advnc
  !call evap_limits (.false.,evap_max_ij(i0,j0),fr_sat_ij(i0,j0))
  !--------------------------------------------------------------------!
  ! Save local vegetation state to global arrays.
  !--------------------------------------------------------------------!
  call veg_save_cell
  !--------------------------------------------------------------------!
  ! Save to global arrays.
  !--------------------------------------------------------------------!
  wbare  (1:nsoil_layers_max,i0,j0) = w  (1:nsoil_layers_max,1) ! m
  wvege  (0:nsoil_layers_max,i0,j0) = w  (0:nsoil_layers_max,2)
  htbare (0:nsoil_layers_max,i0,j0) = ht (0:nsoil_layers_max,1)
  htvege (0:nsoil_layers_max,i0,j0) = ht (0:nsoil_layers_max,2)
  snowbv (1:2,i0,j0) = snowd (1:2)
  betad_ij (i0,j0) = betad
  !--------------------------------------------------------------------!
  wdiffb = 0.0d0
  wdiffv = 0.0d0
  ediffb = 0.0d0
  ediffv = 0.0d0
  if (process_bare) then
    waterb_new = sum (wbare (1:nsoil_layers_max,i0,j0)) + trnf (1) + &
    &            trnff (1) + aevapb
    energyb_new = sum (htbare (1:nsoil_layers_max,i0,j0)) + &
    &             elh * aevapb + &
    &             trnffe (1) + &
    &             trnfe (1) + &
    &             athrm_tot (1) + asnsh_tot (1)
    wdiffb = waterb_new-waterb_old
    ediffb = energyb_new-energyb_old
  endif
  if (process_vege) then
    waterv_new = sum (wvege (0:nsoil_layers_max,i0,j0)) + trnf (2) + &
    &            trnff (2) + aevapw + aevapd + aevapvg
    energyv_new = sum (htvege (0:nsoil_layers_max,i0,j0)) + &
    &             elh * (aevapw + aevapd + aevapvg) + &
    &             trnffe (2) + &
    &             trnfe (2) + &
    &             athrm_tot (2) + asnsh_tot (2)
    wdiffv = waterv_new-waterv_old
    ediffv = energyv_new-energyv_old
  endif
  if (abs (wdiffb) > 1.0d-10) then
    write (*,*) 'wdiffb =',wdiffb,i0,j0
    write (98,*) 'wdiffb =',wdiffb,i0,j0
    write (98,*) trnf(1)
    write (98,*) trnff(1),process_bare
    write (98,*) fb,fv
    write (98,*) pr*dtsrc/float(nisurf)
    stop
  endif
  if (abs (ediffb) > 1.0d-5) then
    write (*,*) 'ediffb =',ediffb,i0,j0
    write (98,*) 'ediffb =',ediffb,i0,j0
    write (*,*) 'energyb_old',energyb_old,energyb_new
    write (*,*) 'flux imbalance    ',ediffb/(dtsrc/float(nisurf))
    write (*,*) 'LE heat flux      ',elh*aevapb/(dtsrc/float(nisurf))
    write (*,*) 'sub runoff heat   ',trnffe(1)/(dtsrc/float(nisurf))
    write (*,*) 'top runoff heat   ',trnfe(1)/(dtsrc/float(nisurf))
    write (*,*) 'sensible heat flux',asnsh_tot(1)/(dtsrc/float(nisurf))
    write (*,*) 'LWu heat flux     ',athrm_tot(1)/(dtsrc/float(nisurf))
    write (*,*) 'LWd heat flux     ',trht
    write (*,*) 'SWd heat flux     ',srht
    stop
  endif
  if (abs (wdiffv) > 1.0d-10) then
    write (*,*) 'wdiffv =',wdiffv,i0,j0
    write (98,*) 'wdiffv =',wdiffv,i0,j0
    write (98,*) 'aevapw ',aevapw
    write (98,*) trnf(2)
    write (98,*) trnff(2),process_vege
    write (98,*) fb,fv
    write (98,*) pr*dtsrc/float(nisurf)
    stop
  endif
  if (abs (ediffv) > 1.0d-5) then
    write (*,*) 'ediffv =',ediffv,i0,j0
    write (98,*) 'ediffv =',ediffv,i0,j0
    write (*,*) 'energyv_old',energyv_old,energyv_new
    write (*,*) 'flux imbalance    ',ediffv/(dtsrc/float(nisurf))
    write (*,*) 'LE heat flux      ',elh*(aevapw + aevapd + &
    &            aevapvg)/(dtsrc/float(nisurf))
    write (*,*) 'sub runoff heat   ',trnffe(2)/(dtsrc/float(nisurf))
    write (*,*) 'top runoff heat   ',trnfe(2)/(dtsrc/float(nisurf))
    write (*,*) 'sensible heat flux',asnsh_tot(2)/(dtsrc/float(nisurf))
    write (*,*) 'LWu heat flux     ',athrm_tot(2)/(dtsrc/float(nisurf))
    write (*,*) 'LWd heat flux     ',trht,t(i0,j0)-tf
    write (*,*) 'SWd heat flux     ',srht,tp(0,2),tp(1,2)
    stop
  endif
  !--------------------------------------------------------------------!
  ! Accumulate monthly soil moistures (m^3/m^3).
  !--------------------------------------------------------------------!
  aij_wtop (i0,j0) = aij_wtop (i0,j0) + &
  & (fb * wbare (1,i0,j0) + fv * wvege (1,i0,j0)) / dz_ij (i0,j0,1)
  !--------------------------------------------------------------------!
  ! Accumulate monthly water saturation (fraction).
  !--------------------------------------------------------------------!
  do nsoil_layer = 1, nsoil_layers
    aij_wtot (i0,j0) = aij_wtot (i0,j0) + &
    &                  (fb * soilsat (nsoil_layer,1) + &
    &                  fv * soilsat (nsoil_layer,2)) * &
    &                  dz (nsoil_layer) / (-zb (nsoil_layers+1))
  enddo
  !--------------------------------------------------------------------!
  ! Accumulate monthly betad (fraction).
  !--------------------------------------------------------------------!
  aij_betad (i0,j0) = aij_betad (i0,j0) + betad
  !--------------------------------------------------------------------!
  ! Accumulate monthly GPP (kgC/m2/mo).
  !--------------------------------------------------------------------!
  aij_gpp (i0,j0) = aij_gpp (i0,j0) + agpp
  !--------------------------------------------------------------------!
  ! Accumulate monthly NPP (kgC/m^2/mo).
  !--------------------------------------------------------------------!
  aij_npp (i0,j0) = aij_npp (i0,j0) + anpp
  !--------------------------------------------------------------------!
  ! Accumulate monthly NEE (kgC/m^2/mo).
  !--------------------------------------------------------------------!
  aij_nee (i0,j0) = aij_nee (i0,j0) + anee
  !--------------------------------------------------------------------!
  ! Accumulate monthly carbon pool masses (kgC/m^2).
  !--------------------------------------------------------------------!
  aij_cm (i0,j0,:) = aij_cm (i0,j0,:) + fv * cland (:)
  !--------------------------------------------------------------------!
  ! Accumulate monthly plant LAI (m^2/m^2).
  !--------------------------------------------------------------------!
  aij_lai (i0,j0,:) = aij_lai (i0,j0,:) + alait (:)
  !--------------------------------------------------------------------!
  ! Accumulate monthly plant fPAR (fraction).
  !--------------------------------------------------------------------!
  aij_fpar (i0,j0,:) = aij_fpar (i0,j0,:) + fpart (:)
  !--------------------------------------------------------------------!
  ! Accumulated monthly soil evaporation and sublimation (kg/m^2/mo.).
  !--------------------------------------------------------------------!
  aij_esoil (i0,j0) = aij_esoil (i0,j0) + rhow * (fb * aevapb + &
  &                                               fv * aevapvg)
  !--------------------------------------------------------------------!
  ! Accumulate monthly canopy interception (kg/m^2/s).
  !--------------------------------------------------------------------!
  aij_intercep (i0,j0) = aij_intercep (i0,j0) + rhow * fv * aevapw
  !--------------------------------------------------------------------!
  ! Accumulate monthly transpiration (kg/m^2/s).
  !--------------------------------------------------------------------!
  aij_transp (i0,j0) = aij_transp (i0,j0) + rhow * fv * aevapd
  !--------------------------------------------------------------------!
  ! Accumulate monthly total runoff (kg/m^2/s).
  !--------------------------------------------------------------------!
  aij_trunoff (i0,j0) = aij_trunoff (i0,j0) + rhow * fb * (trnf (1) + &
  &                     trnff (1))          + rhow * fv * (trnf (2) + &
  &                     trnff (2))
  !--------------------------------------------------------------------!
  ! Accumulate monthly surface runoff (kg/m^2/s).
  !--------------------------------------------------------------------!
  aij_srunoff (i0,j0) = aij_srunoff (i0,j0) + rhow * fb * trnf (1) + & 
  &                                           rhow * fv * trnf (2)
  !--------------------------------------------------------------------!
endif ! if (pearth > 0.0) then
!----------------------------------------------------------------------!
end subroutine proc_surf
!======================================================================!
