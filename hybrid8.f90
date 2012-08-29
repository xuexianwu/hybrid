!======================================================================!
program hybrid8 ! same position as 'earth' in modelE
!----------------------------------------------------------------------!
! Attempt to create full land surface model extracted from modelE.
! Use structure from modelE until have functioning soil scheme.
! Very carefully and slowly build up. Develop as basis for new
! Hybrid based on better physics. Keep options for slow parts to see if
! they matter.
!
! Then, put in growth model and dynamics and look at questions:
!
!    (i) increasing turnover in future forests with greater NPP;
!   (ii) importance of acclimation;
!  (iii) future food production;
!   (iv) historical NPP/NEE trends;
!    (v) grass phenology right;
!   (vi) constrain alpha, beta, gamma
!    (v) develop modelE
!   (vi) 1998 land carbon sink shift?
!  (vii) divergence problem
! (viii) maxval and phenology
!
!----------------------------------------------------------------------!
use control_parameters
use control_definitions_wg
use constants
use global_variables
use soil_variables
use veg_variables
use atmos_variables
use gen_atmos_wg
use physical_parameters_wg
use wgen_params_wg
!----------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------!
integer :: i,ii,jj,kk,ns,iostat,ioerr,nsoil_layer_save,iyr,iday,nyrm
integer :: yr_ooaf,co2_ooaf,icpool
integer :: yr_co2  ! CO2 year                                (Julian yr)
integer :: yr_rcp  ! RCP year                                (Julian yr)
integer :: jdpery_max ! Limit for checking climate input           (day)
real*8 :: lonw,latw,a_lon,b_lon,a_lat,b_lat
real*8 :: sla0,slre,slim,phase,dif
real*8 :: scs0,scsre,scsim
real*8 :: sfv,svh
real*8 :: dtsurf
real*8 :: pearth
real*8 :: z,frup,frdn
real*8 :: ptop,psf,plbot(2),psfmpt,sige(2),dsig
real*8 :: ratio ! Daily actual/potential radiation (ratio)
real*8 :: dtemp ! 24-hr temperature amplitude (degC)
real*8 :: raint,tc,qsat_modele,cloud,tak_ld,emm,tmean,qsat
!----------------------------------------------------------------------!
integer :: kyr_ahp,kday_ahp
real*8 :: dtmin_ahp,dtmax_ahp,rain_ahp
real*8 :: vap_ahp,rad_ahp,rc_ahp,windsp_ahp,pressure_ahp
!----------------------------------------------------------------------!
integer, parameter :: nfiles_in = 10 ! No. climate input files       (n)
!----------------------------------------------------------------------!
! Read run control parameters.
!----------------------------------------------------------------------!
open (10,file='driver.txt',status='old')
read (10,*) dtsrc   ! Source time step (s) = 1 ITU
read (10,*) houri   ! Start hour of model run (hr)
read (10,*) datei   ! Start date of model run (day)
read (10,*) monthi  ! Start month of model run (month)
read (10,*) yeari   ! Start year of model run (Julian yr)
read (10,*) restart ! Restart from restart file? (log)
read (10,*) sub_rsf ! Sub-directory for restart files (strg)
read (10,*) sub_ouf ! Sub-directory for output files (strg)
read (10,*) ite1    ! Length of model run from 1/1/year1 (hr)
read (10,*) lonw    ! Longitude for detailed diagnostics (deg)
read (10,*) latw    ! Latitude for detailed diagnostics (deg)
read (10,*) subd    ! Create sub-daily weather (log)
read (10,*) sing    ! Just one site? (log)
read (10,*) im      ! No. longitudinal grid boxes (n)
read (10,*) jm      ! No. latitudinal grid boxes (n)
read (10,*) daily   ! Direct input of daily climate? (log)
read (10,*) sclim   ! Single-site climate input file (log)
read (10,*) isi     ! isi-mip run? (log)
read (10,*) rcp     ! RCP to use for CO2 mixing ratio (strg)
read (10,*) chkclm  ! Read all climate to check? (log)
read (10,*) ooaf    ! OOAF simulation? (log)
read (10,*) physiol ! Full canopy physiology? (log)
close (10)


!----------------------------------------------------------------------!
! Locations for detailed diagnostics.
!----------------------------------------------------------------------!
open (88,file='OUT/'//sub_ouf//'/fort.88',status='unknown')
open (89,file='OUT/'//sub_ouf//'/fort.89',status='unknown')
open (90,file='OUT/'//sub_ouf//'/fort.90',status='unknown')
open (98,file='OUT/'//sub_ouf//'/fort.98',status='unknown')
open (99,file='OUT/'//sub_ouf//'/fort.99',status='unknown')
write (90,*) ' Starting...'
!----------------------------------------------------------------------!
b_lon = float (im - 1) / (360.0 - 360.0 / float (im))
a_lon = 1.0 - b_lon * (-180.0 + (360.0 / float (im)) / 2.0)
if ((im == 72) .and. (jm == 46)) then
  b_lat = float (jm - 3) / (180.0 / float (jm) - 176.0)
  a_lat = 2.0 - b_lat * (88.0 - (180.0 / float (jm)) / 2.0)
  if (latw > 88.0) then
    a_lat = 1.0
    b_lat = 0.0
  endif
  if (latw < -88.0) then
    a_lat = 46.0
    b_lat = 0.0
  endif
else
  b_lat = float (jm - 1) / (180.0 / float (jm) - 180.0)
  a_lat = 1.0 - b_lat * (90.0 - (180.0 / float (jm)) / 2.0)
endif
olon = nint (b_lon * lonw + a_lon)
olat = nint (b_lat * latw + a_lat)
!write (*,*)(float(513)-a_lon)/b_lon
!write (*,*)(float( 99)-a_lat)/b_lat
!----------------------------------------------------------------------!
! Open daily climate data file if use.
!----------------------------------------------------------------------!
if (sclim) then
  open (21,file='make_weather/daily_generated.txt',status='old')
endif
!----------------------------------------------------------------------!
! No. of times source executed per day (n)
!----------------------------------------------------------------------!
nday = 2 * nint (0.5 * sday / dtsrc) ! Execute nday times per day (n)
!---------------------------------------------------------------------!
ntpoints = nday ! Replaces direct generator control.
!----------------------------------------------------------------------!
! Set sizes of allocatable arrays.
!----------------------------------------------------------------------!
! atmos_variables
!----------------------------------------------------------------------!
allocate (raini_ij    (im,jm)) ! mm/d
allocate (tmin_ij     (im,jm)) ! degC
allocate (tmax_ij     (im,jm)) ! degC
allocate (vap_ij      (im,jm)) ! hPa or g/kg > kg/kg (isi)
allocate (radi_ij     (im,jm)) ! W/m^2
allocate (rc_ij       (im,jm)) ! MJ/m^2/d
allocate (windsp_ij   (im,jm)) ! m/s
allocate (pressure_ij (im,jm)) ! MPa or kPa (isi)
allocate (ldi_ij      (im,jm)) ! W/m^2
allocate (prec        (im,jm)) ! kg/m^2
allocate (pedn        (im,jm)) ! hPa
allocate (t           (im,jm)) ! K
allocate (fsf         (im,jm)) ! W/m^2
allocate (trhr        (im,jm)) ! W/m^2
allocate (q_surf      (im,jm)) ! kg/kg
allocate (am          (im,jm)) ! kg/m^2
allocate (sbeta_can   (im,jm)) ! ratio
allocate (ipar_df_can (im,jm)) ! umol(PAR)/m^2/s
allocate (ipar_dr_can (im,jm)) ! umol(PAR)/m^2/s
allocate (rain_ij     (im,jm,ntpoints))
allocate (pedn_ij     (im,jm,ntpoints))
allocate (tak_ij      (im,jm,ntpoints))
allocate (rad_ij      (im,jm,ntpoints))
allocate (ld_ij       (im,jm,ntpoints))
allocate (qsrf_ij     (im,jm,ntpoints))
allocate (sbeta_ij    (im,jm,ntpoints)) ! ratio
allocate (ipar_df_ij  (im,jm,ntpoints)) ! umol(PAR)/m^2/s
allocate (ipar_dr_ij  (im,jm,ntpoints)) ! umol(PAR)/m^2/s
!----------------------------------------------------------------------!
! soil_variables
!----------------------------------------------------------------------!
allocate (dz_ij  (im,jm,1:nsoil_layers_max)) ! m
allocate (q_ij   (im,jm,ntextures,nsoil_layers_max)) ! fraction
allocate (qk_ij  (im,jm,ntextures,nsoil_layers_max)) ! fraction
allocate (wbare  (1:nsoil_layers_max,im,jm)) ! m
allocate (wvege  (0:nsoil_layers_max,im,jm)) ! m
allocate (htbare (0:nsoil_layers_max,im,jm)) ! J/m^2
allocate (htvege (0:nsoil_layers_max,im,jm)) ! J/m^2
allocate (snowbv (nsurf_max,im,jm))   ! m (water eq.)
allocate (nsnow_layers_ij (nsurf_max,im,jm))   ! n
allocate (fr_snow_ij      (nsurf_max,im,jm))   ! fraction
allocate (wsn_ij (nsnow_layers_max,nsurf_max,im,jm)) ! m (water eq.)
allocate (aij_wtop  (im,jm)) ! Diag. monthly top soil water (m^3/m^3)
allocate (aij_wtot  (im,jm)) ! Diag. monthly total soil water (fraction)
allocate (aij_esoil (im,jm)) ! Diag. monthly soil evap (kg/m^2/s)
allocate (aij_intercep (im,jm)) ! Diag. monthly can int (kg/m^2/s)
allocate (aij_transp   (im,jm)) ! Diag. monthly transpiration (kg/m^2/s)
allocate (aij_trunoff  (im,jm)) ! Diag. monthly total runoff  (kg/m^2/s)
allocate (aij_srunoff  (im,jm)) ! Diag. mon. surface runoff   (kg/m^2/s)
allocate (wtot      (im,jm)) ! Mean grid box soil water (m)
allocate (sl_ij     (im,jm)) ! Slope (units?)
allocate (tearth    (im,jm)) ! Skin temperature (degC)
allocate (w  (0:nsoil_layers_max,nsurf_max)) ! Local soil water (m)
allocate (ht (0:nsoil_layers_max,nsurf_max)) ! Local soil heat (J/m^2)
allocate (dz (nsoil_layers_max))   ! Local soil layer depths (m)
allocate (zb (nsoil_layers_max+1)) ! Local soil layer boundaries (m)
allocate (zc (nsoil_layers_max+1)) ! Local soil layer centres (m)
allocate (snsh  (nsurf_max)) ! Sensible heat flux bare to air (J/s).
allocate (snshs (nsurf_max)) ! Sensible heat from bare sn to air (J/s)
allocate (hsn (nsnow_layers_max,nsurf_max)) ! Heat in snow (J?)
allocate (evap_tot (nsurf_max))
allocate (theta   (0:nsoil_layers_max+1,nsurf_max))
allocate (soilsat (0:nsoil_layers_max+1,nsurf_max))
!----------------------------------------------------------------------!
! Vegetation variables.
!----------------------------------------------------------------------!
allocate (cm       (im,jm,ncpools))
allocate (cland    (ncpools))
allocate (aij_gpp  (im,jm))
allocate (aij_npp  (im,jm))
allocate (aij_nee  (im,jm))
allocate (aij_cm   (im,jm,ncpools))
allocate (aij_lai  (im,jm,nplant))
allocate (aij_fpar (im,jm,nplant))
allocate (fpart    (nplant))
allocate (alaif    (nplant,im,jm))
allocate (alait    (nplant))
allocate (ala      (3,im,jm))
allocate (alaf     (3,nplant,im,jm))
allocate (nl       (ncanopy_layers_max))
allocate (n3       (ncanopy_layers_max))
allocate (fabsb    (ncanopy_layers_max))
!----------------------------------------------------------------------!
allocate (focean (im,jm))
allocate (flake0 (im,jm))
allocate (fland  (im,jm))
allocate (fearth (im,jm))
allocate (flice  (im,jm))
allocate (zatmo  (im,jm))
allocate (hlake  (im,jm))
allocate (f   (1:nsoil_layers_max+1,nsurf_max))
allocate (fh  (1:nsoil_layers_max+1,nsurf_max))
allocate (fch (0:1))
allocate (snowd (nsurf_max))
allocate (tsn1  (nsurf_max))
allocate (wsn (nsnow_layers_max,nsurf_max))
allocate (fr_snow (nsurf_max))
allocate (tp (0:nsoil_layers_max,nsurf_max))
allocate (q (ntextures,nsoil_layers_max))
allocate (qk (ntextures,nsoil_layers_max))
allocate (shc (0:nsoil_layers_max+1,nsurf_max))
allocate (fice (0:nsoil_layers_max,nsurf_max))
allocate (thets (0:nsoil_layers_max+1,nsurf_max))
allocate (thetm (0:nsoil_layers_max+1,nsurf_max))
allocate (ws (0:nsoil_layers_max,nsurf_max))
allocate (h (0:nsoil_layers_max+1,nsurf_max))
allocate (xinfc (nsurf_max))
allocate (d (0:nsoil_layers_max+1,nsurf_max))
allocate (xkus (nsoil_layers_max+1,nsurf_max))
allocate (xk (0:nsoil_layers_max+1,nsurf_max))
allocate (xku (nsoil_layers_max+1,nsurf_max))
allocate (xkh (nsoil_layers_max+1,nsurf_max))
allocate (xkhm(nsoil_layers_max+1,nsurf_max))
allocate (nsnow_layers (nsurf_max))
allocate (rnff (nsoil_layers_max+1,nsurf_max))
allocate (rnf         (nsurf_max))
allocate (dripw       (nsurf_max))
allocate (drips       (nsurf_max))
allocate (trnf        (nsurf_max))
allocate (trnff       (nsurf_max))
allocate (trnfe       (nsurf_max))
allocate (trnffe      (nsurf_max))
allocate (fhsng       (nsurf_max))
allocate (fhsng_scale (nsurf_max))
allocate (thrmsn      (nsurf_max))
allocate (htdripw     (nsurf_max))
allocate (htdrips     (nsurf_max))
allocate (thrm_tot    (nsurf_max))
allocate (snsh_tot    (nsurf_max))
allocate (athrm_tot   (nsurf_max))
allocate (asnsh_tot   (nsurf_max))
allocate (vdata (im,jm,nveg+1))
allocate (vdata_in(im,jm,nveg+1))
allocate (afb (im,jm))
allocate (avh (im,jm))
allocate (acs (3,im,jm))
allocate (afr (nsoil_layers_max,im,jm))
nyrm = ceiling(float(ite1)/8760.0)
allocate (lat (im,jm))
allocate (lon (im,jm))
allocate (veg_ci_ij(im,jm))
allocate (aij_betad (im,jm))
allocate (betad_ij (im,jm))
allocate (qfol_ij  (im,jm))
allocate (fc (nsurf_max))
allocate (evapdl (nsoil_layers_max,nsurf_max))
allocate (betadl (nsoil_layers_max))
allocate (f_swpl (nsoil_layers_max))
allocate (fr (nsoil_layers_max))
allocate (evap_max     (nsurf_max))
allocate (evap_max_wet (nsurf_max))
allocate (evap_max_dry (nsurf_max))
allocate (rg (ntpoints))
allocate (ld (ntpoints))
allocate (ta (ntpoints))
allocate (ppt(ntpoints))
allocate (p  (ntpoints))
allocate (tak(ntpoints))
allocate (qa (ntpoints))
allocate (u  (ntpoints))
allocate (sbeta(ntpoints))
allocate (ipar_dr(ntpoints))
allocate (ipar_df(ntpoints))
!---------------------------------------------------------------------!
! Open climate files.
!---------------------------------------------------------------------!
if ((.not. (sclim)) .and. (.not. (isi))) then
  if ((im == 720) .and. (jm == 360)) then
    open (10,file='make_weather/half/tmin.txt'    ,status='old')
    open (11,file='make_weather/half/tmax.txt'    ,status='old')
    open (12,file='make_weather/half/rain.txt'    ,status='old')
    open (13,file='make_weather/half/vap.txt'     ,status='old')
    open (14,file='make_weather/half/rad.txt'     ,status='old')
    open (15,file='make_weather/half/rc.txt'      ,status='old')
    open (16,file='make_weather/half/windsp.txt'  ,status='old')
    open (17,file='make_weather/half/pressure.txt',status='old')
    open (18,file='make_weather/half/ld.txt'      ,status='old')
    !open (18,file='make_weather/co2.txt'     ,status='old') ! ppm
  endif
  if ((im == 360) .and. (jm == 180) .and. (.not. (ooaf))) then
    open (10,file='make_weather/one/tmin.txt'    ,status='old')
    open (11,file='make_weather/one/tmax.txt'    ,status='old')
    open (12,file='make_weather/one/rain.txt'    ,status='old')
    open (13,file='make_weather/one/vap.txt'     ,status='old')
    open (14,file='make_weather/one/rad.txt'     ,status='old')
    open (15,file='make_weather/one/rc.txt'      ,status='old')
    open (16,file='make_weather/one/windsp.txt'  ,status='old')
    open (17,file='make_weather/one/pressure.txt',status='old')
    open (18,file='make_weather/one/ld.txt'      ,status='old')
    !open (18,file='make_weather/co2.txt'     ,status='old') ! ppm
  endif
  if ((im == 360) .and. (jm == 180) .and. (ooaf)) then
    open (10,file='make_weather/60ka/tmin.txt'    ,status='old')
    open (11,file='make_weather/60ka/tmax.txt'    ,status='old')
    open (12,file='make_weather/60ka/rain.txt'    ,status='old')
    open (13,file='make_weather/60ka/vap.txt'     ,status='old')
    open (14,file='make_weather/60ka/rad.txt'     ,status='old')
    open (15,file='make_weather/60ka/rc.txt'      ,status='old')
    open (16,file='make_weather/60ka/windsp.txt'  ,status='old')
    open (17,file='make_weather/60ka/pressure.txt',status='old')
    open (18,file='make_weather/60ka/ld.txt'      ,status='old')
  endif
  if (rcp == '8p5') then
    open (19,file='isi/RCP85_MIDYEAR_CONCENTRATIONS.DAT',status='old')
    do i = 1, 39
      read (19,*)
    enddo
    do kyr = 1765, yeari-1
      read (19,*) yr_rcp
    enddo
  endif
!---------------------------------------------------------------------!
! Read global climate fields.
!---------------------------------------------------------------------!
  if (chkclm) then
    jdpery_max = jdpery
  else
    jdpery_max = 1
  endif
  do kyr = 1, 1
!---------------------------------------------------------------------!
    do kday = 1, jdpery_max
      write (*,*) 'Reading kyr = ',kyr,'kday = ',kday
      do j0 = 1, jm
        read (10,*) (tmin_ij    (i0,j0),i0=1,im) ! degC
        read (11,*) (tmax_ij    (i0,j0),i0=1,im) ! degC
        read (12,*) (raini_ij   (i0,j0),i0=1,im) ! mm/d
        read (13,*) (vap_ij     (i0,j0),i0=1,im) ! hPa
        read (14,*) (radi_ij    (i0,j0),i0=1,im) ! W/m^2
        read (15,*) (rc_ij      (i0,j0),i0=1,im) ! MJ/m^2/d
        read (16,*) (windsp_ij  (i0,j0),i0=1,im) ! m/s
        read (17,*) (pressure_ij(i0,j0),i0=1,im) ! MPa
        read (18,*) (ldi_ij     (i0,j0),i0=1,im) ! W/m^2
      enddo
    enddo
    read (19,*) yr_rcp,co2eq_ppm,kyoto_co2eq_ppm,co2_ppm
  enddo
!---------------------------------------------------------------------!
! Close climate files.
!---------------------------------------------------------------------!
  do i = 1, nfiles_in
    close (i+9)
  enddo
endif
!----------------------------------------------------------------------!
! See if using isi-mip climate.
!----------------------------------------------------------------------!
if (isi) then
  open (15,file='make_weather/rc_yr.txt'    ,status='old') ! MJ/m2/d
  if (yeari == 1951) then
    ! degC
    open (10,file='isi/historical/tasmin1951-1960.txt',status='old')
    ! degC
    open (11,file='isi/historical/tasmax1951-1960.txt',status='old')
    ! mm/d
    open (12,file='isi/historical/pr1951-1960.txt'    ,status='old')
    ! g/kg
    open (13,file='isi/historical/huss1951-1960.txt'  ,status='old')
    ! W/m^2
    open (14,file='isi/historical/rsds1951-1960.txt'  ,status='old')
    open (16,file='make_weather/windsp.txt',status='old') ! m/s
    ! kPa
    open (17,file='isi/historical/psl1951-1960.txt'   ,status='old')
    ! W/m^2
    open (18,file='isi/historical/rlds1951-1960.txt'  ,status='old')
  elseif (yeari == 1991) then
  open (10,file='isi/hadgem/tasmin_hadgem2-es_rcp8p5_1991-2010.txt', &
  &  status='old') ! degC
  open (11,file='isi/hadgem/tasmax_hadgem2-es_rcp8p5_1991-2010.txt',&
  &  status='old') ! degC
  open (12,file='isi/hadgem/pr_hadgem2-es_rcp8p5_1991-2010.txt'    ,  &
  &  status='old') ! mm/d
  open (13,file='isi/hadgem/rhs_hadgem2-es_rcp8p5_1991-2010.txt'  ,  &
  &  status='old') ! %age
  open (14,file='isi/hadgem/rsds_hadgem2-es_rcp8p5_1991-2010.txt'  ,  &
  &  status='old') ! W/m^2
  open (16,file='make_weather/windsp.txt',status='old') ! m/s
  open (17,file='isi/hadgem/ps_hadgem2-es_rcp8p5_1991-2010.txt'   ,  &
  &  status='old') ! kPa
  open (18,file='isi/hadgem/rlds_hadgem2-es_rcp8p5_1991-2010.txt'  ,  &
  &  status='old') ! W/m^2
  elseif (yeari == 2006) then
    open (10,file='isi/tasmin_ipsl-cm5a-lr_2006-2010.txt',  &
    &  status='old') ! degC
    open (11,file='isi/tasmax_ipsl-cm5a-lr_2006-2010.txt',  &
    &  status='old') ! degC
    open (12,file='isi/pr_ipsl-cm5a-lr_2006-2010.txt'    ,  &
    &  status='old') ! mm/d
    open (13,file='isi/rhs_ipsl-cm5a-lr_2006-2010.txt'  ,  &
    &  status='old') ! %age
    open (14,file='isi/rsds_ipsl-cm5a-lr_2006-2010.txt'  ,  &
    &  status='old') ! W/m^2
    open (16,file='make_weather/windsp.txt',status='old') ! m/s
    open (17,file='isi/ps_ipsl-cm5a-lr_2006-2010.txt'   ,  &
    &  status='old') ! kPa
    open (18,file='isi/rlds_ipsl-cm5a-lr_2006-2010.txt'  ,  &
    &  status='old') ! W/m^2
  elseif (yeari == 2011) then
    open (10,file='isi/tasmin2011-2020.txt',status='old') ! degC
    open (11,file='isi/tasmax2011-2020.txt',status='old') ! degC
    open (12,file='isi/pr2011-2020.txt'    ,status='old') ! mm/d
    open (13,file='isi/huss2011-2020.txt'  ,status='old') ! g/kg
    open (14,file='isi/rsds2011-2020.txt'  ,status='old') ! W/m^2
    open (16,file='make_weather/windsp.txt',status='old') ! m/s
    open (17,file='isi/psl2011-2020.txt'   ,status='old') ! kPa
    open (18,file='isi/rlds2011-2020.txt'  ,status='old') ! W/m^2
  elseif (yeari == 2081) then
  open (10,file='isi/hadgem/tasmin_hadgem2-es_rcp8p5_2081-2099.txt', &
  &  status='old') ! degC
  open (11,file='isi/hadgem/tasmax_hadgem2-es_rcp8p5_2081-2099.txt',&
  &  status='old') ! degC
  open (12,file='isi/hadgem/pr_hadgem2-es_rcp8p5_2081-2099.txt'    ,  &
  &  status='old') ! mm/d
  open (13,file='isi/hadgem/rhs_hadgem2-es_rcp8p5_2081-2099.txt'  ,  &
  &  status='old') ! %age
  open (14,file='isi/hadgem/rsds_hadgem2-es_rcp8p5_2081-2099.txt'  ,  &
  &  status='old') ! W/m^2
  open (16,file='make_weather/windsp.txt',status='old') ! m/s
  open (17,file='isi/hadgem/ps_hadgem2-es_rcp8p5_2081-2099.txt'   ,  &
  &  status='old') ! kPa
  open (18,file='isi/hadgem/rlds_hadgem2-es_rcp8p5_2081-2099.txt'  ,  &
  &  status='old') ! W/m^2
  elseif (yeari == 2091) then
    open (10,file='isi/tasmin_ipsl-cm5a-lr_2091-2099.txt',  &
    &  status='old') ! degC
    open (11,file='isi/tasmax_ipsl-cm5a-lr_2091-2099.txt',  &
    &  status='old') ! degC
    open (12,file='isi/pr_ipsl-cm5a-lr_2091-2099.txt'    ,  &
    &  status='old') ! mm/d
    open (13,file='isi/rhs_ipsl-cm5a-lr_2091-2099.txt'  ,  &
    &  status='old') ! %age
    open (14,file='isi/rsds_ipsl-cm5a-lr_2091-2099.txt'  ,  &
    &  status='old') ! W/m^2
    open (16,file='make_weather/windsp.txt',status='old') ! m/s
    open (17,file='isi/ps_ipsl-cm5a-lr_2091-2099.txt'   ,  &
    &  status='old') ! kPa
    open (18,file='isi/rlds_ipsl-cm5a-lr_2091-2099.txt'  ,  &
    &  status='old') ! W/m^2
  else
    write (*,*) 'No isi climate for yeari ',yeari
    stop
  endif
  if (rcp == '8p5') then
    open (19,file='isi/RCP85_MIDYEAR_CONCENTRATIONS.DAT',status='old')
    do i = 1, 39
      read (19,*)
    enddo
    do kyr = 1765, yeari-1
      read (19,*) yr_rcp
    enddo
  endif
  if (chkclm) then
    jdpery_max = jdpery
  else
    jdpery_max = 1
  endif
  do kday = 1, jdpery_max
    read (15,*) (rc_ij (i0,j0),i0=1,im) ! MJ/m^2/d
  enddo
  do kyr = yeari, yeari+nyrm-1
    do kday = 1, jdpery_max
      write (*,*) 'Reading kyr = ',kyr,'kday = ',kday
      do j0 = 1, jm
        read (10,*) (tmin_ij     (i0,j0),i0=1,im) ! degC
        read (11,*) (tmax_ij     (i0,j0),i0=1,im) ! degC
        read (12,*) (raini_ij    (i0,j0),i0=1,im) ! mm/d
        read (13,*) (vap_ij      (i0,j0),i0=1,im) ! hPa or g/kg (isi)
        read (14,*) (radi_ij     (i0,j0),i0=1,im) ! W/m^2
        read (16,*) (windsp_ij   (i0,j0),i0=1,im) ! m/s
        read (17,*) (pressure_ij (i0,j0),i0=1,im) ! kPa
        read (18,*) (ldi_ij      (i0,j0),i0=1,im) ! W/m^2
      enddo
    enddo
      read (19,*) yr_rcp,co2eq_ppm,kyoto_co2eq_ppm,co2_ppm
  enddo
  do i = 1, nfiles_in
    close (i+9)
  enddo
  write (*,*) 'isi climate input reads OK'
endif
!----------------------------------------------------------------------!
! Surface time step (s).
!----------------------------------------------------------------------!
dtsurf = dtsrc / float (nisurf)
!----------------------------------------------------------------------!
! Time step for ground hydrology (s).
!----------------------------------------------------------------------!
dt = dtsurf
!----------------------------------------------------------------------!
! Open relevant input file for global surface fields.
!----------------------------------------------------------------------!
if ((im == 72) .and. (jm == 46)) then
  open (10,file='make_fields/surface_fields_72x46.txt',&
  &status='old',iostat=ioerr , err=100)
elseif ((im == 360) .and. (jm == 180)) then
  !open (10,file='make_fields/surface_fields_360x180.txt',&
  !&status='old',iostat=ioerr, err=100)
  open (10,file='make_fields/surface_fields_360x180_h_s.txt',&
  &status='old',iostat=ioerr, err=100)
elseif ((im == 720) .and. (jm == 360)) then
  open (10,file='make_fields/surface_fields_720x360_h_s.txt',&
  &status='old',iostat=ioerr, err=100)
else
  write (90,*) 'No input fields available for grid ',im,jm
  write (90,*) 'Simulation terminated'
  stop
endif
!----------------------------------------------------------------------!
! Read all global input fields.
! Read surface fractions (no LI).
!----------------------------------------------------------------------!
do j0 = 1, jm
  read (10,*) (focean (i0,j0), i0 = 1, im)
enddo
!----------------------------------------------------------------------!
do j0 = 1, jm
  read (10,*) (flake0 (i0,j0), i0 = 1, im)
enddo
!----------------------------------------------------------------------!
do j0 = 1, jm
  read (10,*) (fearth (i0,j0), i0 = 1, im)
enddo
!----------------------------------------------------------------------!
do j0 = 1, jm
  read (10,*) (flice (i0,j0), i0 = 1, im)
enddo
!----------------------------------------------------------------------!
do j0 = 1, jm
  read (10,*) (zatmo (i0,j0), i0 = 1, im)
enddo
zatmo (:,:) = zatmo (:,:) * grav ! Geopotential
!----------------------------------------------------------------------!
do j0 = 1, jm
  read (10,*) (hlake (i0,j0), i0 = 1, im)
enddo
!----------------------------------------------------------------------!
do nsoil_layer = 1, nsoil_layers_max
  do j0 = 1, jm
    read (10,*) (dz_ij (i0,j0,nsoil_layer), i0 = 1, im)
  enddo
enddo
!----------------------------------------------------------------------!
do nsoil_layer = 1, nsoil_layers_max
  do ntexture = 1, ntextures
    do j0 = 1, jm
      read (10,*) (q_ij (i0,j0,ntexture,nsoil_layer), i0 = 1, im)
    enddo
  enddo
enddo
!----------------------------------------------------------------------!
do nsoil_layer = 1, nsoil_layers_max
  do ntexture = 1, ntextures
    do j0 = 1, jm
      read (10,*) (qk_ij (i0,j0,ntexture,nsoil_layer), i0 = 1, im)
    enddo
  enddo
enddo
!*******************
qk_ij (:,:,:,:) = q_ij (:,:,:,:)
!----------------------------------------------------------------------!
do j0 = 1, jm
  read (10,*) (sl_ij (i0,j0), i0 = 1, im)
enddo
!----------------------------------------------------------------------!
! Read fraction of grid box of veg type 1-11.
!  1 = bright bare soil (in combination with 10)
!  2 = tundra
!  3 = grassland
!  4 = shrub/grassland
!  5 = tree/grassland
!  6 = deciduous forest
!  7 = evergreen forest
!  8 = rainforest
!  9 = cultivation
! 10 = dark bare soil (in combination with 1)
!----------------------------------------------------------------------!
do kveg = 1, nveg
  do j0 = 1, jm
    read (10,*) (vdata_in (i0,j0,kveg), i0 = 1, im)
  enddo
enddo
vdata (:,:,:) = float (vdata_in (:,:,:)) / 10000.0d0
!----------------------------------------------------------------------!
! Zero-out vdata (11) until it is properly read in.
!----------------------------------------------------------------------!
vdata (:,:,11) = 0.0d0
!----------------------------------------------------------------------!
close (10)
!----------------------------------------------------------------------!
! Set all vegetation C pools to zero.
!----------------------------------------------------------------------!
 cm (:,:,:) = zero
!----------------------------------------------------------------------!
! As in MODELE.f, make sure constraints are satisfied by defining FLAND/
! FEARTH as residual terms (deals with SP => DP problem).
!----------------------------------------------------------------------!
do j0 = 1, jm
  do i0 = 1, im
    !------------------------------------------------------------------!
    if (focean (i0,j0) > 0.0d0) then
      fland (i0,j0) = 1.0d0 - focean (i0,j0) ! Land fraction
      if (flake0 (i0,j0) > 0.0d0) then
        write (90,*) 'Ocean and lake cannot exist in same grid box', &
        &            i0,j0,focean(i0,j0),flake0(i0,j0)
        flake0 (i0,j0) = 0.0d0
      endif
    elseif (flake0 (i0,j0) > 0.0d0) then
      fland (i0,j0) = 1.0d0 - flake0 (i0,j0)
    else
      fland (i0,j0) = 1.0d0
    endif
    !------------------------------------------------------------------!
    ! Ensure no round-off error affects land with ice and earth.
    !------------------------------------------------------------------!
    if (((flice (i0,j0) - fland (i0,j0)) > -1.0d-4) .and. &
    &   (flice (i0,j0) > 0.0d0)) then
      flice (i0,j0) = fland (i0,j0)
      fearth (i0,j0) = 0.0d0
    else
      fearth (i0,j0) = fland (i0,j0) - flice (i0,j0)
    endif
    !------------------------------------------------------------------!
    ! Needed when using soils data from non-GISS source (e.g. HWSD).
    !------------------------------------------------------------------!
    if (dz_ij (i0,j0,1) == 0.0d0) fearth (i0,j0) = 0.0d0
    if (dz_ij (i0,j0,1) >  0.0d0) fearth (i0,j0) = 1.0d0
    !------------------------------------------------------------------!
    ! If subsoil has no texture, limit soil depth to topsoil.
    !------------------------------------------------------------------!
    if (sum(q_ij(i0,j0,1:ntextures,3:nsoil_layers_max)) == 0.0d0) then
      dz_ij (i0,j0,3:nsoil_layers_max) = 0.0d0
    endif
    !------------------------------------------------------------------!
    ! Catch boxes with no climate.
    !------------------------------------------------------------------!
    if (.not. (sclim)) then
      if (raini_ij (i0,j0) < -90.0d0) fearth (i0,j0) = 0.0d0
    endif
    !------------------------------------------------------------------!
    ! Needed with new (SiB) vegetation dataset.
    !-----------------------------------------------------------------!
    if (sum(vdata(i0,j0,1:10)) == 0.0) fearth(i0,j0) = 0.0d0
    !------------------------------------------------------------------!
  enddo
enddo
fland  (2:im, 1) = fland  (1, 1)
fland  (2:im,jm) = fland  (1,jm)
fearth (2:im, 1) = fearth (1, 1)
fearth (2:im,jm) = fearth (1,jm)
flice  (2:im, 1) = flice  (1, 1)
flice  (2:im,jm) = flice  (1,jm)
!----------------------------------------------------------------------!
! Save some of the forcings for analysis.
!----------------------------------------------------------------------!
open (20, file='OUT/'//sub_ouf//'/fearth_field.txt',status='unknown')
do j0 = 1, jm
  write (20,'(720f7.3)') (fearth (i0,j0), i0 = 1, im)
enddo
close (20)
!----------------------------------------------------------------------!
open (20, file='OUT/'//sub_ouf//'/dz_field.txt',status='unknown')
do j0 = 1, jm
  write (20,'(720f7.3)') (dz_ij (i0,j0,1), i0 = 1, im)
enddo
close (20)
!----------------------------------------------------------------------!
open (20, file='OUT/'//sub_ouf//'/q_field.txt',status='unknown')
do j0 = 1, jm
  write (20,'(720f7.3)') (q_ij (i0,j0,1,1), i0 = 1, im)
enddo
close (20)
!----------------------------------------------------------------------!
open (20, file='OUT/'//sub_ouf//'/qk_field.txt',status='unknown')
do j0 = 1, jm
  write (20,'(720f7.3)') (qk_ij (i0,j0,1,1), i0 = 1, im)
enddo
close (20)
!----------------------------------------------------------------------!
open (20, file='OUT/'//sub_ouf//'/sl_field.txt',status='unknown')
do j0 = 1, jm
  write (20,'(720f7.3)') (sl_ij (i0,j0), i0 = 1, im)
enddo
close (20)
!----------------------------------------------------------------------!
! Set some surface properties.
!----------------------------------------------------------------------!
ala  (:,:,:)   = 0.0d0
alaf (:,:,:,:) = 0.0d0
acs (:,:,:) = 0.0d0
avh (:,:)   = 0.0d0
afr (:,:,:) = 0.0d0
acs (1,:,:) = 0.01d0
do j0 = 1, jm
  do i0 = 1, im
    pearth = fearth (i0,j0)
    ! Bare fraction (fraction)
    afb (i0,j0) = vdata (i0,j0,1) + vdata (i0,j0,10)
    if (afb (i0,j0) > 0.999d0) afb (i0,j0) = 1.0d0
    if ((pearth <= 0.0d0) .or. (afb (i0,j0) >= 1.0d0)) cycle
    sfv = 0.0d0
    sla0 =0.0d0
    slre =0.0d0
    slim =0.0d0
    scs0  = 0.0d0
    scsre = 0.0d0
    scsim = 0.0d0
    svh = 0.0d0
    do kveg = 1, nveg - 2
      phase = twopi * laday (kveg) / 365.d0
      fv = vdata (i0,j0,kveg+1)
      sfv = sfv + fv
      svh = svh + fv * vhght (kveg)
      dif = (alamax (kveg) - alamin (kveg))
      sla0 = sla0 + fv * (alamax (kveg) + alamin (kveg))
      slre = slre + fv * dif * cos (phase)
      slim = slim + fv * dif * sin (phase)
      !----------------------------------------------------------------!
      alaf (1,kveg,i0,j0) = 0.5d0 * (alamax (kveg) + alamin (kveg))
      alaf (2,kveg,i0,j0) = 0.5d0 * dif * cos (phase)
      alaf (3,kveg,i0,j0) = 0.5d0 * dif * sin (phase)
      !----------------------------------------------------------------!
      scs0  = scs0  + fv * (alamax (kveg) + alamin (kveg)) / rsar (kveg)
      scsre = scsre + fv * dif * cos (phase) / rsar (kveg)
      scsim = scsim + fv * dif * sin (phase) / rsar (kveg)
    enddo
    !------------------------------------------------------------------!
    ala (1,i0,j0) = 0.5d0 / sfv * sla0
    ala (2,i0,j0) = 0.5d0 / sfv * slre
    ala (3,i0,j0) = 0.5d0 / sfv * slim
    !------------------------------------------------------------------!
    acs (1,i0,j0) = 0.5d0 / sfv * scs0
    acs (2,i0,j0) = 0.5d0 / sfv * scsre
    acs (3,i0,j0) = 0.5d0 / sfv * scsim
    !------------------------------------------------------------------!
    avh (i0,j0) = svh / sfv ! Vegetation height (m)
    !------------------------------------------------------------------!
    ! Compute root fraction afr averaged over vegetation types.
    !------------------------------------------------------------------!
    do nsoil_layer = 1, nsoil_layers_max
      dz (nsoil_layer) = dz_ij (i0,j0,nsoil_layer)
      if (dz (nsoil_layer) <= 0.0d0) goto 320
    enddo
 320 nsoil_layer_save = nsoil_layer - 1
    do kveg = 1, nveg - 2
      fv = vdata (i0,j0,kveg+1)
      z = 0.0d0
      frup = 0.0d0
      do nsoil_layer = 1, nsoil_layer_save
        z = z + dz (nsoil_layer)
        frdn = aroot (kveg) * z ** broot (kveg)
        frdn = min (frdn,one)
        if (nsoil_layer == nsoil_layer_save) frdn = 1.0d0
        afr (nsoil_layer,i0,j0) = afr (nsoil_layer,i0,j0) + &
        &                         fv * (frdn - frup)
        frup = frdn
      enddo
    end do
    do nsoil_layer = 1, nsoil_layer_save
      afr (nsoil_layer,i0,j0) = afr (nsoil_layer,i0,j0) / &
      &                         (1.0d0 - afb (i0,j0))
    enddo
    !------------------------------------------------------------------!
  enddo
enddo
if (sing) then
  do j0 = olat, olat
    do i0 = olon, olon
      if (fearth (i0,j0) == 1.0d0) fearth (i0,j0) = 2.0d0
    enddo
  enddo
  do j0 = 1, jm
    do i0 = 1, im
     if (fearth (i0,j0) .ne. 2.0d0) fearth (i0,j0) = 0.0d0
    enddo
  enddo
  do j0 = olat, olat
    do i0 = olon,olon
      if (fearth (i0,j0) == 2.0d0) fearth (i0,j0) = 1.0d0
    enddo
  enddo
endif
write (90,*) lonw,latw
write (90,*) olon,olat,fearth(olon,olat),afb(olon,olat)
write (90,*) 'fbright bare soil  = ',vdata(olon,olat, 1)
write (90,*) 'ftundra            = ',vdata(olon,olat, 2)
write (90,*) 'fgrassland         = ',vdata(olon,olat, 3)
write (90,*) 'fshrub/grassland   = ',vdata(olon,olat, 4)
write (90,*) 'ftree/grassland    = ',vdata(olon,olat, 5)
write (90,*) 'fdeciduous forest  = ',vdata(olon,olat, 6)
write (90,*) 'fevergreen forest  = ',vdata(olon,olat, 7)
write (90,*) 'frainforest        = ',vdata(olon,olat, 8)
write (90,*) 'fcultivation       = ',vdata(olon,olat, 9)
write (90,*) 'fdark bare soil    = ',vdata(olon,olat,10)
write (90,*)
!----------------------------------------------------------------------!
! Set up table of volumetric moisture contents and calculate
! soil moisture diffusivities and conductivities.
!----------------------------------------------------------------------!
call hl0
!----------------------------------------------------------------------!
! Initialise water on canopy and in soil layers (m)
!----------------------------------------------------------------------!
do i0 = 1, im
  do j0 = 1, jm
    !------------------------------------------------------------------!
    veg_ci_ij (i0,j0) = 0.0d0 ! Internal leaf CO2 (mol/m^3)
    qfol_ij   (i0,j0) = 0.0d0 ! Canopy surface mixing ratio (kg/kg)
    wbare   (:,i0,j0) = 0.0d0 ! Soil water (m)
    wvege   (:,i0,j0) = 0.0d0 ! Soil water (m)
    !------------------------------------------------------------------!
    ! Fraction of grid box with soil (n).
    !------------------------------------------------------------------!
    pearth = fearth (i0,j0)
    !------------------------------------------------------------------!
    if (pearth > 0.0) then
      !----------------------------------------------------------------!
      qfol_ij   (i0,j0) = 3.0d-6   ! Canopy surface mixing ratio (kg/kg)
      veg_ci_ij (i0,j0) = 0.0127d0 ! Internal leaf CO2 (mol/m^3)
      ! Skin temperature (deg C).
      tearth (i0,j0) = (tmin_ij (i0,j0) + tmax_ij (i0,j0)) / 2.0d0
      !----------------------------------------------------------------!
      ! To get saturated water content of each layer.
      !----------------------------------------------------------------!
      call ghinij
      do nsoil_layer = 1, nsoil_layers
        wbare (nsoil_layer,i0,j0) = 0.5d0 * ws (nsoil_layer,1) ! m
        wvege (nsoil_layer,i0,j0) = 0.5d0 * ws (nsoil_layer,2) ! m
      enddo
      call veg_set_cell
      !----------------------------------------------------------------!
      htvege (0,i0,j0) = tearth (i0,j0) * shc (0,2) ! J/m^2
      htbare (0,i0,j0) = htvege (0,i0,j0) ! not used (J/m^2)
      do nsoil_layer = 1, nsoil_layers
        htbare (nsoil_layer,i0,j0) = tearth (i0,j0) * &
        &                            shc (nsoil_layer,1)
        htvege (nsoil_layer,i0,j0) = tearth (i0,j0) * &
        &                            shc (nsoil_layer,2)
      enddo
    endif
  enddo
enddo
!----------------------------------------------------------------------!
! Initialise snow depths (m water equivalent)
!----------------------------------------------------------------------!
do i0 = 1, im
  do j0 = 1, jm
    snowbv (1,i0,j0) = 1.0d0
    snowbv (2,i0,j0) = 1.0d0
  enddo
enddo
!----------------------------------------------------------------------!
! Initialise numbers of snow layers (n).
!----------------------------------------------------------------------!
do i0 = 1, im
  do j0 = 1, jm
    nsnow_layers_ij (1,i0,j0) = 1
    nsnow_layers_ij (2,i0,j0) = 1
  enddo
enddo
!----------------------------------------------------------------------!
! Initialise snow depths in snow layers (m water equivalent).
!----------------------------------------------------------------------!
do i0 = 1, im
  do j0 = 1, jm
    wsn_ij (:,1,i0,j0) = 0.0d0
    wsn_ij (:,2,i0,j0) = 0.0d0
  enddo
enddo
!----------------------------------------------------------------------!
! Initialise fractions of snow on each surface (fraction).
!----------------------------------------------------------------------!
do i0 = 1, im
  do j0 = 1, jm
    fr_snow_ij (1,i0,j0) = 0.0d0
    fr_snow_ij (2,i0,j0) = 0.0d0
  enddo
enddo
!----------------------------------------------------------------------!
if (restart) then
  if (yeari > 0) then
    if (monthi > 1) then
      write (adate(1:7),'(a3,i4.4)') amon(monthi-1),yeari
    else
      write (adate(1:7),'(a3,i4.4)') amon(12),yeari-1
    endif
  else
    if (monthi > 1) then
      write (adate(1:7),'(a3,i4.4)') amon(monthi-1),-yeari/10
    else
      write (adate(1:7),'(a3,i4.4)') amon(12),-(yeari+1)/10
    endif
  endif
  !------------------------------------------------------------------!
  open (10,file='OUT/'//sub_rsf//'/'//adate(1:7)//'.rsf',status='unknown')
  write (90,*) 'Restarting from file ','OUT/'//sub_rsf//'/'//adate(1:7)//'.rsf'
  !------------------------------------------------------------------!
  do j0 = 1, jm
    read (10,*) (wvege           (0,i0,j0)  ,i0=1,im)
    read (10,*) (htvege          (0,i0,j0)  ,i0=1,im)
    read (10,*) (cm              (i0,j0,:)  ,i0=1,im)
    read (10,*) (veg_ci_ij       (i0,j0)    ,i0=1,im)
    read (10,*) (qfol_ij         (i0,j0)    ,i0=1,im)
    read (10,*) (snowbv          (1,i0,j0)  ,i0=1,im)
    read (10,*) (snowbv          (2,i0,j0)  ,i0=1,im)
    read (10,*) (nsnow_layers_ij (1,i0,j0)  ,i0=1,im)
    read (10,*) (nsnow_layers_ij (2,i0,j0)  ,i0=1,im)
    read (10,*) (wsn_ij          (:,1,i0,j0),i0=1,im)
    read (10,*) (wsn_ij          (:,2,i0,j0),i0=1,im)
    read (10,*) (fr_snow_ij      (1,i0,j0)  ,i0=1,im)
    read (10,*) (fr_snow_ij      (2,i0,j0)  ,i0=1,im)
  enddo
  do nsoil_layer = 1, nsoil_layers_max
    do j0 = 1, jm
      read (10,*) (wbare  (nsoil_layer,i0,j0),i0=1,im)
      read (10,*) (wvege  (nsoil_layer,i0,j0),i0=1,im)
      read (10,*) (htbare (nsoil_layer,i0,j0),i0=1,im)
      read (10,*) (htvege (nsoil_layer,i0,j0),i0=1,im)
    enddo
  enddo
  !------------------------------------------------------------------!
  ! Revert to internal heat units (J/m^2).
  !------------------------------------------------------------------!
  htvege (:,:,:) = 1.0d6 * htvege (:,:,:)
  htbare (:,:,:) = 1.0d6 * htbare (:,:,:)
  !------------------------------------------------------------------!
  close (10)
  !------------------------------------------------------------------!
endif
!----------------------------------------------------------------------!
! Maximum iterations to find hl in hydro.
!----------------------------------------------------------------------!
jcm = nint (log (float (nth)) / log (2.0d0))
!----------------------------------------------------------------------!
! Sigma coordinate is pressure/surface pressure. dsig is difference in
! sigma between bottom and top of bottom layer. Values taken from
! RES_M23.f, so presumably 23-layer version.
! Pressure in hPa.
!----------------------------------------------------------------------!
ptop = 150.0d0
psf = 984.0d0
plbot (1) = psf
plbot (2) = 960.0d0
psfmpt = psf - ptop
sige (1) = (plbot (1) - ptop) / psfmpt
sige (2) = (plbot (2) - ptop) / psfmpt
dsig = sige (1) - sige (2)
!----------------------------------------------------------------------!
! Following just for some simple diagnostics.
!----------------------------------------------------------------------!
do i0 = 1, im
  do j0 = 1, jm
    ! Soil layer depths.
    !write (89,'(2i5,6f7.3)') i0,j0,dz_ij(i0,j0,1:6)
    !write (89,'(2i5,5f7.3)') i0,j0,q_ij(i0,j0,1:5,6)
    write (89,'(2i5,10f7.3)') i0,j0,vdata(i0,j0,1:10)
  enddo
enddo
!----------------------------------------------------------------------!
iyear1 = yeari ! default
ihri = houri + hr_in_day * (datei - 1 + jdendofm (monthi - 1) + &
&      jdpery * (yeari - iyear1))
itimei= ihri * nday / hr_in_day
itime = itimei
itimee = ite1
!----------------------------------------------------------------------!
! Main loop.
!----------------------------------------------------------------------!
write (*,*) 'Starting main routine...'
!---------------------------------------------------------------------!
! Open daily weather files.
!---------------------------------------------------------------------!
if ((.not. (sclim)) .and. (.not. (isi))) then
  if ((im == 360) .and. (jm == 180) .and. (.not. (ooaf))) then
    open (10,file='make_weather/one/tmin.txt'    ,status='old')
    open (11,file='make_weather/one/tmax.txt'    ,status='old')
    open (12,file='make_weather/one/rain.txt'    ,status='old')
    open (13,file='make_weather/one/vap.txt'     ,status='old')
    open (14,file='make_weather/one/rad.txt'     ,status='old')
    open (15,file='make_weather/one/rc.txt'      ,status='old')
    open (16,file='make_weather/one/windsp.txt'  ,status='old')
    open (17,file='make_weather/one/pressure.txt',status='old')
    open (18,file='make_weather/one/ld.txt'      ,status='old')
  endif
  if ((im == 720) .and. (jm == 360)) then
    open (10,file='make_weather/half/tmin.txt'    ,status='old')
    open (11,file='make_weather/half/tmax.txt'    ,status='old')
    open (12,file='make_weather/half/rain.txt'    ,status='old')
    open (13,file='make_weather/half/vap.txt'     ,status='old')
    open (14,file='make_weather/half/rad.txt'     ,status='old')
    open (15,file='make_weather/half/rc.txt'      ,status='old')
    open (16,file='make_weather/half/windsp.txt'  ,status='old')
    open (17,file='make_weather/half/pressure.txt',status='old')
    open (18,file='make_weather/half/ld.txt'      ,status='old')
  endif
  if ((im == 360) .and. (jm == 180) .and. (ooaf)) then
    open (10,file='make_weather/60ka/tmin.txt'    ,status='old')
    open (11,file='make_weather/60ka/tmax.txt'    ,status='old')
    open (12,file='make_weather/60ka/rain.txt'    ,status='old')
    open (13,file='make_weather/60ka/vap.txt'     ,status='old')
    open (14,file='make_weather/60ka/rad.txt'     ,status='old')
    open (15,file='make_weather/60ka/rc.txt'      ,status='old')
    open (16,file='make_weather/60ka/windsp.txt'  ,status='old')
    open (17,file='make_weather/60ka/pressure.txt',status='old')
    open (18,file='make_weather/60ka/ld.txt'      ,status='old')
  endif
  if (rcp == '8p5') then
    open (19,file='isi/RCP85_MIDYEAR_CONCENTRATIONS.DAT',status='old')
    do i = 1, 39
      read (19,*)
    enddo
    do kyr = 1765, yeari-1
      read (19,*) yr_rcp
    enddo
  endif
endif
!---------------------------------------------------------------------!
if (isi) then
  open (15,file='make_weather/rc_yr.txt'    ,status='old') ! MJ/m2/d
  if (yeari == 1951) then
    ! degC
    open (10,file='isi/historical/tasmin1951-1960.txt',status='old')
    ! degC
    open (11,file='isi/historical/tasmax1951-1960.txt',status='old')
    ! mm/d
    open (12,file='isi/historical/pr1951-1960.txt'    ,status='old')
    ! g/kg
    open (13,file='isi/historical/huss1951-1960.txt'  ,status='old')
    ! W/m^2
    open (14,file='isi/historical/rsds1951-1960.txt'  ,status='old')
    open (16,file='make_weather/windsp.txt',status='old') ! m/s
    ! kPa
    open (17,file='isi/historical/psl1951-1960.txt'   ,status='old')
    ! W/m^2
    open (18,file='isi/historical/rlds1951-1960.txt'  ,status='old')
  elseif (yeari == 1991) then
  open (10,file='isi/hadgem/tasmin_hadgem2-es_rcp8p5_1991-2010.txt', &
  &  status='old') ! degC
  open (11,file='isi/hadgem/tasmax_hadgem2-es_rcp8p5_1991-2010.txt',&
  &  status='old') ! degC
  open (12,file='isi/hadgem/pr_hadgem2-es_rcp8p5_1991-2010.txt'    ,  &
  &  status='old') ! mm/d
  open (13,file='isi/hadgem/rhs_hadgem2-es_rcp8p5_1991-2010.txt'  ,  &
  &  status='old') ! %age
  open (14,file='isi/hadgem/rsds_hadgem2-es_rcp8p5_1991-2010.txt'  ,  &
  &  status='old') ! W/m^2
  open (16,file='make_weather/windsp.txt',status='old') ! m/s
  open (17,file='isi/hadgem/ps_hadgem2-es_rcp8p5_1991-2010.txt'   ,  &
  &  status='old') ! kPa
  open (18,file='isi/hadgem/rlds_hadgem2-es_rcp8p5_1991-2010.txt'  ,  &
  &  status='old') ! W/m^2
  elseif (yeari == 2006) then
    open (10,file='isi/tasmin_ipsl-cm5a-lr_2006-2010.txt',  &
    &  status='old') ! degC
    open (11,file='isi/tasmax_ipsl-cm5a-lr_2006-2010.txt',  &
    &  status='old') ! degC
    open (12,file='isi/pr_ipsl-cm5a-lr_2006-2010.txt'    ,  &
    &  status='old') ! mm/d
    open (13,file='isi/rhs_ipsl-cm5a-lr_2006-2010.txt'  ,  &
    &  status='old') ! %age
    open (14,file='isi/rsds_ipsl-cm5a-lr_2006-2010.txt'  ,  &
    &  status='old') ! W/m^2
    open (16,file='make_weather/windsp.txt',status='old') ! m/s
    open (17,file='isi/ps_ipsl-cm5a-lr_2006-2010.txt'   ,  &
    &  status='old') ! kPa
    open (18,file='isi/rlds_ipsl-cm5a-lr_2006-2010.txt'  ,  &
    &  status='old') ! W/m^2
  elseif (yeari == 2011) then
    open (10,file='isi/tasmin2011-2020.txt',status='old') ! degC
    open (11,file='isi/tasmax2011-2020.txt',status='old') ! degC
    open (12,file='isi/pr2011-2020.txt'    ,status='old') ! mm/d
    open (13,file='isi/huss2011-2020.txt'  ,status='old') ! g/kg
    open (14,file='isi/rsds2011-2020.txt'  ,status='old') ! W/m^2
    open (16,file='make_weather/windsp.txt',status='old') ! m/s
    open (17,file='isi/psl2011-2020.txt'   ,status='old') ! kPa
    open (18,file='isi/rlds2011-2020.txt'  ,status='old') ! W/m^2
  elseif (yeari == 2081) then
  open (10,file='isi/hadgem/tasmin_hadgem2-es_rcp8p5_2081-2099.txt', &
  &  status='old') ! degC
  open (11,file='isi/hadgem/tasmax_hadgem2-es_rcp8p5_2081-2099.txt',&
  &  status='old') ! degC
  open (12,file='isi/hadgem/pr_hadgem2-es_rcp8p5_2081-2099.txt'    ,  &
  &  status='old') ! mm/d
  open (13,file='isi/hadgem/rhs_hadgem2-es_rcp8p5_2081-2099.txt'  ,  &
  &  status='old') ! %age
  open (14,file='isi/hadgem/rsds_hadgem2-es_rcp8p5_2081-2099.txt'  ,  &
  &  status='old') ! W/m^2
  open (16,file='make_weather/windsp.txt',status='old') ! m/s
  open (17,file='isi/hadgem/ps_hadgem2-es_rcp8p5_2081-2099.txt'   ,  &
  &  status='old') ! kPa
  open (18,file='isi/hadgem/rlds_hadgem2-es_rcp8p5_2081-2099.txt'  ,  &
  &  status='old') ! W/m^2
  elseif (yeari == 2091) then
    open (10,file='isi/tasmin_ipsl-cm5a-lr_2091-2099.txt',  &
    &  status='old') ! degC
    open (11,file='isi/tasmax_ipsl-cm5a-lr_2091-2099.txt',  &
    &  status='old') ! degC
    open (12,file='isi/pr_ipsl-cm5a-lr_2091-2099.txt'    ,  &
    &  status='old') ! mm/d
    open (13,file='isi/rhs_ipsl-cm5a-lr_2091-2099.txt'  ,  &
    &  status='old') ! %age
    open (14,file='isi/rsds_ipsl-cm5a-lr_2091-2099.txt'  ,  &
    &  status='old') ! W/m^2
    open (16,file='make_weather/windsp.txt',status='old') ! m/s
    open (17,file='isi/ps_ipsl-cm5a-lr_2091-2099.txt'   ,  &
    &  status='old') ! kPa
    open (18,file='isi/rlds_ipsl-cm5a-lr_2091-2099.txt'  ,  &
    &  status='old') ! W/m^2
  else
    write (*,*) 'No isi climate for yeari ',yeari
    stop
  endif
  if (rcp == '8p5') then
    open (19,file='isi/RCP85_MIDYEAR_CONCENTRATIONS.DAT',status='old')
    do i = 1, 39
      read (19,*)
    enddo
    do kyr = 1765, yeari-1
      read (19,*) yr_rcp
    enddo
  endif
endif
!---------------------------------------------------------------------!
if (sclim) then
  open (50,file='OUT/'//sub_ouf//'/nee.txt',status='replace') !Laszlo
  open (51,file='OUT/'//sub_ouf//'/evapo.txt',status='replace')
  if (rcp == '8p5') then
    open (19,file='isi/RCP85_MIDYEAR_CONCENTRATIONS.DAT',status='old')
    do i = 1, 39
      read (19,*)
    enddo
    do kyr = 1765, yeari-1
      read (19,*) yr_rcp
    enddo
  endif
endif
!---------------------------------------------------------------------!
write (90,*) ite1
dtpoint = 1 ! Timepoint within day (n).
read (19,*) yr_rcp,co2eq_ppm,kyoto_co2eq_ppm,co2_ppm
!---------------------------------------------------------------------!
! Co2 values for palaeo-runs are in co2.txt.
!---------------------------------------------------------------------!
if (ooaf) then
  open (20,file='co2.txt',status='old')
  read (20,*)
  read (20,*)
  do i = 1, 62
    read (20,*) yr_ooaf,co2_ooaf
    if (yr_ooaf == -yeari/1000) co2_ppm = co2_ooaf
  enddo
  write (*,*) 'co2_ppm ',co2_ppm
  close (20)
endif
!---------------------------------------------------------------------!
! Jump to restart point in climate and CO2 files.
!---------------------------------------------------------------------!
if (restart) then
  do jday = 1, ((datei-1) + jdendofm (monthi-1))
    do j0 = 1, jm
      read (10,*) (tmin_ij    (i0,j0),i0=1,im) ! degC
      read (11,*) (tmax_ij    (i0,j0),i0=1,im) ! degC
      read (12,*) (raini_ij   (i0,j0),i0=1,im) ! mm/d
      read (13,*) (vap_ij     (i0,j0),i0=1,im) ! kPa or %age (isi)
      ! MJ/m^2/d or W/m^2 (isi)
      read (14,*) (radi_ij    (i0,j0),i0=1,im)
      read (15,*) (rc_ij      (i0,j0),i0=1,im) ! MJ/m^2/d
      if (jday == 365) rewind (15)
      read (16,*) (windsp_ij  (i0,j0),i0=1,im) ! m/s
      read (17,*) (pressure_ij(i0,j0),i0=1,im) ! MPa or kPa (isi)
      read (18,*) (ldi_ij     (i0,j0),i0=1,im) ! W/m^2
    enddo
  enddo
endif
!---------------------------------------------------------------------!
raint                = 0.0d0
rain_ij      (:,:,:) = 0.0d0
aij_wtop     (:,:)   = 0.0d0
aij_wtot     (:,:)   = 0.0d0
aij_betad    (:,:)   = 0.0d0
aij_gpp      (:,:)   = 0.0d0
aij_npp      (:,:)   = 0.0d0
aij_nee      (:,:)   = 0.0d0
aij_cm       (:,:,:) = 0.0d0
aij_lai      (:,:,:) = 0.0d0
aij_fpar     (:,:,:) = 0.0d0
aij_esoil    (:,:)   = 0.0d0
aij_intercep (:,:)   = 0.0d0
aij_transp   (:,:)   = 0.0d0
aij_trunoff  (:,:)   = 0.0d0
aij_srunoff  (:,:)   = 0.0d0
jyear                = 0
!---------------------------------------------------------------------!
do while (itime < itimee)
  !--------------------------------------------------------------------!
  ! If this is a new year, read atmospheric CO2 (ppm).
  !--------------------------------------------------------------------!
  if ((yeari + itime / (nday * jdpery)) > jyear) then
    if (ooaf) then
      write (99,*) 'yeari, co2_ppm ',yeari,co2_ppm
    else
      read (19,*) yr_rcp,co2eq_ppm,kyoto_co2eq_ppm,co2_ppm
    endif
    raint = 0.0d0
  endif
  !--------------------------------------------------------------------!
  ! Current year (yr).
  !--------------------------------------------------------------------!
  jyear = yeari + itime / (nday * jdpery)
  !--------------------------------------------------------------------!
  ! Julian day (day).
  !--------------------------------------------------------------------!
  jday = 1 + mod (itime / nday, 365)
  !--------------------------------------------------------------------!
  ! Current month (mo).
  !--------------------------------------------------------------------!
  jmon = 1
  do while (jday > jdendofm (jmon))
    jmon = jmon + 1
  enddo
  !--------------------------------------------------------------------!
  ! Keep note of run.
  !--------------------------------------------------------------------!
  !write (*,*) itime,jyear,jday,jmon
  !--------------------------------------------------------------------!
  cosday = cos (twopi / edpery * float (jday))
  sinday = sin (twopi / edpery * float (jday))
  !--------------------------------------------------------------------!
  ! If this is a new day, read in daily weather and generate
  ! sub-daily weather.
  !--------------------------------------------------------------------!
  if (mod (itime,nday) == 0) then
    if (.not. (sclim)) then
      do j0 = 1, jm
        read (10,*) (tmin_ij    (i0,j0),i0=1,im) ! degC
        read (11,*) (tmax_ij    (i0,j0),i0=1,im) ! degC
        read (12,*) (raini_ij   (i0,j0),i0=1,im) ! mm/d
        read (13,*) (vap_ij     (i0,j0),i0=1,im) ! kPa or %age (isi)
        ! MJ/m^2/d or W/m^2 (isi)
        read (14,*) (radi_ij    (i0,j0),i0=1,im)
        read (15,*) (rc_ij      (i0,j0),i0=1,im) ! MJ/m^2/d
        if (jday == 365) rewind (15)
        read (16,*) (windsp_ij  (i0,j0),i0=1,im) ! m/s
        read (17,*) (pressure_ij(i0,j0),i0=1,im) ! MPa or kPa (isi)
        read (18,*) (ldi_ij     (i0,j0),i0=1,im) ! W/m^2
      enddo
      if (isi) then
        do j0 = 1, jm
          do i0 = 1, im
            tmean = (tmin_ij(i0,j0)+tmax_ij(i0,j0))/2.0 + tf
            qsat = qsat_modele (tmean,lhe,pressure_ij(i0,j0)*10.0)
            vap_ij (i0,j0) = vap_ij (i0,j0) * qsat / 100.0 ! kg/kg
          enddo
        enddo
        ! Convert radiation from W/m^2 to MJ/m2/d.
        radi_ij     (:,:) = radi_ij     (:,:) * 86400.0d0 / 1.0d6
        pressure_ij (:,:) = pressure_ij (:,:) / 1.0d3 ! kPa > MPa
      endif
      do j0 = 1, jm
        do i0 = 1, im
          !------------------------------------------------------------!
          if (fearth (i0,j0) > 0.0d0) then
            !----------------------------------------------------------!
            ! Generate sub_daily weather.
            !----------------------------------------------------------!
            dtmin (jday) = tmin_ij (i0,j0) * 9.0d0 / 5.0d0 + 32.0d0
            dtmax (jday) = tmax_ij (i0,j0) * 9.0d0 / 5.0d0 + 32.0d0
            rc  (jday) = rc_ij   (i0,j0) / 0.04189d0
            rad (jday) = radi_ij (i0,j0) / 0.04189d0
            tsteph = 24.0d0 * dtsrc / sday
            lat (i0,j0) = (float (j0) - a_lat) / b_lat
            latr = rfac * lat (i0,j0)
            pressure (jday) = pressure_ij (i0,j0) * 1.0d6 ! MPa>Pa
            windsp (jday) = windsp_ij (i0,j0) ! m/s
            rain   (jday) = raini_ij  (i0,j0) ! mm/d
            vap    (jday) = vap_ij    (i0,j0) ! kPa or kg/kg (isi)
            !----------------------------------------------------------!
            ! Algorithm to check is OK as used in make_weather.
            !----------------------------------------------------------!
            !if (rc (jday) > 0.0d0) then
            !  cloud = 1.0d0 - rad (jday) / rc (jday)
            !else
            !  cloud = 0.0d0
            !endif
            !cloud = min (1.0d0,cloud)
            !cloud = max (0.0d0,cloud)
            !emm = 0.787812 * (1.0d0 - (cloud ** 6)) + &
            !&     0.971812 * (cloud ** 4)
            !tak_ld = (tmin_ij (i0,j0) + tmax_ij (i0,j0)) / 2.0d0 + tf
            !ldi_ij (i0,j0) = emm * stblz * (tak_ld ** 4)
            !----------------------------------------------------------!
            call daily_climate_forcing_wg
            !----------------------------------------------------------!
            rain_ij    (i0,j0,:) = ppt     (:) ! m/s
            pedn_ij    (i0,j0,:) = p       (:) ! Pa
            tak_ij     (i0,j0,:) = tak     (:) ! K
            qsrf_ij    (i0,j0,:) = qa      (:) ! kg/kg
            rad_ij     (i0,j0,:) = rg      (:) ! W/m^2
	    ld_ij      (i0,j0,:) = ld      (:) ! W/m^2
            sbeta_ij   (i0,j0,:) = sbeta   (:) ! (ratio)
            ipar_df_ij (i0,j0,:) = ipar_df (:) ! umol(PAR)/m^2/s
            ipar_dr_ij (i0,j0,:) = ipar_dr (:) ! umol(PAR)/m^2/s
           !-----------------------------------------------------------!
    if ((i0 == olon) .and. (j0 == olat)) then
      raint=raint+ppt(12)*86400.0d3
      write (*,*) ldi_ij(i0,j0)
      write (*,*) raint
      write (*,*) ppt(12)*86400.0d3,p(12),tak(12)-tf,qa(12),rg(12),&
      &           sum(wvege(:,i0,j0))/3.5d0,co2_ppm,&
      &           tearth(i0,j0)
      write (*,*) 1.0e3*aij_gpp(i0,j0),1.0e3*aij_npp(i0,j0)
    endif
          endif
        enddo
      enddo
    else
      !----------------------------------------------------------------!
      i0 = olon
      j0 = olat
      !----------------------------------------------------------------!
      read (21,*) kyr, kday, tmin_ij(i0,j0), tmax_ij(i0,j0), &
      & raini_ij(i0,j0),&
      & vap_ij(i0,j0),radi_ij(i0,j0),rc_ij(i0,j0),windsp_ij(i0,j0),&
      & pressure_ij(i0,j0),ldi_ij(i0,j0)
      !write (*,*) kyr, kday, tmin_ij(i0,j0), tmax_ij(i0,j0), &
      !& raini_ij(i0,j0),&
      !& vap_ij(i0,j0),radi_ij(i0,j0),rc_ij(i0,j0),windsp_ij(i0,j0),&
      !& pressure_ij(i0,j0),ldi_ij(i0,j0)
      !stop
      !----------------------------------------------------------------!
      dtmin (jday) = tmin_ij (i0,j0) * 9.0d0 / 5.0d0 + 32.0d0
      dtmax (jday) = tmax_ij (i0,j0) * 9.0d0 / 5.0d0 + 32.0d0
      rc  (jday) = rc_ij   (i0,j0) / 0.04189d0
      rad (jday) = radi_ij (i0,j0) / 0.04189d0
      windsp (jday) = windsp_ij (i0,j0) ! m/s
      rain   (jday) = raini_ij  (i0,j0) ! mm/d
      vap    (jday) = vap_ij    (i0,j0) ! kPa
      !----------------------------------------------------------------!
      pressure(jday) = pressure_ij(i0,j0) * 1.0d6 ! MPa>Pa (isi)
      tsteph = 24.0d0 * dtsrc / sday
      lat (i0,j0) = latw
      latr = rfac * lat (i0,j0)
      !----------------------------------------------------------------!
      call daily_climate_forcing_wg
      !----------------------------------------------------------------!
      rain_ij (i0,j0,:) = ppt (:) ! m/s
      pedn_ij (i0,j0,:) = p   (:) ! Pa
      tak_ij  (i0,j0,:) = tak (:) ! K
      qsrf_ij (i0,j0,:) = qa  (:) ! kg/kg
      rad_ij  (i0,j0,:) = rg  (:) ! W/m^2
!      write(*,*) rad_ij (i0,j0,:)	   
!      stop
      ld_ij   (i0,j0,:) = ld  (:) ! W/m^2
      sbeta_ij   (i0,j0,:) = sbeta   (:) ! (ratio)
      ipar_df_ij (i0,j0,:) = ipar_df (:) ! umol(PAR)/m^2/s
      ipar_dr_ij (i0,j0,:) = ipar_dr (:) ! umol(PAR)/m^2/s
      !----------------------------------------------------------------!
    endif
  !  write (*,*) 'Running yr=',jyear,'day=',jday,&
  !  &     100.0*float(itime-ihri+1)/float(itimee+1),'% done'
    if ((i0 == olon) .and. (j0 == olat)) then
      raint=raint+ppt(12)*86400.0d3
  !    write (*,*) 'ldi=', ldi_ij(i0,j0)
  !    write (*,*) 'raint=', raint
  !    write (*,*) 'ppt=', ppt(12)*86400.0d3
  !    write (*,*) 'p=', p(12)
  !    write (*,*) 'tak=', tak(12)-tf
  !    write (*,*) 'qa=', qa(12)
  !    write (*,*) 'rg=', rg(12)
  !    write (*,*) 'wvege=', sum(wvege(:,i0,j0))/3.5d0
  !    write (*,*) 'co2ppm=', co2_ppm
  !    write (*,*) 'tearth=', tearth(i0,j0)
  !    write (*,*) 'gpp=', 1.0e3*aij_gpp(i0,j0)
  !    write (*,*) 'npp=', 1.0e3*aij_npp(i0,j0)      
  !    write (*,*) '******************************************************************'
  !    write (*,*) itime
    endif
  endif ! Beginning of day
  !--------------------------------------------------------------------!
  ! Get precipitation amount in timestep (from kg m/s) (kg/m^2).
  !--------------------------------------------------------------------!
  prec (:,:) = 1.0d3 * rain_ij (:,:,dtpoint) * dtsrc
  !--------------------------------------------------------------------!
  ! Get surface layer atmospheric pressure (hPa from Pa).
  !--------------------------------------------------------------------!
  pedn (:,:) = pedn_ij (:,:,dtpoint) / 1.0d2
  !--------------------------------------------------------------------!
  ! Get surface layer temperature (K).
  !--------------------------------------------------------------------!
  t (:,:) = tak_ij (:,:,dtpoint)
  !--------------------------------------------------------------------!
  ! Get surface SW radiation (W/m^2).
  !--------------------------------------------------------------------!
  fsf (:,:) = rad_ij (:,:,dtpoint) * (1.0d0 - 0.061d0)
  !--------------------------------------------------------------------!
  ! Calculate mass of first atmosphere layer (kg/m^2).
  !--------------------------------------------------------------------!
  am (:,:) = pedn (:,:) * dsig * 1.0d2 * bygrav
  !--------------------------------------------------------------------!
  ! Get surface layer specific humidity (kg/kg).
  !--------------------------------------------------------------------!
  q_surf (:,:) = qsrf_ij (:,:,dtpoint)
  !--------------------------------------------------------------------!
  ! Get surface downward longwave radiation (W/m^2).
  !--------------------------------------------------------------------!
  trhr (:,:) = ld_ij (:,:,dtpoint)
  !--------------------------------------------------------------------!
  ! Sine of solar elevation (ratio).
  !--------------------------------------------------------------------!
  sbeta_can (:,:) = sbeta_ij (:,:,dtpoint)
  !--------------------------------------------------------------------!
  ! Downwelling diffuse PAR (umol/m^2/s).
  !--------------------------------------------------------------------!
  ipar_df_can (:,:) = ipar_df_ij (:,:,dtpoint)
  !--------------------------------------------------------------------!
  ! Downwelling direct PAR (umol/m^2/s).
  !--------------------------------------------------------------------!
  ipar_dr_can (:,:) = ipar_dr_ij (:,:,dtpoint)
  !--------------------------------------------------------------------!
  ! Loop over timesteps, executed nisurf times every hour.
  !--------------------------------------------------------------------!
  do ns = 1, nisurf
    call earth
  enddo
  !--------------------------------------------------------------------!
  itime = itime + 1 ! dtsrc-steps since 1/1/iyear1
  !--------------------------------------------------------------------!
  if (mod (itime, nday) == 0) then
    dtpoint = 1
  else
    dtpoint = dtpoint + 1 ! Timepoint within day (n).
  endif
  !--------------------------------------------------------------------!
  ! Keep note of run at noon each day.
  !--------------------------------------------------------------------!
  write (*,*) 'jyear jday itime ',jyear,jday,itime !Laszlo
  write(*,*) npp !Laszlo
  write(*,*) fb*evap_tot(1)+fv*evap_tot(2) !Laszlo
  write(50,*) npp !Laszlo
  write(51,*) fb*evap_tot(1)+fv*evap_tot(2) !Laszlo

  if (mod (itime+nday/2, nday) == 0) then
  write (90,'(3i8,25f10.4)') jyear,                   & !  1
  & jday,                                             & !  2
  & itime,                                            & !  3
  & sum(1000.0*wbare(1:nsoil_layers_max,olon,olat)),  & !  4
  & sum(1000.0*wvege(1:nsoil_layers_max,olon,olat)),  & !  5
  & prec(olon,olat),                                  & !  6
  & betad_ij(olon,olat),                              & !  7
  & aij_gpp(olon,olat),                               & !  8
  & 1000.0*wvege(1,olon,olat)/70.0d0,                 & !  9
  & trhr(olon,olat),                                  & ! 10
  & q_surf(olon,olat),                                & ! 11
  & am(olon,olat),                                    & ! 12
  & fsf(olon,olat),                                   & ! 13
  & t(olon,olat)-tfrz                                   ! 14
  !write (*,*) 'jyear jday itime ',jyear,jday,itime
  endif
  !--------------------------------------------------------------------!
  !write (88,*) cm (olon,olat,:)
  !--------------------------------------------------------------------!
  ! Write monthly-mean diagnostic and current state for restart.
  !--------------------------------------------------------------------!
  if ((mod (itime, nday) == 0) .and. (jday == jdendofm (jmon))) then
    !------------------------------------------------------------------!
    if (jyear > 0) then
      write (adate(1:7),'(a3,i4.4)') amon(jmon),jyear
    else
      write (adate(1:7),'(a3,i4.4)') amon(jmon),jyear+66000
    endif
    !------------------------------------------------------------------!
    open (20,file='OUT/'//sub_ouf//'/'//adate(1:7)//'.rsf',status='unknown')
    !------------------------------------------------------------------!
    do j0 = 1, jm
      !----------------------------------------------------------------!
      ! Water on canopy (m).
      !----------------------------------------------------------------!
      write (20,'(720f8.4)') (wvege (0,i0,j0),i0=1,im)
      !----------------------------------------------------------------!
      ! Heat in canopy (J/m^2).
      !----------------------------------------------------------------!
      write (20,'(720f10.4)') (htvege (0,i0,j0)/1.0d6,i0=1,im)
      !----------------------------------------------------------------!
      ! Carbon in plant pools (kgC/m^2).
      !----------------------------------------------------------------!
      write (20,'(720f8.4)') (cm (i0,j0,:),i0=1,im)
      !----------------------------------------------------------------!
      ! Internal leaf CO2 (mol/m^3).
      !----------------------------------------------------------------!
      write (20,'(720f8.4)') (veg_ci_ij (i0,j0),i0=1,im)
      !----------------------------------------------------------------!
      ! Canopy surface mixing ratio (kg/kg).
      !----------------------------------------------------------------!
      write (20,'(720f8.4)') (qfol_ij         (i0,j0)    ,i0=1,im)
      !----------------------------------------------------------------!
      write (20,'(720f8.4)') (snowbv          (1,i0,j0)  ,i0=1,im)
      write (20,'(720f8.4)') (snowbv          (2,i0,j0)  ,i0=1,im)
      !----------------------------------------------------------------!
      write (20,'(720f8.4)') (nsnow_layers_ij (1,i0,j0)  ,i0=1,im)
      write (20,'(720f8.4)') (nsnow_layers_ij (2,i0,j0)  ,i0=1,im)
      !----------------------------------------------------------------!
      write (20,'(720f8.4)') (wsn_ij          (:,1,i0,j0),i0=1,im)
      write (20,'(720f8.4)') (wsn_ij          (:,2,i0,j0),i0=1,im)
      !----------------------------------------------------------------!
      write (20,'(720f8.4)') (fr_snow_ij      (1,i0,j0)  ,i0=1,im)
      write (20,'(720f8.4)') (fr_snow_ij      (2,i0,j0)  ,i0=1,im)
      !----------------------------------------------------------------!
    enddo
    do nsoil_layer = 1, nsoil_layers_max
      do j0 = 1, jm
        !--------------------------------------------------------------!
        ! Water in bare soil layers (m).
        !--------------------------------------------------------------!
        write (20,'(720f8.4)') (wbare  (nsoil_layer,i0,j0),i0=1,im)
        !--------------------------------------------------------------!
        ! Water in vegetated soil layers (m).
        !--------------------------------------------------------------!
        write (20,'(720f8.4)') (wvege  (nsoil_layer,i0,j0),i0=1,im)
        !--------------------------------------------------------------!
        ! Heat in bare soil layers (m).
        !--------------------------------------------------------------!
        write (20,'(720f10.4)') (htbare (nsoil_layer,i0,j0)/1.0d6, &
        &                      i0=1,im)
        !--------------------------------------------------------------!
        ! Heat in vegetated soil layers (m).
        !--------------------------------------------------------------!
        write (20,'(720f10.4)') (htvege (nsoil_layer,i0,j0)/1.0d6, &
        &                      i0=1,im)
      enddo
    enddo
    !------------------------------------------------------------------!
    close (20)
    !------------------------------------------------------------------!
    open (20,file='OUT/'//sub_ouf//'/'//adate(1:7)//'.acc',status='unknown')
    !------------------------------------------------------------------!
    do j0 = 1, jm
      do i0 = 1, im
        !--------------------------------------------------------------!
        ! Diag. of global field of monthly-mean soil water (m^3/m^3).
        !--------------------------------------------------------------!
        aij_wtop (i0,j0) = aij_wtop (i0,j0) / &
        &                  float (nd (jmon) * 24 * nisurf)
        !--------------------------------------------------------------!
        ! Diag. of global field of monthly-mean soil water (fraction).
        !--------------------------------------------------------------!
        aij_wtot (i0,j0) = aij_wtot (i0,j0) / &
        &                  float (nd (jmon) * 24 * nisurf)
        !--------------------------------------------------------------!
        ! Diag. of global field of monthly-mean betad (fraction).
        !--------------------------------------------------------------!
        aij_betad (i0,j0) = aij_betad (i0,j0) / &
        &                   float (nd (jmon) * 24 * nisurf)
        !--------------------------------------------------------------!
        ! Diag. of global fields of monthly-mean carbon (kgC/m^2).
        !--------------------------------------------------------------!
        aij_cm (i0,j0,:) = aij_cm (i0,j0,:) / &
        &                  float (nd (jmon) * 24 * nisurf)
        !--------------------------------------------------------------!
        ! Diag. of global fields of plant LAI (m^2/m^2).
        !--------------------------------------------------------------!
        aij_lai (i0,j0,:) = aij_lai (i0,j0,:) / &
        &                   float (nd (jmon) * 24 * nisurf)
        !--------------------------------------------------------------!
        ! Diag. of global fields of plant fPAR (fraction).
        !--------------------------------------------------------------!
        aij_fpar (i0,j0,:) = aij_fpar (i0,j0,:) / &
        &                    float (nd (jmon) * 24 * nisurf)
        !--------------------------------------------------------------!
        ! If find have to set limits, do so here.
        ! f7.3: -9.999 to 99.999
        ! f7.1: -999.9 to 9999.9
        !--------------------------------------------------------------!
        !--------------------------------------------------------------!
      enddo ! Volumetric water content in top layer (m^3/m^3).
      write (20,'(720f7.3)') (aij_wtop  (i0,j0), i0 = 1, im) ! 1
    !------------------------------------------------------------------!
    enddo
    do j0 = 1, jm ! fraction
      write (20,'(720f7.3)') (aij_wtot  (i0,j0), i0 = 1, im) ! 2
    enddo
    do j0 = 1, jm ! fraction
      write (20,'(720f7.3)') (aij_betad (i0,j0), i0 = 1, im) ! 3
    enddo
    do j0 = 1, jm ! kgC/m^2/mo.
      write (20,'(720f7.3)') (aij_gpp   (i0,j0), i0 = 1, im) ! 4
    enddo
    do j0 = 1, jm ! kgC/m^2/mo.
      write (20,'(720f7.3)') (aij_npp   (i0,j0), i0 = 1, im) ! 5
    enddo
    do j0 = 1, jm ! kgC/m^2/mo.
      write (20,'(720f7.3)') (aij_nee   (i0,j0), i0 = 1, im) ! 6
    enddo
    do icpool = 1, ncpools ! 7-10 kgC/m^2
      do j0 = 1, jm
        write (20,'(720f7.3)') (aij_cm (i0,j0,icpool), i0 = 1, im)
      enddo
    enddo
    do kplant = 1, nplant ! 11-18 m^2/m^2
      do j0 = 1, jm
        write (20,'(720f7.3)') (aij_lai (i0,j0,kplant), i0 = 1, im)
      enddo
    enddo
    do kplant = 1, nplant ! 19-26 fraction
      do j0 = 1, jm
        write (20,'(720f7.3)') (aij_fpar (i0,j0,kplant), i0 = 1, im)
      enddo
    enddo
    do j0 = 1, jm ! kg/m^2/mo.
      write (20,'(720f7.1)') (aij_esoil (i0,j0), i0 = 1, im) ! 27
    enddo ! kg/m^2/mo.
    do j0 = 1, jm !kg/m^2/mo.
      write (20,'(720f7.1)') (aij_intercep (i0,j0), i0 = 1, im) ! 28
    enddo
    do j0 = 1, jm !kg/m^2/mo.
      write (20,'(720f7.1)') (aij_transp (i0,j0), i0 = 1, im) ! 29
    enddo
    do j0 = 1, jm !kg/m^2/mo.
      write (20,'(720f7.1)') (aij_trunoff (i0,j0), i0 = 1, im) ! 30
    enddo
    do j0 = 1, jm !kg/m^2/mo.
      write (20,'(720f7.1)') (aij_srunoff (i0,j0), i0 = 1, im) ! 31
    enddo
    !------------------------------------------------------------------!
    ! Fractional coverage of cover types.
    !------------------------------------------------------------------!
    do kveg = 1, nveg ! 32-41
      do j0 = 1, jm
        write (20,'(720f7.3)') (vdata (i0,j0,kveg), i0 = 1, im)
      enddo
    enddo
    !------------------------------------------------------------------!
    close (20)
    !------------------------------------------------------------------!
    aij_wtop     (:,:)   = 0.0d0
    aij_wtot     (:,:)   = 0.0d0
    aij_betad    (:,:)   = 0.0d0
    aij_gpp      (:,:)   = 0.0d0
    aij_npp      (:,:)   = 0.0d0
    aij_nee      (:,:)   = 0.0d0
    aij_cm       (:,:,:) = 0.0d0
    aij_lai      (:,:,:) = 0.0d0
    aij_esoil    (:,:)   = 0.0d0
    aij_intercep (:,:)   = 0.0d0
    aij_transp   (:,:)   = 0.0d0
    aij_trunoff  (:,:)   = 0.0d0
    aij_srunoff  (:,:)   = 0.0d0
    !------------------------------------------------------------------!
  endif
  !-------------------------------------------------------------------!
enddo ! End of main loop (while (itime < itimee))
!---------------------------------------------------------------------!
! Close daily weather files.
!---------------------------------------------------------------------!
do i = 1, nfiles_in
  close (i+9)
enddo
!----------------------------------------------------------------------!
! Close daily climate data file if used.
!----------------------------------------------------------------------!
if (sclim) then
  close (21)
endif
close (50) !Laszlo
close (51) !Laszlo
!----------------------------------------------------------------------!
! Diagnostic of global field of mean soil water (m).
!----------------------------------------------------------------------!
open (20, file='OUT/'//sub_ouf//'/w_field.txt',status='unknown')
do j0 = 1, jm
  do i0 = 1, im
    fb = afb (i0,j0)
    fv = 1.0d0 - fb
    wtot (i0,j0) = fb * sum (wbare (1:nsoil_layers_max,i0,j0)) + &
    &              fv * sum (wvege (1:nsoil_layers_max,i0,j0))
  enddo
  write (20,'(720f7.3)') (wtot (i0,j0), i0 = 1, im)
enddo
close (20)
!----------------------------------------------------------------------!
if (itime >= itimee) then
  write (90,*) 'Terminated normally (reached maximum time)'
  write ( *,*) 'Terminated normally (reached maximum time)'
endif
write (90,*) 'Finished'
write ( *,*) 'Finished'
close (88)
close (89)
close (90)
close (98)
close (99)
stop
!----------------------------------------------------------------------!
100 write (90,*) 'Problem opening global fields file. iostat =',iostat
write (90,*) 'Simulation terminated'
close (88)
close (89)
close (90)
close (98)
close (99)
stop
!----------------------------------------------------------------------!
end program hybrid8
!======================================================================!
