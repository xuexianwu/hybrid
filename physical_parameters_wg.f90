!======================================================================!
MODULE PHYSICAL_PARAMETERS_WG
!----------------------------------------------------------------------!
USE CONTROL_DEFINITIONS_WG
!----------------------------------------------------------------------!
IMPLICIT NONE
save
!----------------------------------------------------------------------!
!  Physical constants.
!----------------------------------------------------------------------!
INTEGER, PARAMETER :: dyr   = 365     ! No. days in year             (n)
INTEGER, PARAMETER :: nmo   = 12      ! No. months in year           (n)
REAL, PARAMETER :: scs  = 1367.0    ! Solar constant (W/m2)
REAL, PARAMETER :: stblz = 5.67e-8  ! Stefan-Boltzmann constant 
                                    ! (W/m2/K4)
!----------------------------------------------------------------------!
!  Internal physical parameters.
!----------------------------------------------------------------------!
REAL :: latr            ! Latitude of simulated site            (radian)
REAL :: rfrac_mm  (nmo) ! Wet days in month                   (fraction)
REAL :: pwet_mm   (nmo) ! Rain per wet day                        (mm/d)
REAL :: tmax_mm   (nmo) ! Mean daily maximum temperature            (oC)
REAL :: tmin_mm   (nmo) ! Mean daily minimum temperature            (oC)
REAL :: solrad_mm (nmo) ! Mean daily solar radiation           (MJ/m2/d)
REAL :: dq_mm     (nmo) ! Atmospheric moisture at 09h00            (MPa)
REAL :: wind_mm   (nmo) ! Mean wind speed                          (m/s)
REAL :: p_mm      (nmo) ! Mean atmospheric pressure                 (Pa)
!----------------------------------------------------------------------!
REAL   , PARAMETER :: R     = 8.3144  ! Gas constant           (J/mol/K)
REAL   , PARAMETER :: Tfrz  = 273.16  ! Freezing temp of water       (K)
INTEGER, DIMENSION(nmo), PARAMETER :: ni = (/  31, 59, 90,120,151,181, &
&                                             212,243,273,304,334,365 /)
! No, days in each month                                         (d/mo)
integer, parameter :: nd (nmo) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
INTEGER, DIMENSION(nmo), PARAMETER :: nii = (/  31, 60, 91,121,152,182,&
&                                              213,244,274,305,335,366/)
!----------------------------------------------------------------------!
END MODULE PHYSICAL_PARAMETERS_WG
!======================================================================!
