!======================================================================!
MODULE CONTROL_DEFINITIONS_WG
!----------------------------------------------------------------------!
IMPLICIT NONE
save
!----------------------------------------------------------------------!
!  Parameters set in 'driver_wg.txt'.
!----------------------------------------------------------------------!
REAL    :: tsteph  ! Timestep for simulation                        (hr)
LOGICAL :: subd    ! True if want sub-daily output                 (log)
!----------------------------------------------------------------------!
!  Derived control parameters.
!----------------------------------------------------------------------!
INTEGER :: ntpoints ! No. time points in day                         (n)
!----------------------------------------------------------------------!
!  Control variables.
!----------------------------------------------------------------------!
INTEGER :: kyr                           ! Year count               (yr)
INTEGER :: kmo                           ! Month count              (mo)
INTEGER :: kday                          ! Day count                 (d)
INTEGER :: ktpoint                       ! Time point count          (n)
INTEGER :: K (4)                         ! Weather generator seed    (x)
INTEGER :: IP                            ! Weather gen variable      (x)
REAL    :: XIM1 (4)                      ! Weather gen variable      (x)
!----------------------------------------------------------------------!
END MODULE CONTROL_DEFINITIONS_WG
!======================================================================!
