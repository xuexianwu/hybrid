!======================================================================!
MODULE WGEN_PARAMS_WG
!----------------------------------------------------------------------!
USE CONTROL_DEFINITIONS_WG
!----------------------------------------------------------------------!
IMPLICIT NONE
save
!----------------------------------------------------------------------!
!  Parameters for daily weather generator, WGEN.
!----------------------------------------------------------------------!
REAL   , DIMENSION(12)  :: pww    ! Probability of wet day
                                  ! following wet day (fraction)
REAL   , DIMENSION(12)  :: pwd    ! Probability of wet day
                                  ! following dry day (fraction)
REAL   , DIMENSION(12)  :: alphap ! Rainfall Gamma
                                  ! distribution shape parameter
REAL   , DIMENSION(12)  :: beta   ! Rainfall Gamma
                                  ! distribution scale parameter
REAL   , DIMENSION(366) :: txm    ! Mean of tmax? (?)
REAL   , DIMENSION(366) :: txs    ! ?
REAL   , DIMENSION(366) :: txm1   ! ?
REAL   , DIMENSION(366) :: txs1   ! ?
REAL   , DIMENSION(366) :: tnm    ! ?
REAL   , DIMENSION(366) :: tns    ! ?
REAL   , DIMENSION(366) :: rmo    ! ?
REAL   , DIMENSION(366) :: rso    ! ?
REAL   , DIMENSION(366) :: rm1    ! ?
REAL   , DIMENSION(366) :: rs1    ! ?
INTEGER                 :: kgen = 1 ! WGEN option
REAL   , DIMENSION(366) :: rc     ! Potential radiation (langleys/day)
INTEGER                 :: idays  ! Generation option
REAL   , DIMENSION(366) :: vps    ! ?
REAL   , DIMENSION(366) :: vpm    ! ?
REAL   , DIMENSION( 12) :: tcfmax ! ?
REAL   , DIMENSION( 12) :: tcfmin ! ?
REAL   , DIMENSION( 12) :: rcf    ! ?
REAL   , DIMENSION( 12) :: cc     ! ?
REAL   , DIMENSION( 12) :: kay    ! ?
!----------------------------------------------------------------------!
!  To replace penn.param.
CHARACTER*69 :: pHEADER
INTEGER :: pNYR,pKGEN,pKTCF,pKRCF
REAL :: pALAT
REAL :: pCC(12),pKAY(12),pPWW(12),pPWD(12),pALPHAP(12),pBETA(12)
REAL :: pTXMD,pATX,pCVTX,pACVTX,pPEAKTX,pTXMW,pTN,pATN
REAL :: pCVTN,pACVTN,pPEAKTN
REAL :: pRMD,pAR,pPEAKSR,pCVRD,pACVRD,pRMW,pCVRW,pACVRW
REAL :: ppVP,pAVP,pPEAKVP,pCVVP,pACVVP
!----------------------------------------------------------------------!
END MODULE WGEN_PARAMS_WG
!======================================================================!
