!======================================================================!
	SUBROUTINE WEASHO_WG
!	SUBROUTINE WEASHO (  PWW,PWD,ALPHAP,BETA,TXM,TXS,TXM1,TXS1,TNM,&
!&                          TNS,RMO,RSO,RM1,RS1,RAIN, dTMAX, dTMIN,RAD,&
!&                          KGEN,RC,IDAYS,VPS,VPM,VAP,TCFMAX,    &
!&                          TCFMIN,RCF,CC,KAY,WINDSP)
!----------------------------------------------------------------------!
!	WEASHO.FOR = "SHORT METHOD" VERSION OF WGEN MODEL OF C.W.RICHARDSON
!				AND D.A. WRIGHT, USDA
!		     SAME AS WGEN, BUT HAS SOME PARAMETERS VARIABLE 
!		                ADDED OR CHANGED
!	USE ESTPAR.FOR (PART OF SIMMETEO.FOR) TO CREATE INPUT FILE FROM 
!                             MONTHLY SUMMARY DATA
!	DEVELOPED BY BAOQUAN LI (AUG. 1988)   
!	THE MINOR FORMAT CHANGES:
!		A) ADD WEATHER FACTORS HUMIDITY & WIND SPEED
!		B) INPUT PARAMETERS ARE F OR C FOR TEMPERATURES,
!		    MM OR IN FOR RAINFALL, 10**6 J/M**2/DAY OR LANGLEY/DAY
!                   FOR SOLRAD, MPAR (VAPOR PRESSURE) OR C/F (DEWPOINT)
!		    FOR HUMIDITY, METER/SEC OR MILE/HOUR FOR WIND SPEED 
!		C) OUTPUT STANDARD UNITS ARE C FOR TEMPS, MM FOR RAINFALL, 
!	    	   10**6 J/M**2/DAY FOR SOLRAD ,MPAR (VAPOR PRESSURE)
!		   FOR HUMIDITY, AND METER/SEC FOR WIND SPEED. OUTPUT
!                  UNITS CAN BE CHOSEN AS THE INPUT UNITS.  
!		D) OUTPUT CONSISTENT WITH SOME PROGRAMS
!		        PREVIOUS MODIFICATIONS BY 
!                       1) I. SUPIT, CABO, WAGENINGEN, THE NETHERLANDS
!       	   	2) JILL SHORE AUBURN,U.C.D.,CA.,U. S. A.
!----------------------------------------------------------------------!
USE WGEN_params_WG
USE CONTROL_DEFINITIONS_WG
USE GEN_ATMOS_WG
use physical_parameters_WG
!----------------------------------------------------------------------!
!INTEGER, PARAMETER :: ndays = 365 ! No. days in year
!INTEGER, PARAMETER :: nmo   = 12  ! No. months in year
!INTEGER :: ntpoints               ! No. time points in day
!INTEGER :: nsites                 ! No. sites in simulation
!INTEGER :: ksite                  ! Site count
!INTEGER :: kyr                    ! Year count
!INTEGER :: kmo                    ! Month count
!INTEGER :: kday                   ! Day count
!INTEGER :: ktpoint                ! Time point count
!INTEGER :: npo_month              ! Number of time points in month
!INTEGER :: l                      ! Canopy layer
!INTEGER :: nit
!INTEGER :: i
!INTEGER, DIMENSION(4) :: K        ! Weather generator seed
!INTEGER :: IP                     ! Weather generator parameter
!INTEGER :: ic   ! Index for Table of coefficients (a)
!INTEGER :: surf ! Surface type (1 = bare; 2 = vegetated)
!INTEGER,DIMENSION(nmo), PARAMETER :: nd = (/  31, 28, 31, 30, 31, 30, &
!&                                             31, 31, 30, 31, 30, 31 /)
!INTEGER,DIMENSION(nmo), PARAMETER :: ni = (/  31, 59, 90,120,151,181, &
!&                                            212,243,273,304,334,365 /)
!INTEGER,DIMENSION(nmo),PARAMETER :: nii = (/  31, 60, 91,121,152,182, &
!&                                            213,244,274,305,335,366 /)
!REAL   , DIMENSION(4) :: XIM1     ! Weather generator parameter
!REAL    :: td                     ! Time (day)
!REAL    :: th                     ! Solar time (hr)
!REAL    :: tsteph                 ! Timestep for simulation (hr)
!REAL    :: tstep                  ! Timestep for simulation (s)
!REAL    :: dts                    ! Timestep for surface updates (s)
!REAL    :: dtr                    ! Timestep for surface updates (s)
!REAL    :: dtm                    ! Timestep for surface updates (s)
!REAL    :: ntop                   ! N at top of canopy (g/m2)
!REAL    :: lstep                  ! Canopy lai step (m2/m2)
!REAL    :: lc                     ! Cumulalative canopy lai
!                                  ! from top (m2/m2)
!REAL, PARAMETER :: eps = 1.0e-44

!REAL   , DIMENSION(12)  :: pww    ! Probability of wet day
!                                  ! following wet day (fraction)
!REAL   , DIMENSION(12)  :: pwd    ! Probability of wet day
!                                  ! following dry day (fraction)
!REAL   , DIMENSION(12)  :: alpha  ! Rainfall Gamma
!                                  ! distribution shape parameter
!REAL   , DIMENSION(12)  :: beta   ! Rainfall Gamma
!                                  ! distribution scale parameter
!REAL   , DIMENSION(366) :: txm    ! Mean of tmax? (?)
!REAL   , DIMENSION(366) :: txs    ! ?
!REAL   , DIMENSION(366) :: txm1   ! ?
!REAL   , DIMENSION(366) :: txs1   ! ?
!REAL   , DIMENSION(366) :: tnm    ! ?
!REAL   , DIMENSION(366) :: tns    ! ?
!REAL   , DIMENSION(366) :: rmo    ! ?
!REAL   , DIMENSION(366) :: rso    ! ?
!REAL   , DIMENSION(366) :: rm1    ! ?
!REAL   , DIMENSION(366) :: rs1    ! ?
!!INTEGER :: kgen = 1 ! WGEN option
!INTEGER :: KGEN
!REAL   , DIMENSION(366) :: rc     ! Potential radiation (langleys/day)
!INTEGER                 :: idays  ! Generation option
!REAL   , DIMENSION(366) :: vps    ! ?
!REAL   , DIMENSION(366) :: vpm    ! ?
!REAL   , DIMENSION( 12) :: tcfmax ! ?
!REAL   , DIMENSION( 12) :: tcfmin ! ?
!REAL   , DIMENSION( 12) :: rcf    ! ?
!REAL   , DIMENSION( 12) :: cc     ! ?
!REAL   , DIMENSION( 12) :: kay    ! ?
!----------------------------------------------------------------------!
!INTEGER :: KGEN
!INTEGER :: IDAYS
!----------------------------------------------------------------------!
!        DIMENSION TXM(366),TXS(366),TXM1(366),TXS1(366),TNM(366)
!	DIMENSION RMO(366),RSO(366),RM1(366),RS1(366),RC(366),RAIN(366)
     	DIMENSION  SR(12),SSTX(12)
!     	DIMENSION  dTMIN(366),RAD(366),NI(12),SR(12),SSTX(12)
      	DIMENSION SSRAD(32),SRAIN(12),STMAX(12),STMIN(12),SRAD(12)
        DIMENSION TM(12),PW(12),TG(12),RM(12)
!        DIMENSION ALPHAP(12),BETA(12),TM(12),PW(12),TG(12),RM(12)
!       DIMENSION PWD(12),ALPHAP(12),BETA(12),TM(12),PW(12),TG(12),RM(12)
        DIMENSION NWET(12),XNW(12),RGC(12)
!        DIMENSION RCF(12),NWET(12),XNW(12),dTMAX(366),PWW(12),RG(12)
        DIMENSION TAMAX(12),TAMIN(12),SSTN(12)
!        DIMENSION TAMAX(12),TAMIN(12),TNS(366),NII(12),SSTN(12)
        DIMENSION TTMAX(12), TTMIN(12)
!        DIMENSION TTMAX(12), TTMIN(12),TCFMAX(12),TCFMIN(12)
	DIMENSION PTMAX(12),PTMIN(12),PRAD(12),PVP(12),RATIO(366)
	DIMENSION QTMX(12),QTMN(12),QRAD(12),QVP(12)
	DIMENSION SVP(12),SSVP(12)
!	DIMENSION VPM(366),VPS(366),SVP(12),SSVP(12),VAP(366)
	DIMENSION SWS(12),SSWS(12)
!	DIMENSION CC(12),KAY(12),SWS(12),SSWS(12),WINDSP(366)
!	REAL  KAY
	CHARACTER*2 ID
	CHARACTER*5 DATE(12)
	CHARACTER*69 HEADER
	CHARACTER*30 INFILE,FSIMUL,FSUMM
	CHARACTER*3 FRUNIT
	CHARACTER*2 PWUNIT
	CHARACTER*1 MXUNIT,MIUNIT
	CHARACTER*11 SRUNIT
	CHARACTER*7 HYUNIT
	CHARACTER*9 WSUNIT
	COMMON /BLK/ FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT,HYUNIT,WSUNIT
	COMMON /BLK1/ INFILE
	COMMON NFS
!
	DATA FSIMUL,FSUMM  /2*'                              '/
	DATA DATE/' JAN.',' FEB.',' MAR.',' APR.',' MAY ',' JUNE', &
     &    ' JULY',' AUG.',' SEP.',' OCT.',' NOV.',' DEC.'/
!
! FIRST LINE OF EXECUTABLE CODE:
!	WRITE (*,59)
!59	FORMAT(//' **************************************************'
!     &		/'   WEASHO:  SHORT-METHOD WEATHER SIMULATION MODEL'/
!     &		 '   Developed by Shu Geng*, Baoquan Li** '/
!     &           '                and Erich Brandstetter* '//
!     &           '       *   University of California,Davis'/
!     &		 '       ** Agricultural University of South China'//
!     &		 '        Based on Richardson & Wrights Wgen'/
!     &		 '                 and Geng & Auburns Wshort'//           
!     &           '   Please note:  This is a working copy of a model'/
!     &		 '   that is still under development. Please do not '/
!     &		 '   distribute it to others without permission.'/
!     &           '****************************************************'/
!     &           ///' Type <ENTER> to continue....')	
!     	READ (*,105)
!105	FORMAT(A1)
!110	WRITE(*, 120)
!120	FORMAT(/' INPUT FILENAME (=OUTPUT FILE FROM ESTPAR): ',$)
!	READ (*, 130, ERR=110) INFILE
130	FORMAT (A30)
!****GLOBPHOT
	!OPEN (UNIT=1, STATUS='OLD', FILE=INFILE)
!****GLOBPHOT
140	continue
!140	WRITE(*,131)
131 	FORMAT(///////////////)
!	WRITE(*, 141)
141	FORMAT(' OUTPUT FILENAME FOR DAILY SIMULATED DATA: ',$)
!      FSIMUL = 'penn.out'
!	READ (*, 130, ERR=140) FSIMUL
!	OPEN (UNIT=3, STATUS='UNKNOWN', FILE=FSIMUL,  &
!     &  ERR=140)
150	continue
!150	WRITE(*, 151)
151	FORMAT(/' OUTPUT FILENAME FOR MONTHLY TABLES:',$)
!      FSUMM = 'penn.mon'
!	READ (*, 130, ERR=150) FSUMM
!	OPEN (UNIT=4, STATUS='UNKNOWN', FILE=FSUMM,  &
!     &	ERR=150)
!
!**********************************************************************
!*    INPUT # 01 - TITLE    
!*       ACOM(I) - LOCATION IDENTIFIDATION OR OTHER USER               *
!*                 COMMENTS.  80 CHARACTER MAXIMUM                     *
!***********************************************************************
!	WRITE(4,95) FSIMUL
!****GLOBPHOT
	!READ(1,98) HEADER
        HEADER = pHEADER
!****GLOBPHOT
98	FORMAT(A69)
!	WRITE(4,98) HEADER
	ID=HEADER(:2)
!	WRITE(4,99)
99	FORMAT(/' INPUT PARAMETERS:'/)
95	FORMAT('SIMULATION OUTPUT IN FILE: ',A30)
!***********************************************************************
!*     INPUT # 02 - NUMBER OF YEARS, GENERATION CODES, AND LATITUDE
!*           NYRSG - YEARS OF DATA TO BE FENERATED
!*           KGEN - GENERATION OPTION CODE
!*                  IF KGEN = 1,RAIN, MAX TEMP, MIN TEMP, AND
!*                  SOLAR RADIATION WILL BE GENERATED
!*                  IF KGEN = 2 OBSERVED RAIN WILL BE USED AND
!*                  MAX TEMP, MIN TEMP, SOLAR RADIATION WILL
!*                  BE GENERATED
!*           ALAT - STATION LATITUDE IN DEGREES
!*           KTCF - TEMP. CORRECTION FACTOR OPTION CODE
!*                  IF KTCF = 0 NO TEMP CORRECTION WILL BE MADE
!*                  IF KTCF = 2 GENERATED MAX TEMP AND
!*                  MIN TEMP. WILL BE CORRECTED BASED ON
!*                  OBSERVED MEAN MONTHLY MAX AND MIN TEMP
!*                  IF KTCF = 1 GENERATED MAX TEMP AND MIN TEMP
!*                  WILL BE CORRECTED BASED ON OBSERVED MEAN
!*                  MONTHLY TEMP
!*           KRCF - RAIN CORRECTION FACTOR OPTION CODE
!*                  IF KRCF = 1 GENERATED RAIN WILL BE CORRECTED
!*                  BASED ON OBSERVED MEAN MONTHLY RAIN
!*                  IF KRCF = 0 NO RAIN CORRECTION WILL BE MADE
!*******************************************************************
!****GLOBPHOT
	!READ(1,*) NYRSG, KGEN, ALAT, KTCF, KRCF
        NYRSG = pNYR
        KGEN  = pKGEN
        ALAT  = pALAT
        KTCF  = pKTCF
        KRCF  = pKRCF
!****GLOBPHOT
!***** CALCULATE MAXIMUM SOLAR RADIATION FOR EACH DAY
      XLAT = ALAT*6.2832/360.
      DO 6 I = 1,366
      XI = I
      SD = ASIN (0.39795 * COS (0.98563 * (XI - 173.0) * 6.2832/360.))
!      SD = 0.4102*SIN(0.0172*(XI-80.25))
      CH = -TAN(XLAT)*TAN(SD)
      IF(CH .GT. 1.0) H = 0.
      IF(CH .GT. 1.0)GO TO 5
      IF( CH .LT. -1.0)H=3.1416
      IF(CH .LT. -1.0) GO TO 5
      H = ACOS(CH)
 5    DD = 1.0+0.0335*SIN(0.0172*(XI+88.2))
!  Made same as 'powerofsun' book. DD allows for variation in earth-sun
!  distance.
! 5    DD = 1.0 + 0.034 * COS (360.0 * XI / 365.25)
!  Langleys/day?
!      dd=1.0
      RC(I)=897.02*DD*((H*SIN(XLAT)*SIN(SD))+(COS(XLAT)*COS(SD)*SIN(H)&
     &))
!----------------------------------------------------------------------!
!  Do not understand why this is here. Commented out as does not match
!  with my calculated potential, hence ratio does not apply to each
!  point.
!      RC(I) = RC(I) * 0.8
!----------------------------------------------------------------------!
 6    CONTINUE
      DO 7 I = 1,12
      TTMAX(I)=0.
      TTMIN(I)=0.
      RM(I) = 0.
 7    CONTINUE
      IF(KGEN .EQ. 2) GO TO 10
!	INPUTS #03,04,05,06 ARE RAINFALL PARAMETERS; OMIT IF KGEN=2
!*****************************************************************************
!      	INPUT # 3  WEIBILL DISTRIBUTION PARAMETER CC(I) (ADDED BY BQL)
!	INPUT # 4  WEIBILL DISTRIBUTION PARAMETER KAY(I) (BQL)
!	INPUT # 5  PROBABILITY OF WET GIVEN WET PWW(I)
!	INPUT # 6  PROBABILITY OF WET GIVEN DRY PWD(I)
!	INPUT # 7  GAMMA DISTRIBUTION SHAPE PARAMETER ALPHAP(I)
!	INPUT # 8  GAMMA DISTRIBUTION SCALE PARAMETER BETA(I)
!	ALL INPUTS CONTAIN 12 MONTHLY VALUES
!****************************************************************************
      IF((NFS.EQ.1).OR.(NFS.EQ.3)) GOTO 1000
!****GLOBPHOT
      !READ(1,103) (CC(I),I=1,12)
      !READ(1,103) (KAY(I),I=1,12)
1000  CONTINUE	
!1000  READ(1,103) (PWW(I),I=1,12)
      !READ(1,103) (PWD(I),I=1,12)
      !READ(1,103) (ALPHAP(I),I=1,12)
      !READ(1,103) (BETA(I),I=1,12)
      CC     (:) = pCC     (:)
      KAY    (:) = pKAY    (:)
      PWW    (:) = pPWW    (:)
      PWD    (:) = pPWD    (:)
      ALPHAP (:) = pALPHAP (:)
      BETA   (:) = pBETA   (:)
 103  FORMAT(12F6.2)
!****GLOBPHOT
!************************************************************************
!*    INPUT # 9 - FOURIER COEFFICIENTS OF MAX TEMP ON DRY DAYS
!*                TXMD  - MEAN OF TMAX - DRY
!*                ATX   - AMPLITUDE OF TMAX - WET OR DRY
!*                CVTX  - MEAN OF COEF. OF VAR. OF TMAX - WET OR DRY
!*                ACVTX - AMPLITUDE OF COEF.OF VAR. OF TMAX - WET OR DRY
!*		  PEAKTX - PERIOD OF TMAX (ADDED BY J. AUBURN) WET OR DRY
!*    INPUT # 10 - FOURIER COEFFICIENTS OF MAX TEMP ON WET DAYS
!*                TXMW  - MEAN OF TMAX - WET
!*    INPUT # 11 - FOURIER COEFFICIENTS OF MIN TEMP
!*        	  TN    - MEAN OF  TMIN - WET OR DRY
!*        	  ATN   - AMPLITUDE OF  TMIN - WET OR DRY
!*       	  CVTN  - MEAN OF COEF. OF VAR. OF TMIN - WET OR DRY
!*     	          ACVTN - AMPLITUDE OF COEF. OF VAR. OF TMIN - WET OR DRY
!*		  PEAKTN - PERIOD OF TMIN (ADDED BY J. AUBURN) WET OR DRY
!*				(NOT YET USED)
!*    INPUT # 12 - FOURIER COEFFICIENTS OF RAD ON DRY DAYS
!*     	          RMD   - MEAN OF RAD - DRY
!*   	          AR    - AMPLITUDE OF RAD - WET OR DRY
!*		  PEAKSR - PERIOD OF SOL RAD,WET OR DRY (ADDED BY J. AUBURN)
!*		  CVRD  - MEAN OF CV OF SOL RAD, DRY (ADDED BY JSA)
!*		  ACVRD - AMPLITUDE OF SOL RAD, DRY (JSA)
!*    INPUT # 13 - FOURIER COEFFICIENTS OF RAD ON WET DAYS
!*       	   RMW  - MEAN OF RAD - WET
!*		   CVRW - MEAN OF CV OF SOL RAD, WET (JSA)
!*		   ACVRW -AMPLITUDE OF CV OF SOL RAD, WET (JSA)
!************************************************************************
!***GLOBPHOT
 10     CONTINUE
! 10     READ(1,'(F7.1,F7.2,2F7.2,F7.1)') TXMD,ATX,CVTX,ACVTX,PEAKTX
	!READ(1,'(F7.1)') TXMW
	!READ(1,'(F7.1,F7.2,2F7.2,F7.1)') TN,ATN,CVTN,ACVTN,PEAKTN
	!READ(1,'(2F7.2,F7.1,2F7.2)') RMD,AR,PEAKSR,CVRD,ACVRD
	!READ(1,'(F7.2,2F7.2)') RMW,CVRW,ACVRW
	TXMD   = pTXMD
	ATX    = pATX
	CVTX   = pCVTX
	ACVTX  = pACVTX
	PEAKTX = pPEAKTX
	TXMW   = pTXMW
	TN     = pTN
	ATN    = pATN
	CVTN   = pCVTN
	ACVTN  = pACVTN
	PEAKTN = pPEAKTN
	RMD    = pRMD
	AR     = pAR
	PEAKSR = pPEAKSR
	CVRD   = pCVRD
	ACVRD  = pACVRD
	RMW    = pRMW
	CVRW   = pCVRW
	ACVRW  = pACVRW
!***GLOBPHOT
!*************************************************************************
!*    INPUT # 14 - FOURIER COEFFICIENTS OF HUMIDITY
!*      	   VP  - MEAN OF HUMIDITY (ADDED BY BAOQUAN LI)
!*		   AVP - AMPLITUDE OF HUMIDITY (BQL)
!*		   PEAKVP - PERIOD OF HUMIDITY (BQL)
!*		   CVVP - MEAN OF CV OF HUMIDITY (BQL)
!*		   ACVVP - AMPLITUDE OF CV OF HUMIDITY (BQL)
!*************************************************************************		
	IF((NFS.EQ.1).OR.(NFS.EQ.2)) GOTO 1005
!****GLOBPHOT
	!READ(1,'(2F7.2,F7.1,2F7.2)') VP,AVP,PEAKVP,CVVP,ACVVP
	VP     = ppVP
	AVP    = pAVP
	PEAKVP = pPEAKVP
	CVVP   = pCVVP
	ACVVP  = pACVVP
!****GLOBPHOT
 101    FORMAT(5(F7.0))
1005	continue
!1005	WRITE(4,801) TXMD,ATX,CVTX,ACVTX,PEAKTX
 801	FORMAT('TXMD  =',F10.3,/,'ATX   =',F10.3,/,'CVTX  =',F10.3,/, &
     &  'ACVTX =',F10.3/'PEAKTX=',F10.3)
!	WRITE(4,802) TXMW
 802	FORMAT('TXMW  =',F10.3)
!        WRITE(4,803) TN,ATN,CVTN,ACVTN,PEAKTN
 803	FORMAT('TN    =',F10.3,/,'ATN   =',F10.3,/,'CVTN  =',F10.3,/, &
     &	'ACVTN =',F10.3,/,'PEAKTN=',F10.3)
!	WRITE(4,804) RMD,AR,PEAKSR,CVRD,ACVRD
 804	FORMAT('RMD   =',F10.3,/,'AR    =',F10.3,/'PEAKSR=',F10.3, &
     &  /'CVRD  =',F10.3,/'ACVRD =',F10.3)
!	WRITE(4,805) RMW,CVRW,ACVRW
805     FORMAT('RMW   =',F10.3,/'CVRW  =',F10.3/'ACVRW =',F10.3)
	IF((NFS.EQ.1).OR.(NFS.EQ.2)) GOTO 1010
!	WRITE(4,806) VP,AVP,CVVP,ACVVP,PEAKVP
806	FORMAT('VP    =',F10.3,/'AVP   =',F10.3,/'CVVP  =',F10.3, &
     &	/,'ACVVP =',F10.3,/'PEAKVP=',F10.3)
!
!   CONVERT SOLRAD FROM 10**6 J/M**2/DAY TO LANGLEYS (ORIG SIMUL IN LANGLEYS)
!   SO WILL STAY IN THOSE UNITS, AND CONVERT BACK TO STD FOR OUTPUT
!
1010	RMD = RMD * 23.87
	RMW = RMW * 23.87
	AR = AR * 23.87
!
      D1 = TXMD - TXMW
      D2 = RMD - RMW
      IF(KTCF  .EQ. 0) GO TO 12
      IF(KTCF  .EQ. 2) GO TO 8
!************************************************************************
!*    INPUT # 15 - MONTHLY VALUES OF ACTUAL MEAN TEMP
!*                 OMIT IF KTCF = 1 OR 2
!*        TM(I) - 12 MONTHLY VALUES OF ACTUAL MEAN TEMP
!*************************************************************************
      READ(1,103)(TM(I),I=1,12)
      GO TO 12
!*************************************************************************
!*    INPUT # 16 - MONTHLY VALUES OF ACTUAL MEAN MAX TEMP
!*                 OMIT IF KTCF = 0 OR 1
!*        TTMAX(I) - 12 MONTHLY VALUES OF ACTUAL MEAN MAX TEMP
!*************************************************************************
 8    READ(1,103) (TTMAX(I),I=1,12)
!*************************************************************************
!*    INPUT # 17 - MONTHLY VALUES OF ACTUAL MEAN MIN TEMP
!*                 OMIT IF KTCF = 0 OR 1
!*        TTMIN(I) - 12 MONTHLY VALUES OF ACTUAL MEAN MIN TEMP
!**************************************************************************
      READ(1,103) (TTMIN(I),I=1,12)
 12   IF(KRCF .EQ. 0) GO TO 13
!**************************************************************************
!*    INPUT # 18 - MONTHLY VALUES OF ACTUAL MEAN RAINFALL
!*                 OMIT IF KRCF = 0
!*        RM(I) = 12 MONTHLY VALUES OF ACTUAL MEAN RAINFALL
!***********************************************************************
      READ(1,103)(RM(I),I=1,12)
	CLOSE(1)
 13	CONTINUE
      IF((NFS.EQ.1).OR.(NFS.EQ.3)) GOTO 1015
!      WRITE(4,698)(CC(I),I=1,12)
 698  FORMAT('CC     ',12F6.3)
!      WRITE(4,699)(KAY(I),I=1,12)
 699  FORMAT('KAY    ',12F6.3)
1015  continue
!1015  WRITE(4,701)(PWW(I),I=1,12)
 701  FORMAT('P(W/W) ',12F6.3)
!      WRITE(4,702) (PWD(I),I=1,12)
 702  FORMAT('P(W/D) ',12F6.3)
!      WRITE(4,703) (ALPHAP(I),I=1,12)
 703  FORMAT('ALPHAP  ',12F6.3)
!      WRITE(4,704) (BETA(I),I=1,12)
 704  FORMAT('BETA   ',12F6.2)
      DO 11 J = 1,366
      XJ = J
      dt_wg = COS(.0172*(XJ-PEAKTX))
      DR = COS(.0172*(XJ-PEAKSR))
      DV = COS(.0172*(XJ-PEAKVP))
      TXM(J) = TXMD+ATX*dt_wg
      XCR1=CVTX+ACVTX*dt_wg
      IF(XCR1 .LT. 0.0) XCR1=0.06
      TXS(J)=TXM(J)*XCR1
      TXM1(J) = TXM(J) - D1
      TXS1(J)=TXM1(J)*XCR1
      TNM(J) = TN + ATN*dt_wg
      XCR2=CVTN+ACVTN*dt_wg
      IF(XCR2 .LT. 0.0) XCR2=0.06
      TNS(J)=TNM(J)*XCR2
      RMO(J) = RMD+AR * DR
      XCR3=CVRD+ACVRD*DR
      IF(XCR3 .LT. 0.0) XCR3=0.06
      RSO(J)=RMO(J)*XCR3
      RM1(J) = RMO(J) - D2
      XCR4=CVRW+ACVRW*DR
      IF(XCR4 .LT. 0.0) XCR4=0.06
      RS1(J)=RM1(J)*XCR4
      VPM(J)=VP+AVP*DV
      XCR5=CVVP+ACVVP*DV
      IF(XCR5 .LT. 0.0) XCR5=0.06
      VPS(J)=VPM(J)*XCR5
 11   CONTINUE
      DO 22 im_wg=1,12
      XNW(im_wg) = 0.
      SR(im_wg) = 0.
      SSTX(im_wg) = 0.
      SSTN(im_wg) = 0.
      SSRAD(im_wg) = 0.
      SSVP(im_wg) = 0.
      SSWS(im_wg) = 0.
      TCFMAX(im_wg) = 0.
      TCFMIN(im_wg) = 0.
      RCF(im_wg) = 1.0
      PW(im_wg) = PWD(im_wg)/(1. -PWW(im_wg)+PWD(im_wg) + EPS)
      S1 = 0.
      S2 = 0.
      S3 = 0.
      NL = NI(im_wg)
      IF(im_wg .EQ. 1) GO TO 14
      NF = NI(im_wg-1) + 1
      GO TO 15
 14   NF = 1
 15   CONTINUE
      ZN = NL - NF + 1
      DO 16 J = NF,NL
      S1 = S1 + TXM(J)/(ZN + EPS)
      S2 = S2 + TXM1(J)/(ZN + EPS)
      S3 = S3 + TNM(J)/(ZN + EPS)
 16   CONTINUE
!*****CALCULATE MONTHLY RAINFALL CORRECTION FACTOR
      RGC(im_wg) = ALPHAP(im_wg)*BETA(im_wg)*ZN*PW(im_wg)
      IF(KRCF .EQ. 0) GO TO 17
      RCF(im_wg) = RM(im_wg)/(RGC(im_wg) + EPS)
 17   IF(KTCF .EQ. 0) GO TO 22
!*****CALCULATE MONTHLY TEMP CORRECTION FACTOR
      IF(KTCF .EQ. 2) GO TO 18
      TMD = (S1 + S3) / 2.
      TMW = (S2 + S3) / 2.
      TG(im_wg) = TMW*PW(im_wg)+TMD*(1-PW(im_wg))
      TCFMAX(im_wg) = TM(im_wg) - TG(im_wg)
      TCFMIN(im_wg) = TCFMAX(im_wg)
      GO TO 22
 18   TAMAX(im_wg)=S2*PW(im_wg) + S1*(1. -PW(im_wg))
      TAMIN(im_wg)=S3
      IF(KTCF .EQ. 0.) GO TO 22
      TCFMAX(im_wg)=TTMAX(im_wg)-TAMAX(im_wg)
      TCFMIN(im_wg)=TTMIN(im_wg)-TAMIN(im_wg)
 22   CONTINUE
      IF(KRCF .EQ. 0) GO TO 52
 !     WRITE(4,712)(RM(I),I=1,12)
 712  FORMAT(10X,'ACT MEAN RAIN',12F7.2)
 !     WRITE(4,713) (RGC(I),I=1,12)
 713  FORMAT(10X,'EST MEAN RAIN',12F7.2)
 !     WRITE(4,714)(RCF(I),I=1,12)
 714  FORMAT(10X,'RAIN CF      ',12F7.3)
 52   IF (KTCF .EQ. 0) GO TO 19
      IF(KTCF .EQ. 2) GO TO 51
 !     WRITE(4,708) (TM(I),I=1,12)
 708  FORMAT(10X,'ACT MEAN TEMP',12F7.1)
 !     WRITE(4,711)(TG(I),I=1,12)
 711  FORMAT(10X,'EST MEAN TEMP',12F7.1)
      GO TO 50
 51   continue
! 51   WRITE(4,722) (TTMAX(I),I=1,12)
 722  FORMAT(10X,'ACT MEAN TMAX',12F7.1)
!      WRITE(4,723) (TTMIN(I),I=1,12)
 723  FORMAT(10X,'ACT MEAN TMIN',12F7.1)
!      WRITE(4,720) (TAMAX(I),I=1,12)
 720  FORMAT(10X,'EST MEAN TMAX',12F7.1)
!      WRITE(4,721) (TAMIN(I),I=1,12)
 721  FORMAT(10X,'EST MEAN TMIN',12F7.1)
 50   continue
! 50   WRITE(4,709) (TCFMAX(I),I=1,12)
 709  FORMAT(10X,'CF. MEAN TMAX',12F7.1)  
!      WRITE(4,724) (TCFMIN(I),I=1,12)
 724  FORMAT(10X,'CF. MEAN TMIN',12F7.1)
 19   XYR = NYRSG
      SYTX = 0.
      SYTN = 0.
      SYRAD = 0.
      SYR = 0.
      SYNW = 0.
      SYVP = 0.
      SYWS = 0.
      DO 40 I = 1,NYRSG
!      WRITE(*,799) I
799   FORMAT(/' SIMULATING YEAR NO. ',I5)
      IYR = I
      IF(KGEN .EQ. 1) GO TO 20
      KK = 0
      IJ = 1
!************************************************************
!*    INPUT # 19 - MEASURED RAINFALL FOR NYRSG
!*                 OMIT IF KGEN = 1
!*        RAIN(I) - ACTUAL RAINFALL DATA - ONE VALUE PER DAY
!*                  FOR NYRSG
!************************************************************
!	WRITE(*,899)
899	FORMAT (' USING ACTUAL RAINFALL; ENTER FILENAME FOR DATA')
!818	WRITE(*, 120)
!	READ (*, 130, ERR=818) INFILE
	OPEN (UNIT=1, STATUS='OLD', FILE=INFILE)
!
!
 21   READ(1,102) IYR,MO,IDAY,RAIN(IJ)
 102  FORMAT(4X,3I2,20X,F10.0)
      IF(KK .EQ. 1) GO TO 24
 20   IDAYS = 365
      IFLG = MOD(IYR,4)
      IF(IFLG .EQ. 0) IDAYS = 366
      KK = 1
      IF(KGEN .EQ. 1) GO TO 28

 24   IJ = IJ + 1
      IF(IJ .LE. IDAYS) GO TO 21
 28   CONTINUE
	CLOSE(UNIT=1)
!      CALL WGEN(PWW,PWD,ALPHAP,BETA,TXM,TXS,TXM1,TXS1,TNM,TNS,RMO,RSO,&
!     &   RM1,RS1,RAIN, dTMAX, dTMIN,RAD,KGEN,RC,IDAYS,NI,NII,VPS,VPM,VAP,&
!     &   TCFMAX,TCFMIN,RCF,CC,KAY,WINDSP)
      DO 23 im_wg = 1,12
      SRAIN(im_wg) = 0.
      STMAX(im_wg) = 0.
      STMIN(im_wg) = 0.
      SRAD(im_wg) = 0.
      NWET(im_wg) = 0.
      SVP(im_wg) = 0.	
      SWS(im_wg) = 0.	
23    CONTINUE
      im_wg = 1
      YTMAX = 0.
      YTMIN = 0.
      YRAD = 0.
      YVP = 0.
      YWS = 0.
      RYR = 0.
      NYWET = 0.
      IDA = 0
      DO 30 J=1,IDAYS
      IDA = IDA + 1
      IF(IDAYS .EQ. 366) GO TO 27
      IF(J .GT.NI(im_wg)) GO TO 251
      GO TO 29
 251  im_wg = im_wg + 1
      IDA = 1
      GO TO 29
 27   IF(J .GT. NII(im_wg)) GO TO 251
 29   CONTINUE
	  RATIO(J) = RAD(J)/(RC(J) + EPS)
!  WRITE OUT RESULTS
!  
!  CONVERT TO OUTPUT UNITS BEFORE PRINTING...
	IF (PWUNIT.EQ.'IN') THEN
	TRAIN=.03937*RAIN(J)
	ELSE
	TRAIN=RAIN(J)
	ENDIF
	IF (MXUNIT.EQ.'C') THEN
	TMAXPR = (dTMAX(J) - 32.) * 5./9.
	ELSE
	TMAXPR = dTMAX(J)
	ENDIF
	IF (MIUNIT.EQ.'C') THEN
	TMINPR = ( dTMIN(J) - 32.) * 5./9.
	ELSE
	TMINPR = dTMIN(J)
	ENDIF
	IF (SRUNIT.EQ.'MJ/M*M/DAY ') THEN
	RADPR = RAD(J) * 0.04189
	ELSE
	RADPR =RAD(J)
	ENDIF
	IF (HYUNIT.EQ.'C(DEWP)')THEN
	VAPPR=237.3*ALOG(VAP(J)/6.108)/17.27/ &
     &        (1.-ALOG(VAP(J)/6.108)/17.27)
	ELSEIF (HYUNIT.EQ.'F(DEWP)') THEN
	VAPPR=237.3*ALOG(VAP(J)/6.108)/17.27/ &
     &        (1.-ALOG(VAP(J)/6.108)/17.27)
	VAPPR=VAPPR*9./5.+32.
	ELSEIF (HYUNIT.EQ.'MPa(VP)') THEN
	VAPPR = VAP(J)
	ELSE
	VAPPR = 0.
	ENDIF
	IF (WSUNIT.EQ.'MILE/HOUR') THEN
	WINDSP(J)=2.237*WINDSP(J)
	ENDIF
	IF(NFS.EQ.1) THEN
	ASSIGN 1001 TO NFST
	ELSEIF(NFS.EQ.2) THEN
	ASSIGN 1002 TO NFST
	ELSEIF(NFS.EQ.3) THEN
	ASSIGN 1003 TO NFST
	ELSE
	ASSIGN 1004 TO NFST
	ENDIF
	GOTO NFST (1001,1002,1003,1004)
1001	continue
!1001	WRITE(3,1200) ID,IYR,im_wg,IDA,J,TRAIN,TMAXPR,TMINPR,RADPR, &
!     &  RATIO(J)
1200	FORMAT(A2,3I2,I3,3F5.1,2F6.2)
	GOTO 800
1002	continue
!1002	WRITE(3,1205) ID,IYR,im_wg,IDA,J,TRAIN,TMAXPR,TMINPR,RADPR, &
!     &  WINDSP(J),RATIO(J)
1205	FORMAT(A2,3I2,I3,3F5.1,3F6.2)
	GOTO 800
1003	continue
!1003	WRITE(3,1205) ID,IYR,im_wg,IDA,J,TRAIN,TMAXPR,TMINPR,RADPR, &
!     &  VAPPR,RATIO(J)
	GOTO 800
1004    continue
!1004    WRITE(3,1210) ID,IYR,im_wg,IDA,J,TRAIN,TMAXPR,TMINPR,RADPR, &
!     &  VAPPR,WINDSP(J),RATIO(J)
1210	FORMAT(A2,3I2,I3,3F5.1,4F6.2)
800   CONTINUE
 25   CONTINUE
      IF(RAIN(J) .LT. 0.005) GO TO 26
      NWET(im_wg) = NWET(im_wg) + 1
      NYWET = NYWET + 1
 26   CONTINUE
      SRAIN(im_wg) = SRAIN(im_wg) + RAIN(J)
      STMAX(im_wg) = STMAX(im_wg) + dTMAX(J)
      STMIN(im_wg) = STMIN(im_wg) +  dTMIN(J)
      SRAD(im_wg) = SRAD(im_wg) + RAD(J)
      SVP(im_wg) = SVP(im_wg) + VAP(J)
      SWS(im_wg) = SWS(im_wg) + WINDSP(J)
      RYR = RYR + RAIN(J)
 30   CONTINUE
      XNM1 = 0.
      DO 35 im_wg = 1,12
      XXN = NI(im_wg)
      XNI = XXN - XNM1
      XNM1 = XXN
      ANW = NWET(im_wg)
      XNW(im_wg) = XNW(im_wg) + ANW/(XYR + EPS)
      SR(im_wg) = SR(im_wg) + SRAIN(im_wg) / (XYR + EPS)
      STMAX(im_wg) = STMAX(im_wg)/ (XNI + EPS)
      SSTX(im_wg) = SSTX(im_wg) + STMAX(im_wg) / (XYR + EPS)
      STMIN(im_wg) = STMIN(im_wg) / (XNI + EPS)
      SSTN(im_wg) = SSTN(im_wg) + STMIN(im_wg) / (XYR + EPS)
      SRAD(im_wg) = SRAD(im_wg) / (XNI + EPS)
      SSRAD(im_wg) = SSRAD(im_wg) + SRAD(im_wg) / (XYR + EPS)
      SVP(im_wg) = SVP(im_wg) / (XNI + EPS)
      SSVP(im_wg) = SSVP(im_wg) + SVP(im_wg) / (XYR + EPS)
      SWS(im_wg) = SWS(im_wg) / (XNI + EPS)
      SSWS(im_wg) = SSWS(im_wg) + SWS(im_wg) / (XYR + EPS)
      YTMAX = YTMAX + STMAX(im_wg) / 12.
      YTMIN = YTMIN + STMIN(im_wg) / 12.
      YVP = YVP + SVP(im_wg) / 12.	
      YWS = YWS + SWS(im_wg) / 12.
      YRAD = YRAD + SRAD(im_wg) / 12.
 35   CONTINUE
      SYTX= SYTX + YTMAX/(XYR + EPS)
      SYTN = SYTN + YTMIN/(XYR + EPS)
      SYRAD = SYRAD + YRAD / (XYR + EPS)
      SYR = SYR + RYR / (XYR + EPS)
      XYNW = NYWET
      SYNW= SYNW + XYNW / (XYR + EPS)
      SYVP = SYVP + YVP / (XYR + EPS)
      SYWS = SYWS + YWS / (XYR + EPS)
!      WRITE(4,201) IYR
 201  FORMAT(//,5X,'SUMMARY OF THE SIMULATION FOR YEAR ',I5,/)
 	IF(NFS.EQ.1) THEN
	ASSIGN 1301 TO NFST
	ELSEIF(NFS.EQ.2) THEN
	ASSIGN 1302 TO NFST
	ELSEIF(NFS.EQ.3) THEN
	ASSIGN 1303 TO NFST
	ELSE
	ASSIGN 1304 TO NFST
	ENDIF
	GOTO NFST (1301,1302,1303,1304)	
 1301	continue
! 1301	WRITE(4,1300) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT
 1300	FORMAT(' MONTH',2X,'WET-DAYS',1X,'RAINFAL',2X,'MAX-TEMP',2X, &
     &  'MIN-TEMP',3X,'SOLAR-RAD'/,8X, &
     &  ' (',A3,')  ',3X,'(',A2,')',6X,'(',A1,')',7X,'(',A1,')', &
     &  4X,'(',A11,')')
	GOTO 1350
 1302	continue
! 1302	WRITE(4,1305) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT,WSUNIT
 1305	FORMAT(' MONTH',1X,'WET-DAYS',2X,'RAINFAL',2X,'MAX-TEMP',2X, &
     &  'MIN-TEMP',3X,'SOLAR-RAD',4X,'WIND-SPD'/,8X, &
     &  '(',A3,')  ',3X,'(',A2,')',6X,'(',A1,')',7X,'(',A1,')', &
     &  4X,'(',A10,')',1X,'(',A9,')')
	GOTO 1350
 1303	continue
! 1303	WRITE(4,1310) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT,HYUNIT
 1310	FORMAT(' MONTH',1X,'WET-DAYS',2X,'RAINFAL',2X,'MAX-TEMP',2X, &
     &  'MIN-TEMP',3X,'SOLAR-RAD',4X,'HUMIDITY',/,8X, &
     &  '(',A3,')  ',3X,'(',A2,')',6X,'(',A1,')',7X,'(',A1,')', &
     &  4X,'(',A10,')',2X,'(',A7,')')
	GOTO 1350
 1304	continue
! 1304	WRITE(4,1315) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT,HYUNIT,WSUNIT
 1315	FORMAT(' MONTH',1X,'WET-DAYS',2X,'RAINFAL',2X,'MAX-TEMP',1X, &
     &  'MIN-TEMP',2X,'SOLAR-RAD',3X,'HUMIDITY',2X,'WIND-SPD'/,8X, &
     &  '(',A3,')  ',3X,'(',A2,')',6X,'(',A1,')',7X,'(',A1,')', &
     &  2X,'(',A10,')',1X,'(',A7,')',1X,'(',A9,')')
!
! CONVERT TO OUTPUT UNITS...
 1350     DO 204 im_wg=1,12
	IF (PWUNIT.EQ.'IN') THEN
	SRAIN(im_wg)=SRAIN(im_wg)*0.03937
	ENDIF
	IF (MXUNIT.EQ.'C') THEN
	PTMAX(im_wg) = (STMAX(im_wg) - 32.) * 5./9.
	ELSE
	PTMAX(im_wg)=STMAX(im_wg)
	ENDIF
	IF (MIUNIT.EQ.'C') THEN
	PTMIN(im_wg) = (STMIN(im_wg) - 32.) * 5./9.
	ELSE
	PTMIN(im_wg) = STMIN(im_wg)
	ENDIF
	IF (SRUNIT.EQ.'MJ/M*M/DAY ') THEN
	PRAD(im_wg) = SRAD(im_wg) * 0.04189
	ELSE
	PRAD(im_wg) = SRAD(im_wg)
	ENDIF
	IF (HYUNIT.EQ.'C(DEWP)')THEN
	PVP(im_wg)=237.3*ALOG(SVP(im_wg)/6.108)/17.27/ &
     &        (1.-ALOG(SVP(im_wg)/6.108)/17.27)
	ELSEIF (HYUNIT.EQ.'F(DEWP)') THEN
	PVP(im_wg)=237.3*ALOG(SVP(im_wg)/6.108)/17.27/ &
     &        (1.-ALOG(SVP(im_wg)/6.108)/17.27)
	PVP(im_wg)=PVP(im_wg)*9./5.+32.
	ELSEIF (HYUNIT.EQ.'MPa(VP)') THEN
	PVP(im_wg)=SVP(im_wg) 
	ELSE
	PVP(im_wg)=0.
	ENDIF
	IF(NFS.EQ.1) THEN
	ASSIGN 1401 TO NFST
	ELSEIF(NFS.EQ.2) THEN
	ASSIGN 1402 TO NFST
	ELSEIF(NFS.EQ.3) THEN
	ASSIGN 1403 TO NFST
	ELSE
	ASSIGN 1404 TO NFST
	ENDIF
	GOTO NFST (1401,1402,1403,1404)
1401	continue
!1401	WRITE(4,1405) DATE(im_wg),NWET(im_wg),SRAIN(im_wg),&
!       &PTMAX(im_wg),PTMIN(im_wg),PRAD(im_wg)
1405	FORMAT(A5,4X,I3,4X,F7.2,4X,F5.1,5X,F5.1,6X,F6.1)	
	GOTO 204
1402	continue
!1402	WRITE(4,1410) DATE(im_wg),NWET(im_wg),SRAIN(im_wg),&
!       &PTMAX(im_wg),PTMIN(im_wg),PRAD(im_wg),SWS(im_wg)
	GOTO 204
1403	continue
!1403	WRITE(4,1410) DATE(im_wg),NWET(im_wg),SRAIN(im_wg),&
!       &PTMAX(im_wg),PTMIN(im_wg),PRAD(im_wg),PVP(im_wg)
1410	FORMAT(A5,4X,I3,4X,F7.2,4X,F5.1,5X,F5.1,7X,F6.1,6X,F6.2)
	GOTO 204
1404	continue
!1404	WRITE(4,1415) DATE(im_wg),NWET(im_wg),SRAIN(im_wg),&
!       &PTMAX(im_wg),PTMIN(im_wg),PRAD(im_wg),PVP(im_wg),SWS(im_wg)
1415	FORMAT(A4,4X,I3,5X,F7.2,3X,F5.1,5X,F5.1,4X,F6.1,5X,F6.2,5X,F6.2)
 204	CONTINUE
	IF (PWUNIT.EQ.'IN') THEN
	RYR=RYR*0.03937
	ENDIF
	IF(NFS.EQ.1) THEN
	ASSIGN 1501 TO NFST
	ELSEIF(NFS.EQ.2) THEN
	ASSIGN 1502 TO NFST
	ELSEIF(NFS.EQ.3) THEN
	ASSIGN 1503 TO NFST
	ELSE
	ASSIGN 1504 TO NFST
	ENDIF
	GOTO NFST (1501,1502,1503,1504)
1501	continue
!1501	WRITE(4,1515) NYWET,RYR
1515	FORMAT(/' TOTAL',2X,I3,5X,F7.2,7X,'__',8X,'__',10X,'__')
	GOTO 1550
1502	continue
!1502	WRITE(4,1520) NYWET,RYR    
        GOTO 1550
1503	continue
!1503	WRITE(4,1520) NYWET,RYR
1520   FORMAT(/' TOTAL',3X,I3,4X,F7.2,7X,'__',8X,'__',11X,'__',10X,'__')
	GOTO 1550
1504	continue
!1504	WRITE(4,1525) NYWET,RYR
1525   FORMAT(/' TOTAL',2X,I3,5X,F7.2,6X,'__',8X,'__',8X,'__',9X,'__', &
     &  9X,'__')
1550	IF (MXUNIT.EQ.'C') THEN
	PYTMAX = (YTMAX - 32.) * 5./9.
	ELSE
	PYTMAX = YTMAX
	ENDIF
	IF (MIUNIT.EQ.'C') THEN
	PYTMIN = (YTMIN - 32.) * 5./9.
	ELSE
	PYTMIN = YTMIN
	ENDIF
	IF (SRUNIT.EQ.'MJ/M*M/DAY ') THEN
	PYRAD = YRAD * 0.04189
	ELSE
	PYRAD = YRAD
	ENDIF
	IF (HYUNIT.EQ.'C(DEWP)')THEN
	PYVP=237.3*ALOG(YVP/6.108)/17.27/ &
     &        (1.-ALOG(YVP/6.108)/17.27)
	ELSEIF (HYUNIT.EQ.'F(DEWP)') THEN
	PYVP=237.3*ALOG(YVP/6.108)/17.27/ &
     &        (1.-ALOG(YVP/6.108)/17.27)
	PYVP=PYVP*9./5.+32.
	ELSEIF (HYUNIT.EQ.'MPa(VP)') THEN
	PYVP = YVP	
	ELSE
	PYVP = 0.
	ENDIF	
	IF(NFS.EQ.1) THEN
	ASSIGN 1601 TO NFST
	ELSEIF(NFS.EQ.2) THEN
	ASSIGN 1602 TO NFST
	ELSEIF(NFS.EQ.3) THEN
	ASSIGN 1603 TO NFST
	ELSE
	ASSIGN 1604 TO NFST
	ENDIF
	GOTO NFST (1601,1602,1603,1604)
1601	continue
!1601	WRITE(4,1615) PYTMAX,PYTMIN,PYRAD
	GOTO 40
1615	FORMAT(' AVER.',3X,'__',10X,'__',4X,F5.1,5X,F5.1,6X,F6.1)
1602	continue
!1602	WRITE(4,1620) PYTMAX,PYTMIN,PYRAD,YWS
	GOTO 40
1603	continue
!1603	WRITE(4,1620) PYTMAX,PYTMIN,PYRAD,PYVP
 1620	FORMAT(' AVER.',3X,'__',10X,'__',4X,F5.1,5X,F5.1,7X,F6.1,6X, &
     &   F6.2)
	GOTO 40
1604	continue
!1604	WRITE(4,350) PYTMAX,PYTMIN,PYRAD,PYVP,YWS
 350	FORMAT(' AVER.',3X,'__',10X,'__',3X,F5.1,5X,F5.1,4X,F6.1,5X, &
     &   F6.2,5X,F6.2)   
  40   CONTINUE
!        WRITE(4,206) NYRS
 206    FORMAT(///,5X,'SUMMARY FOR',I5,' YEARS')
!	WRITE(*,133) I-1
 133	FORMAT(////' * OVERALL SUMMARY OF ',I3,' YEAR(S) :')	
	IF(NFS.EQ.1) THEN
	ASSIGN 1701 TO NFST
	ELSEIF(NFS.EQ.2) THEN
	ASSIGN 1702 TO NFST
	ELSEIF(NFS.EQ.3) THEN
	ASSIGN 1703 TO NFST
	ELSE
	ASSIGN 1704 TO NFST
	ENDIF
	GOTO NFST (1701,1702,1703,1704)	
1701	continue
!1701	WRITE(4,1300) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT
!	WRITE(*,1300) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT	
	GOTO 1750
1702	continue
!1702	WRITE(4,1305) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT,WSUNIT
!	WRITE(*,1305) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT,WSUNIT
	GOTO 1750
1703	continue
!1703	WRITE(4,1310) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT,HYUNIT
!	WRITE(*,1310) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT,HYUNIT
	GOTO 1750
1704	continue
!1704	WRITE(4,1315) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT,HYUNIT,WSUNIT
!	WRITE(*,1315) FRUNIT,PWUNIT,MXUNIT,MIUNIT,SRUNIT,HYUNIT,WSUNIT
1750	DO 208 im_wg=1,12
	IF (PWUNIT.EQ.'IN') THEN
	SR(im_wg)=SR(im_wg)*0.03937
	ENDIF
	IF (MXUNIT.EQ.'C') THEN
	QTMX(im_wg) = (SSTX(im_wg) - 32.) * 5./9.
	ELSE
	QTMX(im_wg) = SSTX(im_wg)
	ENDIF
	IF (MIUNIT.EQ.'C') THEN
	QTMN(im_wg) = (SSTN(im_wg) - 32.) * 5./9.
	ELSE
	QTMN (im_wg)= SSTN(im_wg)
	ENDIF
	IF (SRUNIT.EQ.'MJ/M*M/DAY ') THEN
	QRAD(im_wg) = SSRAD(im_wg) * 0.04189
	ELSE 
	QRAD(im_wg) =SSRAD(im_wg)
	ENDIF
	IF (HYUNIT.EQ.'C(DEWP)')THEN
	QVP(im_wg)=237.3*ALOG(SSVP(im_wg)/6.108)/17.27/ &
     &        (1.-ALOG(SSVP(im_wg)/6.108)/17.27)
	ELSEIF (HYUNIT.EQ.'F(DEWP)') THEN
	QVP(im_wg)=237.3*ALOG(SSVP(im_wg)/6.108)/17.27/ &
     &        (1.-ALOG(SSVP(im_wg)/6.108)/17.27)
	QVP(im_wg)=QVP(im_wg)*9./5.+32.
	ELSEIF (HYUNIT.EQ.'MPa(VP)') THEN
	QVP(im_wg) =SSVP(im_wg)	
	ELSE
	QVP(im_wg) = 0.
	ENDIF
	IF(NFS.EQ.1) THEN
	ASSIGN 1801 TO NFST
	ELSEIF(NFS.EQ.2) THEN
	ASSIGN 1802 TO NFST
	ELSEIF(NFS.EQ.3) THEN
	ASSIGN 1803 TO NFST
	ELSE
	ASSIGN 1804 TO NFST
	ENDIF
	GOTO NFST (1801,1802,1803,1804)	
1801	continue
!1801   WRITE(4,1805)DATE(im_wg),XNW(im_wg),SR(im_wg),QTMX(im_wg),&
!       &QTMN(im_wg),QRAD(im_wg)
!	WRITE(*,1805)DATE(im_wg),XNW(im_wg),SR(im_wg),QTMX(im_wg),&
!       &QTMN(im_wg),QRAD(im_wg)
1805  FORMAT(A5,4X,F3.0,4X,F7.2,4X,F5.1,5X,F5.1,6X,F6.1)
	GOTO 208
1802    continue
!1802    WRITE(4,1815)DATE(im_wg),XNW(im_wg),SR(im_wg),QTMX(im_wg),&
!       &QTMN(im_wg),QRAD(im_wg) &
!     &    ,SSWS(im_wg)   
!	WRITE(*,1815)DATE(im_wg),XNW(im_wg),SR(im_wg),QTMX(im_wg),&
!       &QTMN(im_wg),QRAD(im_wg) &
!     &    ,SSWS(im_wg)	
	GOTO 208
1803	continue
!1803	WRITE(4,1815)DATE(im_wg),XNW(im_wg),SR(im_wg),QTMX(im_wg),&
!       &QTMN(im_wg),QRAD(im_wg) &
!     &    ,QVP(im_wg)   
!	WRITE(*,1815)DATE(im_wg),XNW(im_wg),SR(im_wg),QTMX(im_wg),&
!       &QTMN(im_wg),QRAD(im_wg) &
!     &    ,QVP(im_wg)
1815  FORMAT(A5,4X,F3.0,4X,F7.2,3X,F5.1,5X,F5.1,7X,F6.1,6X,F6.2)
	GOTO 208
1804	continue
!1804	WRITE(4,1820)DATE(im_wg),XNW(im_wg),SR(im_wg),QTMX(im_wg),&
!       &QTMN(im_wg),QRAD(im_wg) &
!     &    ,QVP(im_wg),SSWS(im_wg)   
!	WRITE(*,1820)DATE(im_wg),XNW(im_wg),SR(im_wg),QTMX(im_wg),&
!       &QTMN(im_wg),QRAD(im_wg) &
!     &    ,QVP(im_wg),SSWS(im_wg)
1820  FORMAT(A5,4X,F3.0,4X,F7.2,3X,F5.1,5X,F5.1,5X,F6.1,5X,F6.2,5X,F6.2)
 208	CONTINUE
	IF (PWUNIT.EQ.'IN') THEN
	SYR=SYR*0.03937
	ENDIF
	IF(NFS.EQ.1) THEN
	ASSIGN 1901 TO NFST
	ELSEIF(NFS.EQ.2) THEN
	ASSIGN 1902 TO NFST
	ELSEIF(NFS.EQ.3) THEN
	ASSIGN 1903 TO NFST
	ELSE
	ASSIGN 1904 TO NFST
	ENDIF
	GOTO NFST (1901,1902,1903,1904)	
1901	continue
!1901	WRITE (4,1915) SYNW,SYR
!	WRITE (*,1915) SYNW,SYR
1915	FORMAT(/' TOTAL',2X,F4.0,4X,F7.2,7X,'__',8X,'__',10X,'__')
	GOTO 1950
1902	continue
!1902	WRITE (4,1920) SYNW,SYR
!	WRITE (*,1920) SYNW,SYR
1920	FORMAT(/' TOTAL',2X,F4.0,4X,F7.2,6X,'__',8X,'__',11X,'__',10X, &
     &  '__')
	GOTO 1950
1903	continue
!1903	WRITE (4,1920) SYNW,SYR
!	WRITE (*,1920) SYNW,SYR
	GOTO 1950
1904	continue
!1904	WRITE (4,1925) SYNW,SYR
!	WRITE (*,1925) SYNW,SYR
1925	FORMAT(/' TOTAL',2X,F4.0,4X,F7.2,6X,'__',8X,'__',9X,'__',9X, &
     &  '__',9X,'__')
1950	IF (MXUNIT.EQ.'C') THEN
	QYTMX = (SYTX - 32.) * 5./9.
	ELSE
	QYTMX = SYTX
	ENDIF
	IF (MIUNIT.EQ.'C') THEN
	QYTMN = (SYTN - 32.) * 5./9.
	ELSE
	QYTMN = SYTN
	ENDIF
	IF (SRUNIT.EQ.'MJ/M*M/DAY ') THEN
	QYRAD = SYRAD * 0.04189
	ELSE
	QYRAD = SYRAD
	ENDIF
	IF (HYUNIT.EQ.'C(DEWP)')THEN
	QYVP=237.3*ALOG(SYVP/6.108)/17.27/ &
     &        (1.-ALOG(SYVP/6.108)/17.27)
	ELSEIF (HYUNIT.EQ.'F(DEWP)') THEN
	QYVP=237.3*ALOG(SYVP/6.108)/17.27/ &
     &        (1.-ALOG(SYVP/6.108)/17.27)
	QYVP=QYVP*9./5.+32.
	ELSEIF (HYUNIT.EQ.'MPa(VP)') THEN
	QYVP = SYVP
	ELSE
	QYVP = 0.
	ENDIF
	IF(NFS.EQ.1) THEN
	ASSIGN 2001 TO NFST
	ELSEIF(NFS.EQ.2) THEN
	ASSIGN 2002 TO NFST
	ELSEIF(NFS.EQ.3) THEN
	ASSIGN 2003 TO NFST
	ELSE
	ASSIGN 2004 TO NFST
	ENDIF
	GOTO NFST (2001,2002,2003,2004)
2001	continue
!2001	WRITE(4,2025) QYTMX,QYTMN,QYRAD
!	WRITE(*,2025) QYTMX,QYTMN,QYRAD
2025	FORMAT(' AVER.',4X,'__',9X,'__',4X,F5.1,5X,F5.1,6X,F6.1)
	GOTO 2050
2002	continue
!2002	WRITE(4,2030) QYTMX,QYTMN,QYRAD,SYWS
!	WRITE(*,2030) QYTMX,QYTMN,QYRAD,SYWS
	GOTO 2050
2003	continue
!2003	WRITE(4,2030) QYTMX,QYTMN,QYRAD,QYVP
!	WRITE(*,2030) QYTMX,QYTMN,QYRAD,QYVP
2030	FORMAT(' AVER.',4X,'__',9X,'__',3X,F5.1,5X,F5.1,7X,F6.1,6X,F6.2)
	GOTO 2050
2004	continue
!2004	WRITE(4,2035) QYTMX,QYTMN,QYRAD,QYVP,SYWS
!	WRITE(*,2035) QYTMX,QYTMN,QYRAD,QYVP,SYWS
2035	FORMAT(' AVER.',4X,'__',9X,'__',3X,F5.1,5X,F5.1,5X,F6.1,5X,F6.2, &
     &  5X,F6.2)
2050	continue
!2050	WRITE(*,135) FSIMUL
 135	FORMAT(/' * DAILY DATA IN FILE: ',A30)
!	WRITE(*,137) FSUMM
 137	FORMAT(' * YEARLY SUMMARIES IN FILE: ',A30)
	RETURN
!----------------------------------------------------------------------!
END SUBROUTINE WEASHO_WG
!======================================================================!


