!======================================================================!
!*****THE FOLLOWING SUBROUTINE GENERATES DAILY WEATHER DATA FOR
!*****ONE YEAR.
      SUBROUTINE daily_wg
!----------------------------------------------------------------------!
      USE WGEN_PARAMS_WG
      USE GEN_ATMOS_WG
      USE CONTROL_DEFINITIONS_WG
      use physical_parameters_wg
!----------------------------------------------------------------------!
      DIMENSION A(4,4),B(4,4),E(4),Ri(4),X(4),RR(4)
      COMMON NFS	
      DATA A/0.753,0.188,-0.023,0.094,0.089,0.683,0.136,0.006,0.092, &
     &0.013,0.763,0.025,0.023,0.065,0.030,0.771/
      DATA B/-0.238,0.169,-0.522,0.230,-0.268,-0.291,-0.158,-0.423, &
     &     -0.075,-0.315,0.041,0.249,0.212,-0.087,-0.136,-0.024/
      imo = 1
      DO 50 IDAY=1,IDAYS
      IF(IDAYS .EQ. 366) GO TO 2
      IF(IDAY .GT. NI(imo)) imo = imo + 1
      GO TO 4
 2    IF(IDAY .GT. NII(imo)) imo = imo + 1
 4    CONTINUE
!     CALL RANDN_WG(RN4,K)
      CALL RANDN_WG(RN,K)	
      IF((NFS.EQ.1).OR.(NFS.EQ.3)) GOTO 5	
      WIND = CC(imo)*(-ALOG(RN))**(1.0/(KAY(imo)+1.0E-8))
!     WIND = CC(imo)*(-ALOG(RN4))**(1.0/(KAY(imo)+1.0E-8))
      WINDSP(IDAY) = WIND
!*****WINDSP(IDAY) IS GENERATED WIND SPEED FOR IDAY
5     IF(KGEN .EQ. 2) GO TO 15
!*****DETERMINE WET OR DRY DAY USING MARKOV CHAIN MODEL
!     CALL RANDN_WG(RN,K)
      IF(IP-0) 7,7,10
 7    IF(RN - PWD(imo  ))11,11,8
 8    IP = 0
      RAIN(IDAY) = 0.
      GO TO 18
 10   IF(RN-PWW(imo  ))11,11,8
 11   IP = 1
!*****DETERMINE RAINFALL AMOUNT FOR WET DAYS USING GAMMA 
!*    DISTRIBUTION
      AA = 1./(ALPHAP(imo)+1.0E-8)
      AB = 1./(1.-ALPHAP(imo)+1.0E-8)
      TR1 = EXP(-18.42/(AA+1.0E-8))
      TR2 = EXP(-18.42/(AB+1.0E-8))
      SUM = 0.
      SUM2 = 0.
 12   CALL RANDN_WG(RN1,K)
      CALL RANDN_WG(RN2,K)
      IF(RN1-TR1) 61,61,62
 61   S1 = 0.
      GO TO 63
 62   S1 = RN1**AA
 63   IF(RN2-TR2) 64,64,65
 64   S2 = 0.
      GO TO 66
 65   S2 = RN2**AB
 66   S12 = S1 + S2
!----------------------------------------------------------------------!
!  Had problem because S12 = 0. Tried to fix.
!----------------------------------------------------------------------!
      IF (S12 .EQ. 0.0) THEN
        WRITE (99,*) 'S12=0, corrected'
        S12 = 1.0e-6
      ENDIF
!----------------------------------------------------------------------!
      IF(S12-1.) 13,13,12
 13   Z = S1/(S12+1.0E-8)
      CALL RANDN_WG(RN3,K)
      RAIN(IDAY) = -Z*ALOG(RN3)*BETA(imo)*RCF(imo)
!      IF ((RAIN(IDAY).GE. 0.0) .AND. (RAIN(IDAY).LE. 2000.0)) THEN
!        CONTINUE
!      ELSE
!        WRITE (99,*) 'Problem in WGEN1'
!        WRITE (99,*) 'IDAY=',IDAY,'imo=',imo
!        WRITE (99,*) 'IP=',IP,'RN=',RN
!        WRITE (99,*) 'RAIN(IDAY)=',RAIN(IDAY),'Z=',Z
!        WRITE (99,*) 'RN3=',RN3,'BETA=',BETA(imo),'RCF(imo)=',RCF(imo)
!        WRITE (99,*) 'S1=',S1,'S2=',S2,'S12,=',S12
!        WRITE (99,*) 'RN1=',RN1,'RN2=',RN2,'AA=',AA
!        WRITE (99,*) 'AB=',AB,'TR1=',TR1,'TR2=',TR2
!        WRITE (99,*) 'ALPHAP=',ALPHAP(imo),'RN=',RN,'PWD=',PWD(imo)
!        WRITE (99,*) 'PWW=',PWW(imo)
!        STOP
!      ENDIF
!*****RAIN(IDAY) IS GENERATED RAINFALL FOR IDAY
 15   IF(RAIN(IDAY)) 16,16,17
 16   IP = 0
      GO TO 18
 17   IP = 1
 18   IF(IP-1) 25,26,26
!*****GENERATE  dTMAX, dTMIN,RAD AND VAP FOR IDAY
 25   RM=RMO(IDAY)
      RS = RSO(IDAY)
      TXXM = TXM(IDAY)
      TXXS = TXS(IDAY)
      GO TO 27
 26   RM = RM1(IDAY)
      RS = RS1(IDAY)
      TXXM = TXM1(IDAY)
      TXXS = TXS1(IDAY)
 27   CONTINUE
      VPMM = VPM(IDAY) 
      VPSS = VPS(IDAY)	
      IF ((NFS.EQ.1).OR.(NFS.EQ.2)) THEN
	LL=3
	ELSE
	LL=4
      ENDIF
	DO 30 KO = 1,LL
 131  AA = 0.
      CALL RANDN_WG(RN1,K)
      CALL RANDN_WG(RN2,K)
      V = SQRT(-2.*ALOG(RN1))*COS(6.283185*RN2)
      IF(ABS(V) .GT. 2.5) GO TO 131
      E(KO) = V
 30   CONTINUE
      DO 31 I = 1,LL
      Ri(I) = 0.
      RR(I) = 0.
 31   CONTINUE
      DO 32 I = 1,LL
      DO 32 J = 1,LL
      Ri(I) = Ri(I)+B(I,J)*E(J)
      RR(I) = RR(I) + A(I,J)*XIM1(J)
 32   CONTINUE
      DO 37 KO = 1,LL
      X(KO) = Ri(KO) + RR(KO)
      XIM1(KO) = X(KO)
 37   CONTINUE
      dTMAX(IDAY) = X(1) * TXXS + TXXM
      dTMIN(IDAY) = X(2)*TNS(IDAY)+TNM(IDAY)
      IF( dTMIN(IDAY) .GT. dTMAX(IDAY)) GO TO 38
      GO TO 39
 38   TMM = dTMAX(IDAY)
      dTMAX(IDAY) =  dTMIN(IDAY)
      dTMIN(IDAY) = TMM
 39   CONTINUE
      dTMAX(IDAY)=dTMAX(IDAY)+TCFMAX(imo)
      dTMIN(IDAY)= dTMIN(IDAY)+TCFMIN(imo)
!*****dTMAX(IDAY) IS GENERATED TMAX FOR IDAY
!*****dTMIN(IDAY) IS GENERATED  dTMIN FOR IDAY
      RAD(IDAY) = X(3)*RS+RM
      RMIN = 0.05*RC(IDAY)
      IF(RAD(IDAY) .LT. RMIN) RAD(IDAY) = RMIN
      IF (RAD(IDAY) .GT. RC(IDAY)) RAD(IDAY) = RC(IDAY)
!*****RAD(IDAY) IS GENERATED RAD FOR IDAY
      pressure (iday) = p_mm (imo) ! Set atmospheric pressure (Pa)
      IF ((NFS.EQ.1).OR.(NFS.EQ.2)) GOTO 50
      VAP(IDAY) = X(4)*VPSS+VPMM
!  Added by ADF 17/03/2005.
      IF (VAP (IDAY) .LT. 0.0) VAP (IDAY) = 0.0
      IF ((RAIN(IDAY).GE. 0.0) .AND. (RAIN(IDAY).LE. 5000.0)) THEN
        CONTINUE
      ELSE
        WRITE (99,*) 'Problem in WGEN2'
        WRITE (99,*) 'IDAY=',IDAY,'imo=',imo
        WRITE (99,*) 'IP=',IP,'RN=',RN
        WRITE (99,*) 'RAIN(IDAY)=',RAIN(IDAY),'Z=',Z
        WRITE (99,*) 'RN3=',RN3,'BETA=',BETA(imo),'RCF(imo)=',RCF(imo)
        WRITE (99,*) 'S1=',S1,'S2=',S2,'S12,=',S12
        WRITE (99,*) 'RN1=',RN1,'RN2=',RN2,'AA=',AA
        WRITE (99,*) 'AB=',AB,'TR1=',TR1,'TR2=',TR2
        WRITE (99,*) 'ALPHAP=',ALPHAP(imo),'RN=',RN,'PWD=',PWD(imo)
        WRITE (99,*) 'PWW=',PWW(imo)
        STOP
      ENDIF
!*****VAP(IDAY) IS GENERATED VAP FOR IDAY
     ! WRITE (99,8000) IDAY,RAIN(IDAY),dTMAX(IDAY),dTMIN(IDAY), &
     !&                RAD(IDAY),    &
     !&                VAP(IDAY),WINDSP(IDAY)
 8000 FORMAT (I5,6F10.3)
      !----------------------------------------------------------------!
50    CONTINUE
      RETURN
!----------------------------------------------------------------------!
END SUBROUTINE daily_wg
!======================================================================!
