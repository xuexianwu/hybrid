!======================================================================!
!*****THE FOLLOWING SUBROUTINE GENERATES A UNIFORM RANDOM 
!*****NUMBER ON THE INTERVAL 0 - 1
      SUBROUTINE RANDN_WG(YFL,K)
      DIMENSION K(4)
!      DATA K/2510,7692,2456,3765/
      K(4) = 3*K(4)+K(2)
      K(3) = 3*K(3)+K(1)
      K(2)=3*K(2)
      K(1)=3*K(1)
      I=K(1)/1000
      K(1)=K(1)-I*1000
      K(2)=K(2) + I
      I = K(2)/100
      K(2)=K(2)-100*I
      K(3) = K(3)+I
      I = K(3)/1000
      K(3)=K(3)-I*1000
      K(4)=K(4)+I
      I = K(4)/100
      K(4)=K(4)-100*I
      YFL=(((FLOAT(K(1))*.001+FLOAT(K(2)))*.01+FLOAT(K(3)))*.001 &
     &    +FLOAT(K(4)))*.01
      RETURN
      END SUBROUTINE RANDN_WG
!======================================================================!
