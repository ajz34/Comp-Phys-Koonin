       REAL K
       K=1.
       DK=1.
       TOLK=1.E-05
       CALL INTGRT(K,PHIP)
       PHIOLD=PHIP

10     CONTINUE
         K=K+DK
         CALL INTGRT(K,PHIP)
         IF (PHIP*PHIOLD .LT. 0) THEN
             K=K-DK
             DK=DK/2
         END IF
       IF (ABS(DK) .GT. TOLK) GOTO 10

       EXACT=4.*ATAN(1.)
       PRINT *, ' eigenvalue, error =',K,EXACT-K
       STOP
       END
                    

       SUBROUTINE INTGRT(K,PHIP)
       REAL K
       DATA NSTEP/100/     

       H=1./NSTEP
       PHIM=0.
       PHIZ=.01
       CONST=(K*H)**2/12.

       DO 10 IX=1,NSTEP-1
          PHIP=2*(1.-5.*CONST)*PHIZ -(1.+CONST)*PHIM
          PHIP=PHIP/(1+CONST)
          PHIM=PHIZ
          PHIZ=PHIP
10     CONTINUE

       PRINT *, K,PHIP
       RETURN
       END
