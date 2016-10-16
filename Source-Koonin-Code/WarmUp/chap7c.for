       PARAMETER (NSTEP=20)
       DIMENSION PHI(0:NSTEP)

       H=1./NSTEP

100    PRINT *, ' Enter step size ( .eq. 0 to stop)'
       READ *, DT
       IF (DT .EQ. 0.) STOP
       DTH=DT/H**2
   
       DO 10 IX=0,NSTEP
          X=IX*H
          PHI(IX)=X*(1.-X)
 10    CONTINUE
       CALL NORM(PHI,NSTEP)

       DO 20 ITER=1,100
          POLD=0.
          DO 30 IX=1,NSTEP-1
             PNEW=PHI(IX)+DTH*(POLD+PHI(IX+1)-2*PHI(IX))
             POLD=PHI(IX)
             PHI(IX)=PNEW
30        CONTINUE
          CALL NORM(PHI,NSTEP)

          IF (MOD(ITER,4) .EQ. 0) THEN
            E=0.
            DO 50 IX=1,NSTEP
               E=E+(PHI(IX)-PHI(IX-1))**2
50          CONTINUE
            E=E/(H)
            PRINT *,' iteration = ', ITER, '  energy = ', E
          END IF
20     CONTINUE
       GOTO 100
       END

       SUBROUTINE NORM(PHI,NSTEP)
       DIMENSION PHI(0:NSTEP)
      
       XNORM=0.
       DO 10 IX=0,NSTEP
          XNORM=XNORM+PHI(IX)**2
10     CONTINUE
       XNORM=SQRT(NSTEP/XNORM)
       DO 20 IX=0,NSTEP
          PHI(IX)=PHI(IX)*XNORM
20     CONTINUE
       RETURN
       END
