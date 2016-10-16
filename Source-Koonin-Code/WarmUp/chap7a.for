       PARAMETER (NSTEP=25)
       DIMENSION PHI(0:NSTEP)

       GAUSS(X,T)=EXP(-20.*(X-.5)**2/(1.+80*T))/SQRT(1+80*T)
       EXACT(X,T)=GAUSS(X,T)-GAUSS(X-1.,T)-GAUSS(X+1.,T)

       H=1./NSTEP

50     PRINT *, ' Enter time step and total time (0 to stop)'
       READ *,DT,TIME
       IF (DT .EQ. 0.) STOP
       NITER=TIME/DT
       DTH=DT/H**2

       T=0.
       PHI(0)=0.
       PHI(NSTEP)=0.
       DO 10 IX=1,NSTEP-1
          PHI(IX)=EXACT(IX*H,T)
10     CONTINUE

       DO 20 ITER=1,NITER
          POLD=0.
          DO 30 IX=1,NSTEP-1
             PNEW=PHI(IX)+DTH*(POLD+PHI(IX+1)-2*PHI(IX))
             POLD=PHI(IX)
             PHI(IX)=PNEW
30        CONTINUE
          IF (MOD(ITER,10) .EQ. 0) THEN
             PRINT *, ' iteration = ', ITER, ' time = ',ITER*DT
             T=ITER*DT
             DO 40 IX=1,NSTEP-1
                DIFF=PHI(IX)-EXACT(IX*H,T)
                PRINT *, ' phi = ', PHI(IX),' error = ', DIFF
40           CONTINUE
          END IF
20     CONTINUE

       GOTO 50
       END
