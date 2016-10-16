       PARAMETER (NSTEP=25)
       DIMENSION PHI(0:NSTEP)
       DIMENSION ALPHA(0:NSTEP),BETA(0:NSTEP),GAMMA(0:NSTEP)

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

       AP=-DTH
       AZ=1.+2*DTH
       ALPHA(NSTEP-1)=0.
       GAMMA(NSTEP-1)=-1./AZ
       DO 15 IX=NSTEP-1,1,-1
          ALPHA(IX-1)=GAMMA(IX)*AP
          GAMMA(IX-1)=-1./(AZ+AP*ALPHA(IX-1))
15     CONTINUE

       DO 20 ITER=1,NITER

          BETA(NSTEP-1)=PHI(NSTEP)
          DO 25 IX=NSTEP-1,1,-1
             BETA(IX-1)=GAMMA(IX)*(AP*BETA(IX)-PHI(IX))
25        CONTINUE

          PHI(0)=0.
          DO 30 IX=0,NSTEP-1
             PHI(IX+1)=ALPHA(IX)*PHI(IX)+BETA(IX)
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
