       PARAMETER (NSTEP=20)
       DIMENSION PHI(0:NSTEP),S(0:NSTEP)

50     PRINT *,' Enter omega (.le. 0 to stop)'
       READ *, OMEGA
       IF (OMEGA .LE. 0) STOP

       H=1./NSTEP
 
       DO 10 IX=0,NSTEP
          X=IX*H
          S(IX)=H*H*12*X*X
          PHI(IX)=0.
10     CONTINUE

       DO 20 ITER=1,500
          DO 15 IX=1,NSTEP-1
             PHIP=(PHI(IX-1)+PHI(IX+1)+S(IX))/2
             PHI(IX)=(1.-OMEGA)*PHI(IX)+OMEGA*PHIP 
15        CONTINUE
          IF (MOD(ITER-1,20) .EQ. 0) THEN
             E=0.
             DO 30 IX=1,NSTEP  
                E=E+(((PHI(IX)-PHI(IX-1))/H)**2)/2
                E=E-S(IX)*PHI(IX)/H**2
30           CONTINUE
             E=E*H
             PRINT *,' iteration = ',iter,'  energy = ',E
          END IF
20      CONTINUE
        
        GOTO 50
        END
