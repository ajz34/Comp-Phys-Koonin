      FUNC(X,Y)=-X*Y
30    PRINT *, ' Enter step size ( .le. 0 to stop)'
      READ *, H
      IF (H .LE. 0.) STOP

      NSTEP=3./H
      Y=1.
      DO 10 IX=0,NSTEP-1
         X=IX*H
         Y=Y+H*FUNC(X,Y)
         DIFF=EXP(-0.5*(X+H)**2)-Y
         PRINT *, IX,X+H,Y,DIFF,DIFF/Y
10    CONTINUE

      DO 20 IX=NSTEP,1,-1
         X=IX*H
         Y=Y-H*FUNC(X,Y)
         DIFF=EXP(-0.5*(X-H)**2)-Y
         PRINT *, IX,X-H,Y,DIFF,DIFF/Y
20    CONTINUE
      PRINT *, Y-1.

      GOTO 30
      END
