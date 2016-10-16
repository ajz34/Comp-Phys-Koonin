      FUNC(X,Y)=-X*Y
20    PRINT *, ' Enter step size ( .le. 0 to stop)'
      READ *, H
      IF (H .LE. 0.) STOP

      NSTEP=3./H
      Y=1.
      DO 10 IX=0,NSTEP-1
         X=IX*H
         Y=Y+H*FUNC(X,Y)
         DIFF=EXP(-0.5*(X+H)**2)-Y
         PRINT *, IX,X+H,Y,DIFF
10    CONTINUE
      GOTO 20
      END
