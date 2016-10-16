10     PRINT *,' Enter value of step size ( .le. 0 to stop)'
       READ *, H
       IF (H .LE. 0) STOP
 
       YMINUS=1
       YZERO=1.-H+H**2/2
       NSTEP=6./H
       DO 20 IX=2,NSTEP
          X=IX*H
          YPLUS=YMINUS-2*H*YZERO
          YMINUS=YZERO
          YZERO=YPLUS
          EXACT=EXP(-X)
          PRINT *, X,EXACT,EXACT-YZERO
20     CONTINUE
       GOTO 10
       END
