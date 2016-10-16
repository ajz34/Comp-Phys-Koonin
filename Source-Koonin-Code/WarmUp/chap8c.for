       INTEGER SEED
       DATA SEED/3274927/
       EXACT=4.*ATAN(1.)

20     PRINT *, ' Enter number of points (0 to stop)'
       READ *, N
       IF (N .EQ. 0) STOP

       ICOUNT=0
       DO 10 IX=1,N
          X=RAN(SEED)
          Y=RAN(SEED)
          IF ((X**2+Y**2) .LE. 1.) ICOUNT=ICOUNT+1
10     CONTINUE
       PI4=REAL(ICOUNT)/N
       SIGMA=SQRT(PI4*(1.-PI4)/N)
       PRINT *,' pi = ', 4*PI4,' +- ',4*SIGMA,'   error = ',EXACT-4*PI4

       GOTO 20
       END
