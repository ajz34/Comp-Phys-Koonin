       INTEGER SEED
       DATA SEED/987654321/
       DATA EXACT/.78540/
       FUNC(X)=1./(1.+X**2)

10     PRINT *, ' Enter number of points (0 to stop)'
       READ *, N
       IF (N .EQ. 0) STOP

       SUMF=0.
       SUMF2=0.
       DO 20 IX=1,N
          FX=FUNC(RAN(SEED))
          SUMF=SUMF+FX
          SUMF2=SUMF2+FX**2
20     CONTINUE
       FAVE=SUMF/N
       F2AVE=SUMF2/N
       SIGMA=SQRT((F2AVE-FAVE**2)/N)
       PRINT *,' integral =',FAVE,' +- ',SIGMA,'  error = ',EXACT-FAVE
       
       GOTO 10
       END
