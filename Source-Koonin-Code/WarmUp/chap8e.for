C subroutine to generate the normally distributed number GAUSS
       SUBROUTINE DIST(GAUSS)
       INTEGER SEED
       DATA SEED/39249187/
  
       GAUSS=0.
       DO 10 I=1,12
          GAUSS=GAUSS+RAN(SEED)
10     CONTINUE
       GAUSS=GAUSS-6.
       RETURN
       END
