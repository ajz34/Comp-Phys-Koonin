C subroutine to take a step in the Metropolis algorithm
       SUBROUTINE METROP(X1,X2,WEIGHT,DELTA)
       EXTERNAL WEIGHT
       INTEGER SEED
       DATA SEED/39249187/

       X1T=X1+DELTA*(2*RAN(SEED)-1)
       X2T=X2+DELTA*(2*RAN(SEED)-1)
       RATIO=WEIGHT(X1T,X2T)/WEIGHT(X1,X2)

       IF (RATIO .GT. RAN(SEED)) THEN
           X1=X1T
           X2=X2T
       END IF

       RETURN
       END
