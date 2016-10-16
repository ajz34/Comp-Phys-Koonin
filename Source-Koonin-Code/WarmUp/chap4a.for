20     PRINT *, ' Enter x, l (l .lt. 0 to stop)'
       READ *, X,L

       IF (L .LT. 0) THEN
           STOP
       ELSE IF (L .EQ. 0) THEN
           PL=0.
       ELSE IF (L .EQ. 1) THEN
           PL=X
       ELSE 
           PM=1.
           PZ=X
           DO 10 IL=1,L-1
              PP=((2*IL+1)*X*PZ-IL*PM)/(IL+1)
              PM=PZ
              PZ=PP
10         CONTINUE
           PL=PZ
       END IF
       
       PRINT *,X,L,PL
       GOTO 20
       END    
