       X=1.
       EXACT=COS(X)
10     PRINT *, 'ENTER VALUE OF H (.LE. 0 TO STOP)'
       READ *,  H
       IF (H .LE. 0) STOP
       FPRIME=(SIN(X+H)-SIN(X-H))/(2*H)
       DIFF=EXACT-FPRIME
       PRINT 20,H,DIFF
20     FORMAT (' H=',E15.8,5X,'ERROR=',E15.8)
       GOTO 10
       END
