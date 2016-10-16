      REAL J(0:50)

100   PRINT *, ' Enter maximum value of n (.le. 50; .lt. 0 to stop)'
      READ *, NMAX
      IF (NMAX .LT. 0) STOP
      IF (NMAX .GT. 50) NMAX=50
      PRINT *,' Enter value of x'
      READ *,X  

      J(NMAX)=0.
      J(NMAX-1)=1.E-20
      DO 10 N=NMAX-1,1,-1
         J(N-1)=(2*N/X)*J(N)-J(N+1)
10    CONTINUE

      SUM=J(O)
      DO 20 N=2,NMAX,2
          SUM=SUM+2*J(N)
20    CONTINUE

      DO 30 N=0,NMAX
         J(N)=J(N)/SUM
         PRINT *,N,J(N)
30    CONTINUE
     
      GOTO 100
      END
