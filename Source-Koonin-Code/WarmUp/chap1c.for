       FUNC(X)=X*X-5.
       TOLX=1.E-06
       X=1.
       FOLD=FUNC(X)
       DX=.5
       ITER=0

10     CONTINUE
         ITER=ITER+1
         X=X+DX
         PRINT *,ITER,X,SQRT(5.)-X
         IF ((FOLD*FUNC(X)) .LT. 0) THEN
            X=X-DX
            DX=DX/2
         END IF
       IF (ABS(DX) .GT. TOLX) GOTO 10

       STOP
       END 
