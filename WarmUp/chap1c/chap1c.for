      FUNC(X) = X * X - 5.              ! function whose root is sought
      TOLX = 1.E-6                      ! tolerance for the search
      X = 1.                            ! initial guess
      FOLD = FUNC(X)                    ! initial function
      DX = .5                           ! initial step
      ITER = 0                          ! initialize count
 10   CONTINUE
        ITER = ITER + 1                 ! increment iteration count
        X = X + DX                      ! step X
        PRINT *, ITER, X, SQRT(5.)-X    ! output current values
        IF ((FOLD*FUNC(X)).LT.0) THEN
          X = X - DX                    ! if sign changes, back up and 
          DX = DX / 2                   !   halve the step
        END IF
      IF (ABS(DX).GT.TOLX) GOTO 10
      STOP
      END
