      FUNC(X) = EXP(X)                  ! function to integrate
      EXACT = EXP(1.) - 1.
  30  PRINT *,'ENTER N EVEN (.LT. 2 TO STOP)'
      READ *, N
      IF (N .LT. 2) STOP
      IF (MOD(N,2) .NE. 0) N=N+1
      H = 1./N
      SUM = FUNC(0.)                    ! contribution from X=0
      FAC = 2.                          ! factor for Simpson's rule
      DO 10 I=1,N-1                     ! loop over lattice points
        IF (FAC .EQ. 2.) THEN           ! factros alternate
          FAC = 4.
        ELSE
          FAC = 2.
        END IF
        X = I * H                       ! X at this point
        SUM = SUM + FAC * FUNC(X)       ! contribution to the integral
  10  CONTINUE
      SUM = SUM + FUNC(1.)              ! contribution from X=1
      XINT = SUM * H / 3.
      DIFF = EXACT - XINT
      PRINT 20, N,DIFF
  20  FORMAT (5X,'N=',I5,5X,'ERROR=',E15.8)
      GOTO 30                           ! get another value of N
      END
