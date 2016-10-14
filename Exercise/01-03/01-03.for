      PROGRAM EX0103

      IMPLICIT NONE

*----------- VARIABLE DECLARATION -------------------------------------

* THE SAME AS EX0102
* USE D5P AS THE EVALUATING METHOD (HOWEVER REMOVE THE COMPARASON WITH
* THE TRUE VALUE)
* DFUN1 REFERS TO THE INTEGRAL FROM 0 - 1/2; WHILE DFUN2 1/2 - 1
* DEV REFERS TO THE DEVIATION BETWEEN THE APPROXIMATE AND THE TRUE VALUE

* FUNCTION DECLARATION
      DOUBLE PRECISION  D5P, DFUN1, DFUN2
      EXTERNAL          D5P, DFUN1, DFUN2
* INTEGRAL GRID CALCULATION
      DOUBLE PRECISION  DMIN, DMAX
      INTEGER           N
* RESULTS
      DOUBLE PRECISION  DSOL1, DSOL2, DSOL, DEV

*----------------------------------------------------------------------

* INTEGRAL FROM 0 TO 1/2 : DSOL1
      
      N = 512
      DMIN = 0.0D0
      DMAX = 0.5D0 ** (1.0D0 / 3.0D0)
      DSOL1 = D5P(N, DFUN1, DMIN, DMAX)

* INTEGRAL FROM 1/2 TO 1 : DSOL2

      N = 512
      DMIN = 0.0D0
      DMAX = 0.5D0 ** (2.0D0 / 3.0D0)
      DSOL2 = D5P(N, DFUN2, DMIN, DMAX)

* FINAL RESULT

      DSOL = DSOL1 + DSOL2
      DEV = DSOL - 4.0D0 * DASIN(1.0D0) / DSQRT(3.0D0)

      PRINT *, 'INTEGRAL FROM 0   TO 0.5 : ', DSOL1
      PRINT *, 'INTEGRAL FROM 0.5 TO 1   : ', DSOL2
      PRINT *, 'TOTAL INTEGRAL VALUE     : ', DSOL
      PRINT *, 'DEVIATION FROM EXACT     : ', DEV
      
      END PROGRAM EX0103

*----------------------------------------------------------------------

      FUNCTION D5P(N, DFUNC, DMIN, DMAX)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      DOUBLE PRECISION  D5P, DFUNC, DMIN, DMAX
      INTEGER           FLR
      EXTERNAL          DFUNC, FLR
* DUMMY 
      INTEGER           I
* ERROR CALCULATION VARIABLES
      DOUBLE PRECISION  DH
* CONVERT REAL TO INTEGER FOR THE DUMMY
      INTEGER           K, N

* DEFINING THE INTEGRAL GRIDS
      K = FLR(REAL(DBLE(N) / 4.0D0))
      DH = (DMAX - DMIN) / (4.0D0 * DBLE(K)) 

* INITIAL AND FINAL POINT CALCULATION
      D5P = 7.0D0 * DFUNC(DMIN) - 7.0D0 * DFUNC(DMAX)

* INTERNAL GRID POINT CALCULATION
      DO 10 I = 1, K
        D5P = D5P + 32.0D0 * DFUNC(DMIN + DBLE(I*4-3) * DH)
        D5P = D5P + 12.0D0 * DFUNC(DMIN + DBLE(I*4-2) * DH)
        D5P = D5P + 32.0D0 * DFUNC(DMIN + DBLE(I*4-1) * DH)
        D5P = D5P + 14.0D0 * DFUNC(DMIN + DBLE(I*4) * DH)
  10  CONTINUE

      D5P = D5P * DH * 2.0D0 / 45.0D0

      RETURN
      END FUNCTION D5P

*----------------------------------------------------------------------

      FUNCTION DFUN1(DVAR)

      IMPLICIT NONE

      DOUBLE PRECISION  DFUN1, DVAR

      DFUN1 = 3.0D0 * (1 - DVAR * DVAR * DVAR) ** (-1.0D0 / 3.0D0)

      RETURN
      END FUNCTION DFUN1

*----------------------------------------------------------------------

      FUNCTION DFUN2(DVAR)

      IMPLICIT NONE

      DOUBLE PRECISION  DFUN2, DVAR

      DFUN2 = 1.5D0 * (1 - DVAR ** 1.5D0) ** (-2.0D0 / 3.0D0)

      RETURN
      END FUNCTION DFUN2

*----------------------------------------------------------------------

      FUNCTION FLR(SVAR)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      REAL              SVAR
      INTEGER           FLR

* MAIN PROGRAM

      FLR = INT(SVAR)

      IF (SVAR .LT. 0) FLR = FLR - 1

      RETURN
      END FUNCTION FLR

*----------------------------------------------------------------------

