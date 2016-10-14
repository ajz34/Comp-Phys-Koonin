      PROGRAM TB0103
      
      IMPLICIT NONE

*----------- VARIABLE DECLARATION -------------------------------------

* FUNCTION DECLARATION
*   SxP: APPROXIMATE FUNCTIONS BY x POINTS
*   SFUNC: THE FUNCTION WE ARE DISCUSSING
*   SINTG: INTEGRAL OF SFUNC
*   D...: ALL WORKING IN THE DOUBLE PRECISION
*   S...: ALL WORKING IN THE SINGLE PRECISION
      REAL              S2P, S3P, S4P, S5P, SFUNC, SINTG
      EXTERNAL          S2P, S3P, S4P, S5P, SFUNC, SINTG
      DOUBLE PRECISION  D2P, D3P, D4P, D5P, DFUNC, DINTG
      EXTERNAL          D2P, D3P, D4P, D5P, DFUNC, DINTG

* ERROR CALCULATION VARIABLES
*   N: LOADED INTERNAL POINTS
*   SH: LOADED POINT INTERVALS
*   SMIN: INTEGRAL LOWER LIMIT
*   SMAX: INTEGRAL UPPER LIMIT
      REAL              SH, SMIN, SMAX
      DOUBLE PRECISION  DH, DMIN, DMAX
      INTEGER           N

* DUMMY
      INTEGER           I

*----------------------------------------------------------------------

* FILE OPERATION

      OPEN (11, FILE='01-02.txt')
      
 1000 FORMAT ( 1X, I4, F12.7, 4(F12.6) )  
 1001 FORMAT ( 1X, A4, A12, 4(A12) )
 1002 FORMAT ( A1, A32, A32, A1 )

* SINGLE PRECISION

      WRITE(11,*)
      WRITE(11,*) 'Evaluated function: Integral of sin(x)'
      WRITE(11,*)
      WRITE(11,1002) '-','-----------------------------Sin',
     $               'gle-----------------------------','-'
      WRITE(11,1001) '    ', '            ', ' Trapezoidal', 
     $               '     Simpson', ' Simpson 3/8', '       Boole'
      WRITE(11,1001) '   N', '           h', '    Eq.(1.9)', 
     $               '   Eq.(1.12)', '  Eq.(1.13a)', '  Eq.(1.13b)'
      WRITE(11,1002) '-','--------------------------------',
     $               '--------------------------------','-'

      N = 2
      SMAX = 1.0E0
      SMIN = 0.0E0
      DO 10 I = 1, 6
        N = N * 2
        SH = (SMAX - SMIN) / REAL(N)
        WRITE(11,1000) N, SH,
     $                 S2P(N, SFUNC, SINTG, SMIN, SMAX), 
     $                 S3P(N, SFUNC, SINTG, SMIN, SMAX),
     $                 S4P(N, SFUNC, SINTG, SMIN, SMAX),
     $                 S5P(N, SFUNC, SINTG, SMIN, SMAX)
  10  CONTINUE

      WRITE(11,1002) '-','--------------------------------',
     $               '--------------------------------','-'

* DOUBLE PRECISION

      WRITE(11,*)
      WRITE(11,1002) '-','-----------------------------Dou',
     $               'ble-----------------------------','-'
      WRITE(11,1001) '    ', '            ', ' Trapezoidal', 
     $               '     Simpson', ' Simpson 3/8', '       Boole'
      WRITE(11,1001) '   N', '           h', '    Eq.(1.9)', 
     $               '   Eq.(1.12)', '  Eq.(1.13a)', '  Eq.(1.13b)'
      WRITE(11,1002) '-','--------------------------------',
     $               '--------------------------------','-'

      N = 2
      DMAX = 1.0D0
      DMIN = 0.0D0
      DO 20 I = 1, 6
        N = N * 2
        DH = (DMAX - DMIN) / DBLE(N)
        WRITE(11,1000) N, DH,
     $                 D2P(N, DFUNC, DINTG, DMIN, DMAX), 
     $                 D3P(N, DFUNC, DINTG, DMIN, DMAX),
     $                 D4P(N, DFUNC, DINTG, DMIN, DMAX),
     $                 D5P(N, DFUNC, DINTG, DMIN, DMAX)
  20  CONTINUE

      WRITE(11,1002) '-','--------------------------------',
     $               '--------------------------------','-'

      END

*----------------------------------------------------------------------

      FUNCTION S2P(N, SFUNC, SINTG, SMIN, SMAX)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      REAL              S2P, SFUNC, SINTG, SMIN, SMAX
      EXTERNAL          SFUNC, SINTG
* DUMMY 
      INTEGER           I
* ERROR CALCULATION VARIABLES
      REAL              SH
* CONVERT REAL TO INTEGER FOR THE DUMMY
      INTEGER           K, N

      SH = (SMAX - SMIN) / REAL(N) 

      S2P = SFUNC(SMIN) + SFUNC(SMAX)

      K = N-1

      DO 10 I = 1, K
        S2P = S2P + 2.0E0 * SFUNC(SMIN + REAL(I) * SH)
  10  CONTINUE

      S2P = S2P * SH / 2.0E0 - SINTG(SMAX) + SINTG(SMIN)

      RETURN
      END FUNCTION S2P

*----------------------------------------------------------------------

      FUNCTION S3P(N, SFUNC, SINTG, SMIN, SMAX)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      REAL              S3P, SFUNC, SINTG, SMIN, SMAX
      INTEGER           FLR
      EXTERNAL          SFUNC, SINTG, FLR
* DUMMY 
      INTEGER           I
* ERROR CALCULATION VARIABLES
      REAL              SH
* CONVERT REAL TO INTEGER FOR THE DUMMY
      INTEGER           K, N

* DEFINING THE INTEGRAL GRIDS
      K = FLR(REAL(N) / 2.0E0)
      SH = (SMAX - SMIN) / (2.0E0 * REAL(K)) 

* INITIAL AND FINAL POINT CALCULATION
      S3P = SFUNC(SMIN) - SFUNC(SMAX)

* INTERNAL GRID POINT CALCULATION
      DO 10 I = 1, K
        S3P = S3P + 4.0E0 * SFUNC(SMIN + REAL(I*2-1) * SH)
        S3P = S3P + 2.0E0 * SFUNC(SMIN + REAL(I*2) * SH)
  10  CONTINUE

      S3P = S3P * SH / 3.0E0 - SINTG(SMAX) + SINTG(SMIN)

      RETURN
      END FUNCTION S3P

*----------------------------------------------------------------------

      FUNCTION S4P(N, SFUNC, SINTG, SMIN, SMAX)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      REAL              S4P, SFUNC, SINTG, SMIN, SMAX
      INTEGER           FLR
      EXTERNAL          SFUNC, SINTG, FLR
* DUMMY 
      INTEGER           I
* ERROR CALCULATION VARIABLES
      REAL              SH
* CONVERT REAL TO INTEGER FOR THE DUMMY
      INTEGER           K, N

* DEFINING THE INTEGRAL GRIDS
      K = FLR(REAL(N) / 3.0E0)
      SH = (SMAX - SMIN) / (3.0E0 * REAL(K)) 

* INITIAL AND FINAL POINT CALCULATION
      S4P = SFUNC(SMIN) - SFUNC(SMAX)

* INTERNAL GRID POINT CALCULATION
      DO 10 I = 1, K
        S4P = S4P + 3.0E0 * SFUNC(SMIN + REAL(I*3-2) * SH)
        S4P = S4P + 3.0E0 * SFUNC(SMIN + REAL(I*3-1) * SH)
        S4P = S4P + 2.0E0 * SFUNC(SMIN + REAL(I*3) * SH)
  10  CONTINUE

      S4P = S4P * SH * 3.0E0 / 8.0E0 - SINTG(SMAX) + SINTG(SMIN)

      RETURN
      END FUNCTION S4P

*-----------------------------------------------------------------------

      FUNCTION S5P(N, SFUNC, SINTG, SMIN, SMAX)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      REAL              S5P, SFUNC, SINTG, SMIN, SMAX
      INTEGER           FLR
      EXTERNAL          SFUNC, SINTG, FLR
* DUMMY 
      INTEGER           I
* ERROR CALCULATION VARIABLES
      REAL              SH
* CONVERT REAL TO INTEGER FOR THE DUMMY
      INTEGER           K, N

* DEFINING THE INTEGRAL GRIDS
      K = FLR(REAL(N) / 4.0E0)
      SH = (SMAX - SMIN) / (4.0E0 * REAL(K)) 

* INITIAL AND FINAL POINT CALCULATION
      S5P = 7.0E0 * SFUNC(SMIN) - 7.0E0 * SFUNC(SMAX)

* INTERNAL GRID POINT CALCULATION
      DO 10 I = 1, K
        S5P = S5P + 32.0E0 * SFUNC(SMIN + REAL(I*4-3) * SH)
        S5P = S5P + 12.0E0 * SFUNC(SMIN + REAL(I*4-2) * SH)
        S5P = S5P + 32.0E0 * SFUNC(SMIN + REAL(I*4-1) * SH)
        S5P = S5P + 14.0E0 * SFUNC(SMIN + REAL(I*4) * SH)
  10  CONTINUE

      S5P = S5P * SH * 2.0E0 / 45.0E0 - SINTG(SMAX) + SINTG(SMIN)

      RETURN
      END FUNCTION S5P

*-----------------------------------------------------------------------

      FUNCTION D2P(N, DFUNC, DINTG, DMIN, DMAX)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      DOUBLE PRECISION  D2P, DFUNC, DINTG, DMIN, DMAX
      EXTERNAL          DFUNC, DINTG
* DUMMY 
      INTEGER           I
* ERROR CALCULATION VARIABLES
      DOUBLE PRECISION  DH
* CONVERT REAL TO INTEGER FOR THE DUMMY
      INTEGER           K, N

      DH = (DMAX - DMIN) / DBLE(N) 

      D2P = DFUNC(DMIN) + DFUNC(DMAX)

      K = N-1

      DO 10 I = 1, K
        D2P = D2P + 2.0D0 * DFUNC(DMIN + DBLE(I) * DH)
  10  CONTINUE

      D2P = D2P * DH / 2.0D0 - DINTG(DMAX) + DINTG(DMIN)

      RETURN
      END FUNCTION D2P

*----------------------------------------------------------------------

      FUNCTION D3P(N, DFUNC, DINTG, DMIN, DMAX)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      DOUBLE PRECISION  D3P, DFUNC, DINTG, DMIN, DMAX
      INTEGER           FLR
      EXTERNAL          DFUNC, DINTG, FLR
* DUMMY 
      INTEGER           I
* ERROR CALCULATION VARIABLES
      DOUBLE PRECISION  DH
* CONVERT REAL TO INTEGER FOR THE DUMMY
      INTEGER           K, N

* DEFINING THE INTEGRAL GRIDS
      K = FLR(REAL(DBLE(N) / 2.0D0))
      DH = (DMAX - DMIN) / (2.0D0 * DBLE(K)) 

* INITIAL AND FINAL POINT CALCULATION
      D3P = DFUNC(DMIN) - DFUNC(DMAX)

* INTERNAL GRID POINT CALCULATION
      DO 10 I = 1, K
        D3P = D3P + 4.0D0 * DFUNC(DMIN + DBLE(I*2-1) * DH)
        D3P = D3P + 2.0D0 * DFUNC(DMIN + DBLE(I*2) * DH)
  10  CONTINUE

      D3P = D3P * DH / 3.0D0 - DINTG(DMAX) + DINTG(DMIN)

      RETURN
      END FUNCTION D3P

*----------------------------------------------------------------------

      FUNCTION D4P(N, DFUNC, DINTG, DMIN, DMAX)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      DOUBLE PRECISION  D4P, DFUNC, DINTG, DMIN, DMAX
      INTEGER           FLR
      EXTERNAL          DFUNC, DINTG, FLR
* DUMMY 
      INTEGER           I
* ERROR CALCULATION VARIABLES
      DOUBLE PRECISION  DH
* CONVERT REAL TO INTEGER FOR THE DUMMY
      INTEGER           K, N

* DEFINING THE INTEGRAL GRIDS
      K = FLR(REAL(DBLE(N) / 3.0D0))
      DH = (DMAX - DMIN) / (3.0D0 * DBLE(K)) 

* INITIAL AND FINAL POINT CALCULATION
      D4P = DFUNC(DMIN) - DFUNC(DMAX)

* INTERNAL GRID POINT CALCULATION
      DO 10 I = 1, K
        D4P = D4P + 3.0D0 * DFUNC(DMIN + DBLE(I*3-2) * DH)
        D4P = D4P + 3.0D0 * DFUNC(DMIN + DBLE(I*3-1) * DH)
        D4P = D4P + 2.0D0 * DFUNC(DMIN + DBLE(I*3) * DH)
  10  CONTINUE

      D4P = D4P * DH * 3.0D0 / 8.0D0 - DINTG(DMAX) + DINTG(DMIN)

      RETURN
      END FUNCTION D4P

*-----------------------------------------------------------------------

      FUNCTION D5P(N, DFUNC, DINTG, DMIN, DMAX)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      DOUBLE PRECISION  D5P, DFUNC, DINTG, DMIN, DMAX
      INTEGER           FLR
      EXTERNAL          DFUNC, DINTG, FLR
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

      D5P = D5P * DH * 2.0D0 / 45.0D0 - DINTG(DMAX) + DINTG(DMIN)

      RETURN
      END FUNCTION D5P

*-----------------------------------------------------------------------

      FUNCTION SFUNC(SVAR)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      REAL              SFUNC, SVAR

      SFUNC = SIN(SVAR)

      RETURN
      END FUNCTION SFUNC

*-----------------------------------------------------------------------

      FUNCTION SINTG(SVAR)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      REAL              SINTG, SVAR

      SINTG = -COS(SVAR)

      RETURN
      END FUNCTION SINTG

*----------------------------------------------------------------------

      FUNCTION DFUNC(DVAR)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      DOUBLE PRECISION  DFUNC, DVAR

      DFUNC = DSIN(DVAR)

      RETURN
      END FUNCTION DFUNC

*-----------------------------------------------------------------------

      FUNCTION DINTG(DVAR)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      DOUBLE PRECISION  DINTG, DVAR

      DINTG = -DCOS(DVAR)

      RETURN
      END FUNCTION DINTG

*-----------------------------------------------------------------------

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














































