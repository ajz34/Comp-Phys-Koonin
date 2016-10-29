      PROGRAM TB0104

      IMPLICIT NONE

*----------- VARIABLE DECLARATION --------------------------------------

* IMPORTANT SUBROUTINE DESCRIPTION
*   DNTAP: DOUBLE PRECISION ANALYTICALLY NEWTON ROOT FINDER (PROCESS)
*   DSCP: DOUBLE PRECISION SECANT ROOT FINDER (PROCESS)
*   DSRP: DOUBLE PRECISION SEARCH ROOT FINDER (PROCESS)
*   (PROCESS MEANS THEY ONLY DO ONE ITERATION; WE SHOULD WRAP THEM BY 
*   THE LOOPS)
*   THESE SUBROUTINES ARE EXPECTED TO RETURN ROOT VALUES AND UPDATED 
*   METHOD SPECIFIC PARAMETERS.

* FUNCTION DECLARATION
*   DFUNC: THE FUNCTION WE ARE DISCUSSING
*   DDERV: THE 1ST ANALYTIC DERVIATIVE OF THE FUNCTION 'DFUNC'
      DOUBLE PRECISION  DFUNC, DDERV
      EXTERNAL          DFUNC, DDERV

* CALCULATION VARIABLES
*   ITR: ITERATION NUMBER
*   ITRMAX: MAXIMUN ITERATION NUMBER
*   ITRTHR: ITERATION THREAHOLD
*   DSRPR1: SEARCH METHOD PARAMETER 1 --- CURRENT GUESS
*   DSRPR2: SEARCH METHOD PARAMETER 2 --- CURRENT STEP
*   DNTPR: NEWTON METHOD PARAMETER --- CURRENT GUESS
*   DSCPR1: SECANT METHOD PARAMETER 1 --- CURRENT LOWER GUESS
*   DSCRP2: SECANT METHOD PARAMETER 2 --- CURRENT UPPER GUESS
*   DSOL: EXPECTED SOLUTION FOR THE PROBLEM
      DOUBLE PRECISION  DSRPR1, DSRPR2, DNTPR, DSCPR1, DSCPR2, ITRTHR,
     $                  DSOL
      INTEGER           ITR, ITRMAX

*---------- MAIN PROGRAM -----------------------------------------------

* FILE OPERATION

      OPEN (11, FILE='01-06.txt')
 1000 FORMAT ( 1X, I5, 4X, 3(F12.6) )
 1001 FORMAT ( A46 )
 1002 FORMAT ( 1X, A9, 3(A12) )

      WRITE(11,*) 'On the book, the deviation of 5th iteration for' 
      WRITE(11,*) 'secant method is actually wrong'
      WRITE(11,*)
      WRITE(11,*) 'I do not know if different compliers behave diff-'
      WRITE(11,*) 'erently facing NaN. It is possible letting iteration'
      WRITE(11,*) 'numbers reach a very large number.'
      WRITE(11,*)

* CALCULATION

      WRITE(11,1001) '----------------------------------------------'
      WRITE(11,1002) 'Iteration', '      Search', '      Newton',
     $               '      Secant'
      WRITE(11,1002) '         ', '            ', '   Eq.(1.14)',
     $               '   Eq.(1.15)'
      WRITE(11,1001) '----------------------------------------------'

      ITR = 0
      ITRMAX = 100
      ITRTHR = 1.0D-6
      DSRPR1 = 1.0D0
      DSRPR2 = -0.33D0
      DNTPR = 1.085D0
      DSCPR1 = -0.33D0
      DSCPR2 = 1.0D0
      DSOL = 0.0D0
      WRITE(11,1000) ITR, DSRPR1-DSOL, DNTPR-DSOL, DSCPR2-DSOL
  10  CONTINUE
        ITR = ITR + 1
        CALL DSRP(DFUNC, DSRPR1, DSRPR2)
        CALL DNTAP(DFUNC, DDERV, DNTPR)
        CALL DSCP(DFUNC, DSCPR1, DSCPR2)
        WRITE(11,1000) ITR, DSRPR1-DSOL, DNTPR-DSOL, DSCPR2-DSOL
      IF (ITR.LT.ITRMAX.AND.MAX(DABS(DSRPR2), 
     $    DABS(DNTPR-DSQRT(5.0D0)), DABS(DSCPR2-DSCPR1))
     $    .GT.ITRTHR) GOTO 10
      
      WRITE(11,1001) '----------------------------------------------'

      END PROGRAM TB0104

*-----------------------------------------------------------------------

      SUBROUTINE DSRP(DFUNC, DSRPR1, DSRPR2)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      DOUBLE PRECISION  DFUNC, DSRPR1, DSRPR2
      EXTERNAL          DFUNC

      DSRPR1 = DSRPR1 + DSRPR2
      IF (DFUNC(DSRPR1) * DFUNC(DSRPR1-DSRPR2).LE.0.0D0) THEN
        DSRPR2 = - DSRPR2 / 2.0D0
      END IF 

      END SUBROUTINE DSRP

*-----------------------------------------------------------------------

      SUBROUTINE DNTAP(DFUNC, DDERV, DNTPR)

      IMPLICIT NONE

* FUNCTION VARIABLE DECLARATION
      DOUBLE PRECISION  DFUNC, DDERV, DNTPR
      EXTERNAL          DFUNC, DDERV

      DNTPR = DNTPR - DFUNC(DNTPR) / DDERV(DNTPR)

      END SUBROUTINE DNTAP

*-----------------------------------------------------------------------

      SUBROUTINE DSCP(DFUNC, DSCPR1, DSCPR2)

      IMPLICIT NONE
* FUNCTION VARIABLE DECLARATION
      DOUBLE PRECISION  DFUNC, DSCPR1, DSCPR2
      EXTERNAL          DFUNC
* TEMPORARY VARIABLE STORAGE
      DOUBLE PRECISION  DTMP

      DTMP = DSCPR2 - DFUNC(DSCPR2) * (DSCPR2 - DSCPR1) / (DFUNC(DSCPR2)
     $       - DFUNC(DSCPR1))
      DSCPR1 = DSCPR2
      DSCPR2 = DTMP

      END SUBROUTINE

*-----------------------------------------------------------------------

      FUNCTION DFUNC(DVAR)

      IMPLICIT NONE

      DOUBLE PRECISION  DVAR, DFUNC

      DFUNC = DTANH(DVAR)

      END FUNCTION DFUNC

*-----------------------------------------------------------------------

      FUNCTION DDERV(DVAR)

      IMPLICIT NONE

      DOUBLE PRECISION  DVAR, DDERV

      DDERV = 2.0D0 / (1.0D0 + DCOSH(2.0D0 * DVAR))

      END FUNCTION DDERV











