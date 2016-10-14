      PROGRAM TB0101

      IMPLICIT NONE

* DOUBLE PRECISION TABLE
      DOUBLE PRECISION  DH(15), D4(15), D5(15)
* SINGLE PRECISION TABLE
      REAL              SH(15), S4(15), S5(15)

* DUMMIES
      INTEGER I

* DATA DECLARATION
      DATA DH / 0.5D-0, 0.2D-0, 0.1D-0, 0.5D-1, 0.2D-1, 0.1D-1,
     $          0.5D-2, 0.2D-2, 0.1D-2, 0.5D-3, 0.2D-3, 0.1D-3,
     $          0.5D-4, 0.2D-4, 0.1D-4 /
      DATA SH / 0.5E-0, 0.2E-0, 0.1E-0, 0.5E-1, 0.2E-1, 0.1E-1,
     $          0.5E-2, 0.2E-2, 0.1E-2, 0.5E-3, 0.2E-3, 0.1E-3,
     $          0.5E-4, 0.2E-4, 0.1E-4 /
    
* FILE OF WRITTEN FILE
      OPEN (11, FILE = '01-01.txt')

*----------------------------------------------------------------------

* ERROR CALCULATION
      DO 10 I = 1, 15
*   SINGLE ERROR CALCULATION        
        S4(I) = ( SIN(1.0E0 + SH(I)) - 2.0E0 * SIN(1.0E0)
     $            + SIN(1.0E0 - SH(I)) ) 
     $          / ( SH(I) * SH(I) ) + SIN(1.0E0)
        S5(I) = ( - SIN(1.0E0 - 2.0E0 * SH(I)) 
     $            + 16.0E0 * SIN(1.0E0 - SH(I))
     $            - 30.0E0 * SIN(1.0E0)
     $            + 16.0E0 * SIN(1.0E0 + SH(I))
     $            - SIN(1.0E0 + 2.0E0 * SH(I)) )
     $          / ( 12.0E0 * SH(I) * SH(I) ) + SIN(1.0E0)
*   DOUBLE ERROR CALCULATION  
        D4(I) = ( DSIN(1.0D0 + DH(I)) - 2.0D0 * DSIN(1.0D0)
     $            + DSIN(1.0D0 - DH(I)) ) 
     $          / ( DH(I) * DH(I) ) + DSIN(1.0D0)
        D5(I) = ( - DSIN(1.0D0 - 2.0D0 * DH(I)) 
     $            + 16.0D0 * DSIN(1.0D0 - DH(I))
     $            - 30.0D0 * DSIN(1.0D0)
     $            + 16.0D0 * DSIN(1.0D0 + DH(I))
     $            - DSIN(1.0D0 + 2.0D0 * DH(I)) )
     $          / ( 12.0D0 * DH(I) * DH(I) ) + DSIN(1.0D0)

  10  CONTINUE

* DOUBLE ERROR TABLE

      WRITE (11,*)
      WRITE (11,35) '2nd Derivate Single Precision Error Table'
      WRITE (11,*)

      WRITE (11,15) '----------' , '---------------', '---------------' 
      WRITE (11,15) '     h    ' , '        4-Point', '        5-Point'
      WRITE (11,15) '----------' , '---------------', '---------------' 

      DO 20 I = 1, 15
        WRITE (11,25) SH(I), S4(I), S5(I)
  20  CONTINUE
      
      WRITE (11,15) '----------' , '---------------', '---------------' 

* DOUBLE ERROR TABLE

      WRITE (11,*)
      WRITE (11,35) '2nd Derivate Double Precision Error Table'
      WRITE (11,*)

      WRITE (11,15) '----------' , '---------------', '---------------' 
      WRITE (11,15) '     h    ' , '        4-Point', '        5-Point'
      WRITE (11,15) '----------' , '---------------', '---------------' 

      DO 30 I = 1, 15
        WRITE (11,25) DH(I), D4(I), D5(I)
  30  CONTINUE
      
      WRITE (11,15) '----------' , '---------------', '---------------' 

* FORMAT
  15  FORMAT (10A, 2(15A))
  25  FORMAT (F10.5, 2(3X, F12.8))
  35  FORMAT (50A)

      END PROGRAM TB0101
