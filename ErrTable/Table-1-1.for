      PROGRAM TB0101

      IMPLICIT NONE

* DOUBLE PRECISION TABLE
      DOUBLE PRECISION  DH(15), D3(15), D2F(15), D2B(15), D5(15)
* SINGLE PRECISION TABLE
      REAL              SH(15), S3(15), S2F(15), S2B(15), S5(15)

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
      OPEN (11, FILE = 'Table-1-1.txt')

*----------------------------------------------------------------------

* ERROR CALCULATION
      DO 10 I = 1, 15
*   SINGLE ERROR CALCULATION        
        S3(I) = ( SIN(1.0E0 + SH(I)) - SIN(1.0E0 - SH(I)) ) 
     $          / ( 2.0E0 * SH(I) ) - COS(1.0E0)
        S2F(I) = ( SIN(1.0E0 + SH(I)) - SIN(1.0E0) ) 
     $           / SH(I) - COS(1.0E0)
        S2B(I) = ( SIN(1.0E0) - SIN(1.0E0 - SH(I)) )
     $           / SH(I) - COS(1.0E0)
        S5(I) = ( SIN(1.0E0 - 2.0E0 * SH(I)) 
     $            - 8.0E0 * SIN(1.0E0 - SH(I))
     $            + 8.0E0 * SIN(1.0E0 + SH(I))
     $            - SIN(1.0E0 + 2.0E0 * SH(I)) )
     $          / ( 12.0E0 * SH(I) ) - COS(1.0E0)
*   DOUBLE ERROR CALCULATION
        D3(I) = ( DSIN(1.0D0 + DH(I)) - DSIN(1.0D0 - DH(I)) ) 
     $          / ( 2.0D0 * DH(I) ) - DCOS(1.0D0)
        D2F(I) = ( DSIN(1.0D0 + DH(I)) - DSIN(1.0D0) ) 
     $           / DH(I) - DCOS(1.0D0)
        D2B(I) = ( DSIN(1.0D0) - DSIN(1.0D0 - DH(I)) )
     $           / DH(I) - DCOS(1.0D0)
        D5(I) = ( DSIN(1.0D0 - 2.0D0 * DH(I)) 
     $            - 8.0D0 * DSIN(1.0D0 - DH(I))
     $            + 8.0D0 * DSIN(1.0D0 + DH(I))
     $            - DSIN(1.0D0 + 2.0D0 * DH(I)) )
     $          / ( 12.0D0 * DH(I) ) - DCOS(1.0D0)

  10  CONTINUE

* DOUBLE ERROR TABLE

      WRITE (11,*)
      WRITE (11,35) 'Single Precision Error Table'
      WRITE (11,*)

      WRITE (11,15) '-------'    , '-----------', '-----------', 
     $              '-----------', '-----------'
      WRITE (11,15) '       '    , '  Symmetric', '    Forward', 
     $              '   Backward', '  Symmetric'
      WRITE (11,15) '   h   '    , '    3-Point', '    2-Point',
     $              '    2-Point', '    5-Point'
      WRITE (11,15) '       '    , '  Eq.(1.3b)', '  Eq.(1.4a)', 
     $              '  Eq.(1.4b)', '   Eq.(1.5)'
      WRITE (11,15) '-------'    , '-----------', '-----------', 
     $              '-----------', '-----------'

      DO 20 I = 1, 15
        WRITE (11,25) SH(I), S3(I), S2F(I), S2B(I), S5(I)
  20  CONTINUE
      
      WRITE (11,15) '-------'    , '-----------', '-----------', 
     $              '-----------', '-----------'

* DOUBLE ERROR TABLE

      WRITE (11,*)
      WRITE (11,35) 'Double Precision Error Table'
      WRITE (11,*)

      WRITE (11,15) '-------'    , '-----------', '-----------', 
     $              '-----------', '-----------'
      WRITE (11,15) '       '    , '  Symmetric', '    Forward', 
     $              '   Backward', '  Symmetric'
      WRITE (11,15) '   h   '    , '    3-Point', '    2-Point',
     $              '    2-Point', '    5-Point'
      WRITE (11,15) '       '    , '  Eq.(1.3b)', '  Eq.(1.4a)', 
     $              '  Eq.(1.4b)', '   Eq.(1.5)'
      WRITE (11,15) '-------'    , '-----------', '-----------', 
     $              '-----------', '-----------'

      DO 30 I = 1, 15
        WRITE (11,25) DH(I), D3(I), D2F(I), D2B(I), D5(I)
  30  CONTINUE
      
      WRITE (11,15) '-------'    , '-----------', '-----------', 
     $              '-----------', '-----------'

* FORMAT
  15  FORMAT (7A, 4(11A))
  25  FORMAT (F7.5, 4(2X, F9.6))
  35  FORMAT (10X, 30A ,11X)

      END PROGRAM TB0101
