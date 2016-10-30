      Program TB0201

      Implicit None

*   Initial, Final Point
      Double Precision  xMin
      Parameter         (xMin = 0.0D0)
      Double Precision  xMax
      Parameter         (xMax = 1.0D1)
      Double Precision  yMin
      Parameter         (yMin = 1.0D0)
*   Data Dimension
      Integer           DNum
      Parameter         (DNum = 9)
*   Data Storage
      Integer           GridN(DNum)
      Double Precision  Gap(DNum), y4El(DNum), y5El(DNum), y4Ty(DNum),
     $                  y5Ty(DNum), y4Ip(DNum), y5Ip(DNum)
*   Dummy
      Integer           I
*   Function
      Double Precision  Func, Ext
      External          Func, Ext

      Open(11, File = '02-03.txt')
 1000 Format(F6.3,6F14.7)
 1010 Format(A6,6A14)
 1020 Format(A6,3A28)

*   Data Declaration
      Data GridN /20,50,100,200,500,1000,2000,5000,10000/
      Do 10 I = 1, DNum
        Gap(I) = (xMax - xMin) / DBLE(GridN(I))
  10  Continue

      Do 110 I = 1, DNum
        Call El(y4El(I), y5El(I), xMin, xMax, yMin, GridN(I), Gap(I))
        Call Ty(y4Ty(I), y5Ty(I), xMin, xMax, yMin, GridN(I), Gap(I))
        Call Ip(y4Ip(I), y5Ip(I), xMin, xMax, yMin, GridN(I), Gap(I))
  110 Continue

      Write(11,1020) '------', '----------------------------',
     $               '----------------------------',
     $               '----------------------------'
      Write(11,1020) ' ', 'Euler Method     ', 
     $               'Taylor Series    ', 'Implicit Method   '
      Write(11,1020) ' ', 'Eq.(2.6)       ', 
     $               'Eq.(2.10)      ', 'Eq.(2.18)      '
      Write(11,1010) 'h  ', 'y(4)  ', 'y(5)  ',  'y(4)  ', 
     $               'y(5)  ', 'y(4)  ', 'y(5)  ' 
      Write(11,1020) '------', '----------------------------',
     $               '----------------------------',
     $               '----------------------------'
      Do 20 I = 1, DNum
        Write(11,1000) Gap(I), y4El(I), y5El(I), y4Ty(I), y5Ty(I),
     $                 y4Ip(I), y5Ip(I)
  20  Continue
      Write(11,1020) '------', '----------------------------',
     $               '----------------------------',
     $               '----------------------------'

      End Program TB0201

*-----------------------------------------------------------------------
*
      Subroutine El(y4, y5, xMin, xMax, yMin, GridN, Gap)

      Implicit None

      Double Precision  y4, y5, xMin, xMax, yMin, Gap
      Integer           GridN
*   Dummy
      Integer           I
*   Value for n+1, n+0
      Double Precision  xn1, yn1, xn0, yn0
*   Function
      Double Precision  Func, Ext
      External          Func, Ext
      
      xn0 = xMin
      yn0 = yMin
      Do 10 I = 1, GridN
        xn1 = DBLE(I) / DBLE(GridN) * xMax +
     $        (1.0D0 - DBLE(I) / DBLE(GridN)) * xMin
        yn1 = yn0 + (xn1 - xn0) * Func(xn0,yn0)
        If(DAbs(xn1-4.0D0).LT.1.0D-10) y4 = yn1 - Ext(xn1)
        If(DAbs(xn1-5.0D0).LT.1.0D-10) y5 = yn1 - Ext(xn1)
        xn0 = xn1
        yn0 = yn1
  10  Continue

      End Subroutine El

*-----------------------------------------------------------------------
*
      Subroutine Ty(y4, y5, xMin, xMax, yMin, GridN, Gap)

      Implicit None

      Double Precision  y4, y5, xMin, xMax, yMin, Gap
      Integer           GridN
*   Dummy
      Integer           I
*   Value for n+1, n+0
      Double Precision  xn1, yn1, xn0, yn0
*   Function
      Double Precision  Func, Ext, FuncDx, FuncDy
      External          Func, Ext, FuncDx, FuncDy
      
      xn0 = xMin
      yn0 = yMin
      Do 10 I = 1, GridN
        xn1 = DBLE(I) / DBLE(GridN) * xMax +
     $        (1.0D0 - DBLE(I) / DBLE(GridN)) * xMin
        yn1 = yn0 + (xn1 - xn0) * Func(xn0,yn0) + 0.5D0 * (xn1 - xn0) * 
     $        (xn1 - xn0) * (FuncDx(xn0,yn0) + Func(xn0,yn0) * 
     $        FuncDy(xn0,yn0))
        If(DAbs(xn1-4.0D0).LT.1.0D-10) y4 = yn1 - Ext(xn1)
        If(DAbs(xn1-5.0D0).LT.1.0D-10) y5 = yn1 - Ext(xn1)
        xn0 = xn1
        yn0 = yn1
  10  Continue

      End Subroutine Ty

*-----------------------------------------------------------------------
*
      Subroutine Ip(y4, y5, xMin, xMax, yMin, GridN, Gap)

      Implicit None

      Double Precision  y4, y5, xMin, xMax, yMin, Gap
      Integer           GridN
*   Dummy
      Integer           I
*   Value for n+1, n+0
      Double Precision  xn1, yn1, xn0, yn0
*   Function
      Double Precision  Func, Ext, FuncG
      External          Func, Ext, FuncG
      
      xn0 = xMin
      yn0 = yMin
      Do 10 I = 1, GridN
        xn1 = DBLE(I) / DBLE(GridN) * xMax +
     $        (1.0D0 - DBLE(I) / DBLE(GridN)) * xMin
        yn1 = yn0 * (1.0D0 + 0.5D0 * (xn1 - xn0) * FuncG(xn0)) /
     $        (1.0D0 - 0.5D0 * (xn1 - xn0) * FuncG(xn1))
        If(DAbs(xn1-4.0D0).LT.1.0D-10) y4 = yn1 - Ext(xn1)
        If(DAbs(xn1-5.0D0).LT.1.0D-10) y5 = yn1 - Ext(xn1)
        xn0 = xn1
        yn0 = yn1
  10  Continue

      End Subroutine Ip

*-----------------------------------------------------------------------
*     
      Function Func(x,y)

      Implicit None
      Double Precision  Func, x, y

      Func = - x * y

      End Function Func

*-----------------------------------------------------------------------
*     
      Function FuncDx(x,y)

      Implicit None
      Double Precision  FuncDx, x, y

      FuncDx = - y

      End Function FuncDx

*-----------------------------------------------------------------------
*     
      Function FuncDy(x,y)

      Implicit None
      Double Precision  FuncDy, x, y

      FuncDy = - x

      End Function FuncDy

*-----------------------------------------------------------------------
*     
      Function FuncG(x)

      Implicit None
      Double Precision  FuncG, x

      FuncG = - x

      End Function FuncG

*-----------------------------------------------------------------------
*
      Function Ext(x)

      Implicit None
      Double Precision  Ext, x

      Ext = DExp(-x*x/2.0D0)

      End Function Ext
