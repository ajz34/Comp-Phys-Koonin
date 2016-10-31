      Program Ex0204

      Implicit None

*   Initial, Final Point
      Double Precision  xMin
      Parameter         (xMin = 0.0D0)
      Double Precision  xMax
      Parameter         (xMax = 3.0D0)
      Double Precision  yMin
      Parameter         (yMin = 1.0D0)
*   Data Dimension
      Integer           DNum
      Parameter         (DNum = 9)
*   Data Storage
      Integer           GridN(DNum)
      Double Precision  Gap(DNum), y1RK1(DNum), y3RK1(DNum),yrRK1(DNum),
     $                  y1RK2(DNum), y3RK2(DNum), yrRK2(DNum), 
     $                  y1RK3(DNum), y3RK3(DNum), yrRK3(DNum)
*   Dummy
      Integer           I
*   Function
      Double Precision  Func, Ext
      External          Func, Ext

      Open(11, File = '02-04.txt')
 1000 Format(F6.3,9F14.7)
 1010 Format(A6,9A14)
 1020 Format(A6,3A42)

*   Data Declaration
      Data GridN /6,15,30,60,150,300,600,1500,3000/
      Do 10 I = 1, DNum
        Gap(I) = (xMax - xMin) / DBLE(GridN(I))
  10  Continue

      Do 110 I = 1, DNum
        Call RK1(y1RK1(I), y3RK1(I), yrRK1(I), xMin,xMax,yMin, GridN(I),
     $          Gap(I))
        Call RK2(y1RK2(I), y3RK2(I), yrRK2(I), xMin,xMax,yMin, GridN(I),
     $          Gap(I))
        Call RK3(y1RK3(I), y3RK3(I), yrRK3(I), xMin,xMax,yMin, GridN(I),
     $          Gap(I))
  110 Continue

      Write(11,1020) '------', 
     $               '------------------------------------------',
     $               '------------------------------------------',
     $               '------------------------------------------'
      Write(11,1020) ' |', 
     $               'Runge-Kutta 2-order       |', 
     $               'Runge-Kutta 3-order       |', 
     $               'Runge-Kutta 4-order       |'
      Write(11,1020) ' |', 
     $               'Eq.(2.22)            |', 
     $               'Eq.(2.24)            |', 
     $               'Eq.(2.25)            |'
      Write(11,1010) 'h |', 
     $               'y(1)  ', 'y(3)  ', 'Retrun |',
     $               'y(1)  ', 'y(3)  ', 'Return |',
     $               'y(1)  ', 'y(3)  ', 'Retrun |'
      Write(11,1020) '------',
     $               '------------------------------------------',
     $               '------------------------------------------',
     $               '------------------------------------------'
      Do 20 I = 1, DNum
        Write(11,1000) Gap(I), y1RK1(I), y3RK1(I), yrRK1(I), y1RK2(I),
     $                 y3RK2(I), yrRK2(I), y1RK3(I), y3RK3(I), yrRK3(I)
  20  Continue
      Write(11,1020) '------',
     $               '------------------------------------------',
     $               '------------------------------------------',
     $               '------------------------------------------'

      End Program Ex0204

*-----------------------------------------------------------------------
*
      Subroutine RK1(y1, y3, yr, xMin, xMax, yMin, GridN, Gap)

      Implicit None

      Double Precision  y1, y3, yr, xMin, xMax, yMin, Gap
      Integer           GridN
*   Dummy
      Integer           I
*   Value for n+1, n+0
      Double Precision  xn1, yn1, xn0, yn0
*   Function
      Double Precision  Func, Ext
      External          Func, Ext
*   Internal Variable
      Double Precision  k1, k2
     
      xn0 = xMin
      yn0 = yMin
      Do 10 I = 1, GridN
        xn1 = DBLE(I) / DBLE(GridN) * xMax +
     $        (1.0D0 - DBLE(I) / DBLE(GridN)) * xMin
        k1 = Gap * Func(xn0, yn0)
        k2 = Gap * Func(xn0 + 0.5D0*Gap, yn0 + 0.5D0*k1)
        yn1 = yn0 + k2
        If(DAbs(xn1-1.0D0).LT.1.0D-10) y1 = yn1 - Ext(xn1)
        If(DAbs(xn1-3.0D0).LT.1.0D-10) y3 = yn1 - Ext(xn1)
        xn0 = xn1
        yn0 = yn1
  10  Continue

      Do 20 I = 1, GridN
        xn1 = DBLE(I) / DBLE(GridN) * xMin +
     $        (1.0D0 - DBLE(I) / DBLE(GridN)) * xMax
        k1 = -Gap * Func(xn0, yn0)
        k2 = -Gap * Func(xn0 - 0.5D0*Gap, yn0 + 0.5D0*k1)
        yn1 = yn0 + k2
        xn0 = xn1
        yn0 = yn1
  20  Continue
      yr = Ext(xn1) - yn1

      End Subroutine RK1

*-----------------------------------------------------------------------
*
      Subroutine RK2(y1, y3, yr, xMin, xMax, yMin, GridN, Gap)

      Implicit None

      Double Precision  y1, y3, yr, xMin, xMax, yMin, Gap
      Integer           GridN
*   Dummy
      Integer           I
*   Value for n+1, n+0
      Double Precision  xn1, yn1, xn0, yn0
*   Function
      Double Precision  Func, Ext
      External          Func, Ext
*   Internal Variable
      Double Precision  k1, k2, k3
     
      xn0 = xMin
      yn0 = yMin
      Do 10 I = 1, GridN
        xn1 = DBLE(I) / DBLE(GridN) * xMax +
     $        (1.0D0 - DBLE(I) / DBLE(GridN)) * xMin
        k1 = Gap * Func(xn0, yn0)
        k2 = Gap * Func(xn0 + 0.5D0*Gap, yn0 + 0.5D0*k1)
        k3 = Gap * Func(xn0 + Gap, yn0 - k1 + 2.0D0*k2)
        yn1 = yn0 + 1.0D0/6.0D0 * (k1 + 4.0D0*k2 + k3)
        If(DAbs(xn1-1.0D0).LT.1.0D-10) y1 = yn1 - Ext(xn1)
        If(DAbs(xn1-3.0D0).LT.1.0D-10) y3 = yn1 - Ext(xn1)
        xn0 = xn1
        yn0 = yn1
  10  Continue

      Do 20 I = 1, GridN
        xn1 = DBLE(I) / DBLE(GridN) * xMin +
     $        (1.0D0 - DBLE(I) / DBLE(GridN)) * xMax
        k1 = -Gap * Func(xn0, yn0)
        k2 = -Gap * Func(xn0 - 0.5D0*Gap, yn0 + 0.5D0*k1)
        k3 = -Gap * Func(xn0 - Gap, yn0 - k1 + 2.0D0*k2)
        yn1 = yn0 + 1.0D0/6.0D0 * (k1 + 4.0D0*k2 + k3)
        xn0 = xn1
        yn0 = yn1
  20  Continue
      yr = Ext(xn1) - yn1

      End Subroutine RK2

*-----------------------------------------------------------------------
*
      Subroutine RK3(y1, y3, yr, xMin, xMax, yMin, GridN, Gap)

      Implicit None

      Double Precision  y1, y3, yr, xMin, xMax, yMin, Gap
      Integer           GridN
*   Dummy
      Integer           I
*   Value for n+1, n+0
      Double Precision  xn1, yn1, xn0, yn0
*   Function
      Double Precision  Func, Ext
      External          Func, Ext
*   Internal Variable
      Double Precision  k1, k2, k3, k4
     
      xn0 = xMin
      yn0 = yMin
      Do 10 I = 1, GridN
        xn1 = DBLE(I) / DBLE(GridN) * xMax +
     $        (1.0D0 - DBLE(I) / DBLE(GridN)) * xMin
        k1 = Gap * Func(xn0, yn0)
        k2 = Gap * Func(xn0 + 0.5D0*Gap, yn0 + 0.5D0*k1)
        k3 = Gap * Func(xn0 + 0.5D0*Gap, yn0 + 0.5D0*k2)
        k4 = Gap * Func(xn0 + Gap, yn0 + k3)
        yn1 = yn0 + 1.0D0/6.0D0 * (k1 + 2.0D0*k2 + 2.0D0*k3 + k4)
        If(DAbs(xn1-1.0D0).LT.1.0D-10) y1 = yn1 - Ext(xn1)
        If(DAbs(xn1-3.0D0).LT.1.0D-10) y3 = yn1 - Ext(xn1)
        xn0 = xn1
        yn0 = yn1
  10  Continue

      Do 20 I = 1, GridN
        xn1 = DBLE(I) / DBLE(GridN) * xMin +
     $        (1.0D0 - DBLE(I) / DBLE(GridN)) * xMax
        k1 = -Gap * Func(xn0, yn0)
        k2 = -Gap * Func(xn0 - 0.5D0*Gap, yn0 + 0.5D0*k1)
        k3 = -Gap * Func(xn0 - 0.5D0*Gap, yn0 + 0.5D0*k2)
        k4 = -Gap * Func(xn0 - Gap, yn0 + k3)
        yn1 = yn0 + 1.0D0/6.0D0 * (k1 + 2.0D0*k2 + 2.0D0*k3 + k4)
        xn0 = xn1
        yn0 = yn1
  20  Continue
      yr = Ext(xn1) - yn1

      End Subroutine RK3

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
