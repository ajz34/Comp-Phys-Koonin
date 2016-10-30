      Program Ex0202

      Implicit None

*   Guess Number
      Integer           N
      Parameter         (N = 30)
*   Guess Minimum
      Double Precision  xMin
      Parameter         (xMin = 0.0D0)
*   Guess Maximum
      Double Precision  xMax
      Parameter         (xMax = 3.0D0)
*   Y Value of Guess Minimum
      Double Precision  yMin
      Parameter         (yMin = 1.0D0)
*   X Interval
      Double Precision  H
      Parameter         (H = (xMax - xMin) / Dble(N))
*   Data Storage
      Double Precision  y(N), x(N), yExt(N), yDif(N), yRat(N)
      Double Precision  y0, xTmp, yTmp
*   Dummy
      Integer           I
*   Diff Func
      Double Precision  Func
      External          Func

      Open(11, File = '02-02.txt')
 1000 Format(5F14.7)
 1010 Format(5A14)
 1020 Format(A50)

*   Initial Value
      Do 10 I = 1, N
        x(I) = Dble(I)/Dble(N) * xMax + 
     $         (Dble(1.0D0) - Dble(I) / Dble(N)) * xMin
        yExt(I) = DExp(-x(I)*x(I)/2.0D0)
  10  Continue

*-- Euler's Method

      y0 = yMiN
      y(1) = y0 + H * Func(xMin,y0)
      yDif(1) = y(1) - yExt(1) 
      yRat(1) = yDif(1) / yExt(1)

      Do 20 I = 2, N
        y(I) = y(I-1) + H * Func(x(I-1), y(I-1))
        yDif(I) = y(I) - yExt(I)
        yRat(I) = yDif(I) / yExt(I)
  20  Continue

      Write(11,*)
      Write(11,1020) 'Euler Method         '
      Write(11,*)
      Write(11,1010) '--------------','--------------','--------------',
     $               '--------------','--------------'
      Write(11,1010) 'x', 'Exact y', 'Calc y', 'Error Value', 
     $               'Error Ratio'
      Write(11,1010) '--------------','--------------','--------------',
     $               '--------------','--------------'
      Do 30 I = 1, N
        Write(11,1000) x(I), yExt(I), y(I), yDif(I), yRat(I)
  30  Continue
      Write(11,1010) '--------------','--------------','--------------',
     $               '--------------','--------------'
      Write(11,*)

      yTmp = y(N)
      Do 40 I = 1, N
        xTmp = Dble(I)/Dble(N) * xMin + H +
     $         (Dble(1.0D0) - Dble(I) / Dble(N)) * xMax
        yTmp = yTmp - H * Func(xTmp, yTmp)
  40  Continue
      Write(11,'(A18,F14.7)') 'Recursion Value: ', yTmp
      Write(11,'(A18,F14.7)') 'Recursion Error: ', yMin - yTmp
      Write(11,*)

*-- Adams-Bashforth's Method (1th extrapolate)

      Do 110 I = 1, N
        x(I) = Dble(I)/Dble(N) * xMax + 
     $         (Dble(1.0D0) - Dble(I) / Dble(N)) * xMin
        yExt(I) = DExp(-x(I)*x(I)/2.0D0)
  110 Continue

      y0 = yMiN
      y(1) = y0 + H * Func(xMin,y0)
      yDif(1) = y(1) - yExt(1) 
      yRat(1) = yDif(1) / yExt(1)

      Do 120 I = 2, 2
        y(I) = y(I-1) + H * Func(x(I-1), y(I-1))
        yDif(I) = y(I) - yExt(I)
        yRat(I) = yDif(I) / yExt(I)
  120 Continue

      Do 130 I = 3, N
        y(I) = y(I-1) + H * (1.5D0 * Func(x(I-1), y(I-1)) - 
     $         0.5D0 * Func(x(I-2), y(I-2)))
        yDif(I) = y(I) - yExt(I)
        yRat(I) = yDif(I) / yExt(I)
  130 Continue

      Write(11,*)
      Write(11,1020) '1th Adams-Bashforth Method  '
      Write(11,*)
      Write(11,1010) '--------------','--------------','--------------',
     $               '--------------','--------------'
      Write(11,1010) 'x', 'Exact y', 'Calc y', 'Error Value', 
     $               'Error Ratio'
      Write(11,1010) '--------------','--------------','--------------',
     $               '--------------','--------------'
      Do 140 I = 1, N
        Write(11,1000) x(I), yExt(I), y(I), yDif(I), yRat(I)
  140 Continue
      Write(11,1010) '--------------','--------------','--------------',
     $               '--------------','--------------'
      Write(11,*)

      Do 150 I = 1, N
        x(I) = Dble(I)/Dble(N) * xMin + 
     $         (Dble(1.0D0) - Dble(I) / Dble(N)) * xMax
        yExt(I) = DExp(-x(I)*x(I)/2.0D0)
  150 Continue
 
      y0 = y(N)
      y(1) = y0 - H * Func(xMax,y0)
      yDif(1) = y(1) - yExt(1) 
      yRat(1) = yDif(1) / yExt(1)
      
      Do 160 I = 2, 2
        y(I) = y(I-1) - H * Func(x(I-1), y(I-1))
        yDif(I) = y(I) - yExt(I)
        yRat(I) = yDif(I) / yExt(I)
  160 Continue
 
      Do 170 I = 3, N
        y(I) = y(I-1) - H * (1.5D0 * Func(x(I-1), y(I-1)) - 
     $         0.5D0 * Func(x(I-2), y(I-2)))
        yDif(I) = y(I) - yExt(I)
        yRat(I) = yDif(I) / yExt(I)
  170 Continue
      
      Write(11,'(A18,F14.7)') 'Recursion Value: ', y(N)
      Write(11,'(A18,F14.7)') 'Recursion Error: ', yMin - y(N)
      Write(11,*)

*-- Adams-Bashforth's Method (3th extrapolate)

      Do 210 I = 1, N
        x(I) = Dble(I)/Dble(N) * xMax + 
     $         (Dble(1.0D0) - Dble(I) / Dble(N)) * xMin
        yExt(I) = DExp(-x(I)*x(I)/2.0D0)
  210 Continue

      y0 = yMiN
      y(1) = y0 + H * Func(xMin,y0)
      yDif(1) = y(1) - yExt(1) 
      yRat(1) = yDif(1) / yExt(1)

      Do 220 I = 2, 4
        y(I) = y(I-1) + H * Func(x(I-1), y(I-1))
        yDif(I) = y(I) - yExt(I)
        yRat(I) = yDif(I) / yExt(I)
  220 Continue

      Do 230 I = 5, N
        y(I) = y(I-1) + H / 24.0D0 * ( 55.0D0 * Func(x(I-1),y(I-1)) - 
     $  59.0D0 * Func(x(I-2),y(I-2)) + 37.0D0 * Func(x(I-3),y(I-3)) - 
     $   9.0D0 * Func(x(I-4),y(I-4)) )
        yDif(I) = y(I) - yExt(I)
        yRat(I) = yDif(I) / yExt(I)
  230 Continue

      Write(11,*)
      Write(11,1020) '3th Adams-Bashforth Method  '
      Write(11,*)
      Write(11,1010) '--------------','--------------','--------------',
     $               '--------------','--------------'
      Write(11,1010) 'x', 'Exact y', 'Calc y', 'Error Value', 
     $               'Error Ratio'
      Write(11,1010) '--------------','--------------','--------------',
     $               '--------------','--------------'
      Do 240 I = 1, N
        Write(11,1000) x(I), yExt(I), y(I), yDif(I), yRat(I)
  240 Continue
      Write(11,1010) '--------------','--------------','--------------',
     $               '--------------','--------------'
      Write(11,*)

      Do 250 I = 1, N
        x(I) = Dble(I)/Dble(N) * xMin + 
     $         (Dble(1.0D0) - Dble(I) / Dble(N)) * xMax
        yExt(I) = DExp(-x(I)*x(I)/2.0D0)
  250 Continue
 
      y0 = y(N)
      y(1) = y0 - H * Func(xMax,y0)
      yDif(1) = y(1) - yExt(1) 
      yRat(1) = yDif(1) / yExt(1)
      
      Do 260 I = 2, 4
        y(I) = y(I-1) - H * Func(x(I-1), y(I-1))
        yDif(I) = y(I) - yExt(I)
        yRat(I) = yDif(I) / yExt(I)
  260 Continue
 
      Do 270 I = 5, N
        y(I) = y(I-1) - H / 24.0D0 * ( 55.0D0 * Func(x(I-1),y(I-1)) - 
     $  59.0D0 * Func(x(I-2),y(I-2)) + 37.0D0 * Func(x(I-3),y(I-3)) - 
     $   9.0D0 * Func(x(I-4),y(I-4)) )
        yDif(I) = y(I) - yExt(I)
        yRat(I) = yDif(I) / yExt(I)
  270 Continue
     
      Write(11,'(A18,F14.7)') 'Recursion Value: ', y(N)
      Write(11,'(A18,F14.7)') 'Recursion Error: ', yMin - y(N)
      Write(11,*)

      End Program Ex0202

*-----------------------------------------------------------------------
*
      Function Func(x,y)

      Implicit None

      Double Precision  Func, x, y

      Func = -x*y

      Return
      End Function Func

