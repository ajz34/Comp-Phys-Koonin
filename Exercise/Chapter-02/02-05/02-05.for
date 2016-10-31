      Program Ex0205
      
      Implicit None

*   Pi
      Double Precision  Pi
      Common /CalPar/   Pi

*   Grid Number and Record Number
      Integer           NGrid
      Parameter         (NGrid = 10000)
      Integer           NRec
      Parameter         (NRec = 50)
      Integer           NSep
      Parameter         (NSep = (NGrid - Mod(NGrid,NRec)) / NRec)
*   Dummy
      Integer           I, J
*   Initial Value and Calculation Range
      Double Precision  tMin, tMax, yMin, pMin
      Common /InitV/    tMin, tMax, yMin, pMin
*   Variable Storage
      Double Precision  t(NGrid), y(NGrid), p(NGrid), yExt(NGrid), 
     $                  pExt(NGrid)
*   t Axis Grid Separation
      Double Precision  H
*   Internal Variables
      Double Precision  k1, k2, k3, k4
*   For Initial Calculation
      Double Precision  y0, p0, t0
*   Functions
      Double Precision  Fdy, Fdp, FyExt, FpExt
      External          Fdy, Fdp, FyExt, FpExt

      Pi = DACos(-1.0D0)
      tMin = 0.0D0
      tMax = 10.0D0
      yMin = 1.5D0
      pMin = 0.5D0

      Open(11, File = '02-05.txt')
 1000 Format(F6.3, 6F14.7)
 1010 Format(A6,6A14)

*-- Go Forward

      t(1) = (1.0D0 - DBLE(1)/DBLE(NGrid)) * tMin + 
     $       DBLE(1)/DBLE(NGrid) * tMax
      t0 = tMin
      H = t(1) - t0
      y0 = yMin
      p0 = pMin
      
      k1 = H * Fdp(t0, y0, p0)
      k2 = H * Fdp(t0 + 0.5D0*H, y0, p0 + 0.5D0*k1)
      k3 = H * Fdp(t0 + 0.5D0*H, y0, p0 + 0.5D0*k2)
      k4 = H * Fdp(t0 + H, y0, p0 + k3)
      p(1) = p0 + 1.0D0/6.0D0* (k1 + 2.0D0*k2 + 2.0D0*k3 + k4)
      p0 = p(1)
      k1 = H * Fdy(t0, y0, p0)
      k2 = H * Fdy(t0 + 0.5D0*H, y0 + 0.5D0*k1, p0)
      k3 = H * Fdy(t0 + 0.5D0*H, y0 + 0.5D0*k2, p0)
      k4 = H * Fdy(t0 + H, y0 + k3, p0)
      y(1) = y0 + 1.0D0/6.0D0* (k1 + 2.0D0*k2 + 2.0D0*k3 + k4)
      y0 = y(1)
      t0 = t(1)

      Do 10 I = 2, NGrid
        t(I) = (1.0D0 - DBLE(I)/DBLE(NGrid)) * tMin + 
     $         DBLE(I)/DBLE(NGrid) * tMax
        H = t(I) - t0
      
        k1 = H * Fdp(t0, y0, p0)
        k2 = H * Fdp(t0 + 0.5D0*H, y0, p0 + 0.5D0*k1)
        k3 = H * Fdp(t0 + 0.5D0*H, y0, p0 + 0.5D0*k2)
        k4 = H * Fdp(t0 + H, y0, p0 + k3)
        p(I) = p0 + 1.0D0/6.0D0* (k1 + 2.0D0*k2 + 2.0D0*k3 + k4)
        p0 = p(I)
        k1 = H * Fdy(t0, y0, p0)
        k2 = H * Fdy(t0 + 0.5D0*H, y0 + 0.5D0*k1, p0)
        k3 = H * Fdy(t0 + 0.5D0*H, y0 + 0.5D0*k2, p0)
        k4 = H * Fdy(t0 + H, y0 + k3, p0)
        y(I) = y0 + 1.0D0/6.0D0* (k1 + 2.0D0*k2 + 2.0D0*k3 + k4)
        y0 = y(I)
        t0 = t(I)
  10  Continue
      
      Do 50 I = 1, NGrid
        yExt(I) = FyExt(t(I))
        pExt(I) = FpExt(t(I))
  50  Continue

      Write(11,*)
      Write(11,'(38x,A10)') 'Go Forward'
      Write(11,*)
      Write(11,1010) '------', 
     $               '--------------','--------------','--------------',
     $               '--------------','--------------','--------------'
      Write(11,1010) 't  ', 'y(t)   ',  'y Exact ', 'y Diff ',
     $               'p(t)   ',  'p Exact ', 'p Diff '
      Write(11,1010) '------', 
     $               '--------------','--------------','--------------',
     $               '--------------','--------------','--------------'
      Do 20 J = 1, NRec
        I = J * NSep
        Write(11,1000) t(I), y(I), yExt(I), y(I)-yExt(I),
     $                 p(I), pExt(I), p(I)-pExt(I)
  20  Continue
      Write(11,1010) '------', 
     $               '--------------','--------------','--------------',
     $               '--------------','--------------','--------------'

*-- Return Back

      t(1) = (1.0D0 - DBLE(1)/DBLE(NGrid)) * tMax + 
     $       DBLE(1)/DBLE(NGrid) * tMin
      H = t(1) - t0
      
      k1 = H * Fdp(t0, y0, p0)
      k2 = H * Fdp(t0 + 0.5D0*H, y0, p0 + 0.5D0*k1)
      k3 = H * Fdp(t0 + 0.5D0*H, y0, p0 + 0.5D0*k2)
      k4 = H * Fdp(t0 + H, y0, p0 + k3)
      p(1) = p0 + 1.0D0/6.0D0* (k1 + 2.0D0*k2 + 2.0D0*k3 + k4)
      p0 = p(1)
      k1 = H * Fdy(t0, y0, p0)
      k2 = H * Fdy(t0 + 0.5D0*H, y0 + 0.5D0*k1, p0)
      k3 = H * Fdy(t0 + 0.5D0*H, y0 + 0.5D0*k2, p0)
      k4 = H * Fdy(t0 + H, y0 + k3, p0)
      y(1) = y0 + 1.0D0/6.0D0* (k1 + 2.0D0*k2 + 2.0D0*k3 + k4)
      y0 = y(1)
      t0 = t(1)

      Do 30 I = 2, NGrid
        t(I) = (1.0D0 - DBLE(I)/DBLE(NGrid)) * tMax + 
     $         DBLE(I)/DBLE(NGrid) * tMin
        H = t(I) - t0
      
        k1 = H * Fdp(t0, y0, p0)
        k2 = H * Fdp(t0 + 0.5D0*H, y0, p0 + 0.5D0*k1)
        k3 = H * Fdp(t0 + 0.5D0*H, y0, p0 + 0.5D0*k2)
        k4 = H * Fdp(t0 + H, y0, p0 + k3)
        p(I) = p0 + 1.0D0/6.0D0* (k1 + 2.0D0*k2 + 2.0D0*k3 + k4)
        p0 = p(I)
        k1 = H * Fdy(t0, y0, p0)
        k2 = H * Fdy(t0 + 0.5D0*H, y0 + 0.5D0*k1, p0)
        k3 = H * Fdy(t0 + 0.5D0*H, y0 + 0.5D0*k2, p0)
        k4 = H * Fdy(t0 + H, y0 + k3, p0)
        y(I) = y0 + 1.0D0/6.0D0* (k1 + 2.0D0*k2 + 2.0D0*k3 + k4)
        y0 = y(I)
        t0 = t(I)
  30  Continue
      
      Do 60 I = 1, NGrid
        yExt(I) = FyExt(t(I))
        pExt(I) = FpExt(t(I))
  60  Continue

      Write(11,*)
      Write(11,'(38x,A11)') 'Return Back'
      Write(11,*)
      Write(11,1010) '------', 
     $               '--------------','--------------','--------------',
     $               '--------------','--------------','--------------'
      Write(11,1010) 't  ', 'y(t)   ',  'y Exact ', 'y Diff ',
     $               'p(t)   ',  'p Exact ', 'p Diff '
      Write(11,1010) '------', 
     $               '--------------','--------------','--------------',
     $               '--------------','--------------','--------------'
      Do 40 J = 1, NRec
        I = J * NSep
        Write(11,1000) t(I), y(I), yExt(I), y(I)-yExt(I),
     $                 p(I), pExt(I), p(I)-pExt(I)
  40  Continue
      Write(11,1010) '------', 
     $               '--------------','--------------','--------------',
     $               '--------------','--------------','--------------'

      End Program Ex0205

*-----------------------------------------------------------------------
*
      Function Fdy(t, y, p)

      Implicit None

      Double Precision  Fdy, t, y, p

      Fdy = p

      End Function Fdy

*-----------------------------------------------------------------------
*
      Function Fdp(t, y, p)

      Implicit None

*   Pi
      Double Precision  Pi
      Common /CalPar/   Pi

      Double Precision  Fdp, t, y, p

      Fdp = - 4.0D0 * Pi*Pi * y

      End Function Fdp

*-----------------------------------------------------------------------
*
      Function FyExt(t)

      Implicit None

*   Pi
      Double Precision  Pi
      Common /CalPar/   Pi
*   Initial Value and Calculation Range
      Double Precision  tMin, tMax, yMin, pMin
      Common /InitV/    tMin, tMax, yMin, pMin

      Double Precision  FyExt, t
 
      FyExt = yMin * DCos(2.0D0*Pi * t) + pMin * DSin(2.0D0*Pi * t) /
     $        (2.0D0*Pi)

      End Function FyExt

*-----------------------------------------------------------------------
*
      Function FpExt(t)

      Implicit None

*   Pi
      Double Precision  Pi
      Common /CalPar/   Pi
*   Initial Value and Calculation Range
      Double Precision  tMin, tMax, yMin, pMin
      Common /InitV/    tMin, tMax, yMin, pMin

      Double Precision  FpExt, t

      FpExt = pMin * DCos(2.0D0*Pi * t) - 2.0D0*Pi * DSin(2.0D0*Pi * t)
     $        * yMin

      End Function FpExt

