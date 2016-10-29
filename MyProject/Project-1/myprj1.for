*----------------------------------------------------------------------_
*
      Program MyPrj1
*
*     Purpose: Solve the vibration problem
*
      Implicit None
*
*-- Common Definitions
*
*   Write file index
      Integer           IOut
      Common /IO/       IOut
*   Calculation Parameter
      Double Precision  e, rScale, rMin
      Common /CalPar/   e, rScale, rMin
*   Constants
      Double Precision  Pi
      Common /Const/    Pi
*   Error Information
      Integer           IStat
      Common /ErStat/   IStat
*   Global Integral grid number
      Integer           IGridN
      Common /Grid/     IGridN
*   Zero Tolerence (Two values refers to the error or warning situation)
      Double Precision  Thr0, ThrW0
      Common /ZeroTh/   Thr0, ThrW0
*   Deviation Threahold
*!!   Non-FORTRAN77 syntax: An effort to make the variables stored into
*!!   the memories of the same length
      Integer*8         ItrMax
      Double Precision  DevX, DevY
      Common /DevThr/   ItrMax, DevX, DevY
*
*-- Internal Functions
*
      Double Precision  Itg5P1, Itg5P2, RtSec2, FI
      External          Itg5P1, ITg5P2, RtSec2, FI
      Integer           IFlr
      External          IFlr
*
*-- Internal Variables
*
*   Numbers of B should be calculated
      Integer           bNum
      Parameter         (bNum=40)
      Double Precision  bMax, bMin, bH
*   B as variable
      Double Precision  b
*   Stored Variables
      Double Precision  bV(bNum), TV(bNum), rMinV(bNum)
*   Dummy
      Integer           I
*
*-- Debug Functions
      Double Precision  BCrit, RCrit, FC, FIC1, FIC2, FCRt
      External          BCrit, RCrit, FC, FIC1, FIC2, FCRt
*
*-----------------------------------------------------------------------

*-- Initialize Common Variables

*   IO
      IOut = 11
*   CalPar
      e = 0.3D0
      RScale = 10.0D0
*   Const
      Pi = DACos(-1.0D0)
*   ErStat
      IStat = 0
*   Grid
      IGridN = 50000
*   ZeroTh
      Thr0 = 1.0D-7
      ThrW0 = 1.0D-10
*   DevThr
      ItrMax = 100
      DevX = 1.0D-10
      DevY = 1.0D-10

*-- File I/O 
      Open(IOut, File='myprj1.txt')
 1000 Format (1x, F10.4, F15.7, F16.7)
 1010 Format (1x, A10, A15, A16)

*-- Internal Initial variables
      bMin = 0.5D0
      bMax = 2.5D0
      bH = (bMax - bMin) / DBLE(bNum - 1)

*-- Main Program
      Do 10 I = 1, bNum
        b = bMin + DBLE(I-1) * bH
        bV(I) = b
        rMin = FCRt(b)
        rMinV(I) = rMin
        TV(I) = FI(b) / Pi * 1.8D2
  10  Continue

      Write(IOut,*) 
      Write(IOut,*) 'Impact Factor Minimum = ', bMin
      Write(IOut,*) 'Impact Factro Maximum = ', bMax
      Write(IOut,*) 'Energy                = ', e
      Write(IOut,*) 'Integral Upper Limit  = ', rScale
      Write(IOut,*) 
      Write(IOut,*) 'X Threshold           = ', DevX
      Write(IOut,*) 'Y Threshold           = ', DevY
      Write(IOut,*) 'Integral Grid Numbers = ', IGridN
      Write(IOut,*) 'Maximun Root Search   = ', ItrMax
      Write(IOut,*) 
      Write(IOut,*) 'bCrit                 = ', BCrit(e)
      Write(IOut,*) 'rCrit                 = ', RCrit(e)
      Write(IOut,*) 
      Write(IOut,1010) '----------', '---------------', 
     $                 '----------------'
      Write(IOut,1010) 'Impact', 'rMin', 'Theta'
      Write(IOut,1010) '----------', '---------------', 
     $                 '----------------'
      Do 20 I = 1, bNum
        Write(IOut,1000) bV(I), rMinV(I), TV(I)
  20  Continue
      Write(IOut,1010) '----------', '---------------', 
     $                 '----------------'

      Stop
      End Program MyPrj1

************************************************************************
*                                                                      *
*     USER DEFINE SUBROUTINES                                          *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*
      Function FI(b)
*
*     Purpose: The integral result (=Theta)
*
      Implicit None
*
*-- Common Variables
*
*   Calculation Parameter
      Double Precision  e, rScale, rMin
      Common /CalPar/   e, rScale, rMin
*   Global Integral grid number
      Integer           IGridN
      Common /Grid/     IGridN
*
*-- Variables
*
      Double Precision  FI, b
*
*-- Internal Variables
*
*   Integral Upper Limit
      Double Precision  rMax
*
*-- External Functions
*
*   Integral Core Function
      Double Precision  FIC1, FIC2
      External          FIC1, FIC2
*   Integral Routine
      Double Precision  Itg5P2
      External          Itg5P2
*   Root of FC (rMin)
      Double Precision  FCRt
      External          FCRt
*
*-----------------------------------------------------------------------

      rMax = rScale

      FI = Itg5P2(FIC2, b, b, rMax, IGridN) 
     $   - Itg5P2(FIC1, b, rMin, rMax, IGridN)

      Call Err
      Return
      End Function FI

*-----------------------------------------------------------------------
*
      Function FCRt(b)
*
*     Purpose: Find rMin (Biggest Root of FC)
*
      Implicit None
*
*-- Common Variables
*
*   Write file index
      Integer           IOut
      Common /IO/       IOut
*   Calculation Parameter
      Double Precision  e, rScale, rMin
      Common /CalPar/   e, rScale, rMin
*   Error Information
      Integer           IStat
      Common /ErStat/   IStat
*
*-- Variables
*
      Double Precision  FCRt, b
*
*-- External Functions
*
*   Critical Values BCrit, RCrit
      Double Precision  BCrit, RCrit
      External          BCrit, RCrit
*   Initial Guess for FC06
      Double Precision  FC06Rt
      External          FC06Rt
*   Initial Guess for FC02
      Double Precision  FC02Rt
*   Secant Search for Root
      Double Precision  RtSec2
      External          RtSec2
*   The Function We Wish to Find Root
      Double Precision  FC
      External          FC
*
*-----------------------------------------------------------------------

      FC02Rt = b

      If (e.LE.0.75D0) then
        If (b.GE.BCrit(e)) then
          FCRt = RtSec2(FC, b, 0.0D0, FC02Rt+0.1D-2, FC02Rt)
          If (FCRt.LT.RCrit(e)) then
            Write(IOut,*) 'Error Info in Function [FCRt]'
            Write(IOut,*) 'You Possibly find a small root rather than'
            Write(IOut,*) 'the greatest one.'
            IStat = -1
            Call Err
          End If
        Else
          FCRt = RtSec2(FC, b, 0.0D0, FC06Rt(e), 1.0D0)
        End If
      Else
        If (b.GT.2.0D0**(1.0D0/6.0D0)) then
          FCRt = RtSec2(FC, b, 0.0D0, FC02Rt+0.1D-2, FC02Rt)
          If (IStat.EQ.-1) then
            IStat = 0
            FCRt = RtSec2(FC, b, 0.0D0, FC06Rt(e), 1.0D0)
            If (IStat.EQ.-1) then
              Write(IOut,*) 'Error Info in Function [FCRt]'
              Write(IOut,*) 'Cannot Find Root for FC.'
              IStat = 2
              Call Err
            End If
          End If
        Else
          FCRt = RtSec2(FC, b, 0.0D0, FC06Rt(e), 1.0D0)
          If (IStat.EQ.-1) then
            IStat = 0
            FCRt = RtSec2(FC, b, 0.0D0, FC02Rt+0.1D-2, FC02Rt)
            If (IStat.EQ.-1) then
              Write(IOut,*) 'Error Info in Function [FCRt]'
              Write(IOut,*) 'Cannot Find Root for FC.'
              IStat = 2
              Call Err
            End If
          End If
        End If
      End If

      Return
      End Function FCRt

*-----------------------------------------------------------------------
*
      Function FIC1(b,r)
*
*     Purpose: The integral core (for the first function)
*
      Implicit None
*
*-- Common Variables
*
*   Calculation Parameter
      Double Precision  e, rScale, rMin
      Common /CalPar/   e, rScale, rMin
*   Deviation Threahold
*!!   Non-FORTRAN77 syntax: An effort to make the variables stored into
*!!   the memories of the same length
      Integer*8         ItrMax
      Double Precision  DevX, DevY
      Common /DevThr/   ItrMax, DevX, DevY
*
*-- Variables
*
      Double Precision  FIC1, b, r
      Double Precision  Tmp, FC, FCRt, FC1
      External          FC, FCRt, FC1
*
*-----------------------------------------------------------------------

      Tmp = FC(b,r)
      Call Test0(Tmp)
      If (Tmp.LE.DevX.OR.(r-rMin).LE.DevX) then
        Tmp = FC1(b,r)
      Else
        Tmp = Tmp / (r-rMin)
      End If
      FIC1 = (2.0D0 * b / (r*r) ) / DSqrt(Tmp)

      Return
      End Function FIC1

*-----------------------------------------------------------------------
*
      Function FIC2(b,r)
*
*     Purpose: The integral core (for the second function calculate Pi)
*
      Implicit None
*
*-- Variables
*
      Double Precision  FIC2, b, r
      Double Precision  Tmp, FC
      External          FC
*
*-----------------------------------------------------------------------

      Tmp = (r+b) / (r*r)
      FIC2 = (2.0D0 * b / (r*r) ) / DSqrt(Tmp)

      Return
      End Function FIC2

*-----------------------------------------------------------------------
*
      Function FC(b,r)
*
*     Purpose: Part of Integral Core 
*
      Implicit None
*
*-- Common Definitions
*
*   Calculation Parameter
      Double Precision  e, rScale, rMin
      Common /CalPar/   e, rScale, rMin
*
*-- Variables
*
      Double Precision  FC, b, r
*
*-----------------------------------------------------------------------

      FC = 1.0D0 - (b*b) / (r*r) 
     $     - (4.0D0 * r**(-12.0D0) - 4.0D0 * r**(-6.0D0)) / e

      Return
      End Function FC

*-----------------------------------------------------------------------
*
      Function FC1(b,r)
*
*     Purpose: 1st Derivative (by r) of FC
*
      Implicit None
*
*-- Common Definitions
*
*   Calculation Parameter
      Double Precision  e, rScale, rMin
      Common /CalPar/   e, rScale, rMin
*
*-- Variables
*
      Double Precision  FC1, b, r
*
*-----------------------------------------------------------------------

      FC1 = 2.0D0 * (b*b) / (r*r*r) 
     $     + (48.0D0 * r**(-13.0D0) - 24.0D0 * r**(-7.0D0)) / e

      Return
      End Function FC1

*-----------------------------------------------------------------------
*
      Function FC06RT(e)
*
*     Purpose: Calculate Positive Root of 1-4/e(r^12-r^6)==0
*
      Implicit None
*
*-- Variables
*
      Double Precision  FC06RT, e
*
*-----------------------------------------------------------------------

      FC06RT = (0.5D0 * (1.0D0 + DSQRT(1.0D0 + e))) ** (-1.0D0/6.0D0)

      End Function FC06RT

*-----------------------------------------------------------------------
*
      Function RCrit(e)
*
      Implicit None
*
*-- Variables
*
      Double Precision  RCrit, e
*
*-----------------------------------------------------------------------

      RCrit = ((4.0D0 + 2.0D0 * DSqrt(4.0D0 - 5.0D0 * e)) / e)
     $        ** (1.0D0/6.0D0)

      End Function RCrit

*-----------------------------------------------------------------------
*
      Function BCrit(e)
*
      Implicit None
*
*-- Variables
*
      Double Precision  BCrit, e
*
*-----------------------------------------------------------------------

      BCrit = (2.0D0**(2.0D0/3.0D0) * DSqrt(3.0D0) / 5.0D0) 
     $        * ((2.0D0 + DSqrt(4.0D0 - 5.0D0 * e)) / e) 
     $          ** (1.0D0/6.0D0)
     $        * DSqrt((2.0D0 - DSqrt(4.0D0 - 5.0D0 * e) 
     $          + 5.0D0 * e) / e)

      End Function BCrit

************************************************************************
*                                                                      *
*     Program Utilities                                                *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*
      Subroutine Err
*
*     Purpose: Call Err Message and Stop the Program
*
      Implicit None
*
*-- Common Definitions
*
*   Write file index
      Integer           IOut
      Common /IO/       IOut
*   Error Information
      Integer           IStat
      Common /ErStat/   IStat
*
*-----------------------------------------------------------------------

      If (IStat.GT.0) then
        Write(IOut,*) 'Error Return!'
        Write(IOut,*) 'Return Code: ', IStat
        Stop
      Else If (IStat.EQ.-1) then
        Write(IOut,*) 'Subroutine Error!'
        Write(IOut,*) 'Return Code: ', IStat
        Write(IOut,*) 'Continue running program!'
      End If

      End Subroutine

************************************************************************
*                                                                      *
*     MATHEMATICAL BACKGROUND UTILITIES                                *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*
      Function RtSec1(Func, Num, RtLo, RtUp)
*
*     Purpose: Find the root with secant method, where the first two
*              guesses are provided mannually.
*              The root-finding problem should be expressed as:
*                Func(Root) == Num
*
      Implicit None
*
*-- Common Blocks
*   Error Information
      Integer           IStat
      Common /ErStat/   IStat
*   Deviation Threahold
      Integer*8         ItrMax
      Double Precision  DevX, DevY
      Common /DevThr/   ItrMax, DevX, DevY
*
*-- Variables
*
      Double Precision  RtSec1, Func, Num, RtLo, RtUp
      External          Func
*
*   Func    The function we are concerning
*   Num     The right hand side of the equation
*   RtLo    The minimun guess
*   RtUp    The maximun guess
*
*   Note: Should have (Func(RtLo)-Num) * (Func(RtUp)-Num) < 0
*         or we probably have no or even roots inside within the maximun
*         guess or the minimun guess.
*
*-- Internal Variables
*
*   Iteration Dummy
      Integer           Itr
*   Internal maximun and minimun guess
      Double Precision  RtLoI, RtUpI
*   Internal Result
      Double Precision  FncLo, FncUp, FncSec
*
*-----------------------------------------------------------------------

*   Initialize internal guess
      RtLoI = RtLo
      RtUpI = RtUp
      Itr = 0
      FncLo = Func(RtLoI) - Num
      FncUp = Func(RtUpI) - Num

*   Find the root
  10  Continue
        Itr = Itr + 1
        RtSec1 = ( RtLoI * FncUp - RtUpI * FncLo ) / ( FncUp - FncLo )
        FncSec = Func(RtSec1) - Num
        RtLoI = RtUpI
        FncLo = FncUp
        RtUpI = RtSec1
        FncUp = FncSec
*   Converge Judgement
      If ( (Itr.LT.ItrMax).And.( DAbs(RtUpI-RtLoI).GT.DevX.Or.
     $      (DAbs(FncUp).GT.DevY.And.(FncUp-FncLo).GT.DevY))) then
*     Not converged
        Go To 10
      Else If ( (Itr.GE.ItrMax).And.( DAbs(RtUpI-RtLoI).GT.DevX.Or.
     $      (DAbs(FncUp).GT.DevY.And.(FncUp-FncLo).GT.DevY))) then
*     Not converged but reached the maximun iteration
        Print *, 'Error Info in Function [RtSec1]'
        Print *, 'Not converged but reached the maximun iteration'
        Print *, 'Last Iteration Information:'
        Print *, 'Upper Guess: ', RtUpI, '  Lower Guess: ', RtLoI, 
     $           'Upper Func: ', FncUp, '  Lower Guess: ', FncLo
        IStat = -1
      Else
        Return
      End If

      Call Err
      Return
      End Function RtSec1

*-----------------------------------------------------------------------
*
      Function RtSec2(Func, VM, Num, RtLo, RtUp)
*
*     Purpose: Find the root with secant method, where the first two
*              guesses are provided mannually.
*              The root-finding problem should be expressed as:
*                Func(VM, Root) == Num
*
      Implicit None
*
*-- Common Blocks
*   Error Information
      Integer           IStat
      Common /ErStat/   IStat
*   Deviation Threahold
      Integer*8         ItrMax
      Double Precision  DevX, DevY
      Common /DevThr/   ItrMax, DevX, DevY
*
*-- Variables
*
      Double Precision  RtSec2, Func, VM, Num, RtLo, RtUp
      External          Func
*
*   Func    The function we are concerning
*   Num     The right hand side of the equation
*   RtLo    The minimun guess
*   RtUp    The maximun guess
*
*   Note: Should have (Func(VM, RtLo)-Num) * (Func(VM, RtUp)-Num) < 0
*         or we probably have no or even roots inside within the maximun
*         guess or the minimun guess.
*
*-- Internal Variables
*
*   Iteration Dummy
      Integer           Itr
*   Internal maximun and minimun guess
      Double Precision  RtLoI, RtUpI
*   Internal Result
      Double Precision  FncLo, FncUp, FncSec
*
*-----------------------------------------------------------------------

*   Initialize internal guess
      RtLoI = RtLo
      RtUpI = RtUp
      Itr = 0
      FncLo = Func(VM, RtLoI) - Num
      FncUp = Func(VM, RtUpI) - Num

*   Find the root
  10  Continue
        Itr = Itr + 1
        RtSec2 = ( RtLoI * FncUp - RtUpI * FncLo ) / ( FncUp - FncLo )
        FncSec = Func(VM, RtSec2) - Num
        RtLoI = RtUpI
        FncLo = FncUp
        RtUpI = RtSec2
        FncUp = FncSec
*   Converge Judgement
      If ( (Itr.LT.ItrMax).And.( DAbs(RtUpI-RtLoI).GT.DevX.Or.
     $      (DAbs(FncUp).GT.DevY.And.(FncUp-FncLo).GT.DevY))) then
*     Not converged
        Go To 10
      Else If ( (Itr.GE.ItrMax).And.( DAbs(RtUpI-RtLoI).GT.DevX.Or.
     $      (DAbs(FncUp).GT.DevY.And.(FncUp-FncLo).GT.DevY))) then
*     Not converged but reached the maximun iteration
        Print *, 'Error Info in Function [RtSec2]'
        Print *, 'Not converged but reached the maximun iteration'
        Print *, 'Last Iteration Information:'
        Print *, 'Upper Guess: ', RtUpI, '  Lower Guess: ', RtLoI, 
     $           'Upper Func: ', FncUp, '  Lower Guess: ', FncLo
        IStat = -1
      Else
        Return
      End If

      Call Err
      Return
      End Function RtSec2

*-----------------------------------------------------------------------
*
      Function Itg5P1(Func, IntgLo, IntgUp, IGridN)
*
*     Purpose: Integral using Bode's method (Func(Var))
*
*-- Variables
*
      Implicit None
*
      Double Precision  Itg5P1, Func, IntgLo, IntgUp
      Integer           IGridN
      External          Func
*
*   Func      The function integrated
*   IntgUp    Upper limit of the integral
*   InteLo    Lower limit of the integral
*
*-- Internal Variables
*
*   Decleared elsewere
      Integer           IFlr
      External          IFlr
*   Dummy
      Integer           I
*   Iteration Number
      Integer           Itr
*   Grid Seperation
      Double Precision  H
*
*-----------------------------------------------------------------------

*   Redefining Integral grids and make the separation
      Itr = IFlr(DBLE(IGridN) / 4.0D0)
      H = (IntgUp - IntgLo) / (4.0D0 * DBLE(Itr))

*   Initial and final point calculation
      Itg5P1 = 7.0D0 * (Func(IntgLo) - Func(IntgUp))

*   Internal grid point calculation
      Do 10 I = 1, Itr
        Itg5P1 = Itg5P1 + 32.0D0 * Func(IntgLo + DBLE(I*4-3) * H)
     $                  + 12.0D0 * Func(IntgLo + DBLE(I*4-2) * H)
     $                  + 32.0D0 * Func(IntgLo + DBLE(I*4-1) * H)
     $                  + 14.0D0 * Func(IntgLo + DBLE(I*4  ) * H)
  10  Continue

*   Coefficient multiplication
      Itg5P1 = Itg5P1 * H * 2.0D0 / 45.0D0

      Call Err
      Return
      End Function 

*-----------------------------------------------------------------------
*
C     Function Itg5P2(Func, VM, IntgLo, IntgUp, IGridN)
*
*     Purpose: Integral using Bode's method (Func(VM, VarMain,VarInt))
*
*-- Variables
*
C     Implicit None
*
C     Double Precision  Itg5P2, Func, IntgLo, IntgUp, VM
C     Integer           IGridN
C     External          Func
*
*   Func      The function integrated
*   IntgUp    Upper limit of the integral
*   InteLo    Lower limit of the integral
*   VM        Main Variable but not used in the integral
*
*-- Internal Variables
*
*   Decleared elsewere
C     Integer           IFlr
C     External          IFlr
*   Dummy
C     Integer           I
*   Iteration Number
C     Integer           Itr
*   Grid Seperation
C     Double Precision  H
*
*-----------------------------------------------------------------------

*   Redefining Integral grids and make the separation
C     Itr = IFlr(DBLE(IGridN) / 4.0D0)
C     H = (IntgUp - IntgLo) / (4.0D0 * DBLE(Itr))

*   Initial and final point calculation
C     Itg5P2 = 7.0D0 * (Func(VM, IntgLo) - Func(VM, IntgUp))

*   Internal grid point calculation
C     Do 10 I = 1, Itr
C       Itg5P2 = Itg5P2 + 32.0D0 * Func(VM, IntgLo + DBLE(I*4-3) * H)
C    $                  + 12.0D0 * Func(VM, IntgLo + DBLE(I*4-2) * H)
C    $                  + 32.0D0 * Func(VM, IntgLo + DBLE(I*4-1) * H)
C    $                  + 14.0D0 * Func(VM, IntgLo + DBLE(I*4  ) * H)
C 10  Continue

*   Coefficient multiplication
C     Itg5P2 = Itg5P2 * H * 2.0D0 / 45.0D0

C     Call Err
C     Return
C     End Function 

*-----------------------------------------------------------------------
*
      Function Itg5P2(Func, VM, IntgLo, IntgUp, IGridN)
*
*     Purpose: Integral using Bode's method (Func(VM, VarMain,VarInt))
*
*-- Variables
*
      Implicit None
*
      Double Precision  Itg5P2, Func, IntgLo, IntgUp, VM
      Integer           IGridN
      External          Func
*
*   Func      The function integrated
*   IntgUp    Upper limit of the integral
*   InteLo    Lower limit of the integral
*   VM        Main Variable but not used in the integral
*
*-- Internal Variables
*
*   Decleared elsewere
      Integer           IFlr
      External          IFlr
*   Dummy
      Integer           I
*   Iteration Number
      Integer           Itr
*   Grid Seperation
      Double Precision  H
*
*-----------------------------------------------------------------------

*   Redefining Integral grids and make the separation
      Itr = IFlr(DBLE(IGridN) / 4.0D0)
      H = DSqrt(IntgUp - IntgLo) / (4.0D0 * DBLE(Itr))

*   Initial and final point calculation
      Itg5P2 = 7.0D0 * (- 2.0D0 * DSqrt(IntgUp-IntgLo)*Func(VM, IntgUp))

*   Internal grid point calculation
      Do 10 I = 1, Itr
        Itg5P2 = Itg5P2 + 32.0D0 * Func(VM, (DBLE(I*4-3) * H)**2.0D0 +
     $                    IntgLo)
     $                  + 12.0D0 * Func(VM, (DBLE(I*4-2) * H)**2.0D0 +
     $                    IntgLo)
     $                  + 32.0D0 * Func(VM, (DBLE(I*4-1) * H)**2.0D0 +
     $                    IntgLo)
     $                  + 14.0D0 * Func(VM, (DBLE(I*4  ) * H)**2.0D0 +
     $                    IntgLo)
  10  Continue

*   Coefficient multiplication
      Itg5P2 = Itg5P2 * H * 4.0D0 / 45.0D0

      Call Err
      Return
      End Function 

*-----------------------------------------------------------------------
*
      Function Dv15P1(Func, Var, Dev)
*
*     Purpose: Calculate the first derviation of a function with Boole's
*              rule (For only one variable function)
*
      Implicit None
*
*-- Variables
      Double Precision  Dv15P1, Func, Var, Dev
      External          Func
*
*   Dev     The seperation of the calculation points
*
*-----------------------------------------------------------------------

      Dv15P1 = (Func(Var-2.0D0*Dev) - 8.0D0 * Func(Var-Dev)
     $        - Func(Var+2.0D0*Dev) + 8.0D0 * Func(Var+Dev))
     $         / (12.0D0 * Dev)

      End Function Dv15P1

*-----------------------------------------------------------------------
*
      Function Dv25P1(Func, Var, Dev)
*
*     Purpose: Calculate the second derviation of a function with 
*              Boole's rule (For only one variable function)
*
      Implicit None
*
*-- Variables
      Double Precision  Dv25P1, Func, Var, Dev
      External          Func
*
*   Dev     The seperation of the calculation points
*
*-----------------------------------------------------------------------

      Dv25P1 = (- Func(Var-2.0D0*Dev) + 16.0D0 * Func(Var-Dev)
     $          - Func(Var+2.0D0*Dev) + 16.0D0 * Func(Var+Dev)
     $          - 30.0D0 * Func(Var)) / (12.0D0 * Dev * Dev)

      End Function Dv25P1

*-----------------------------------------------------------------------
*      
      Function IFlr(Var)
*
*     Purpose: Evaluate the floor of a double precision float number
*
*-- Variables 
*
      Implicit None
*
      Integer           IFlr
      Double Precision  Var
*
*-----------------------------------------------------------------------

      IFlr = Int(Var)

      If (Var.LT.0.0D0) IFlr = IFlr - 1
      
      Return
      End Function IFlr

*-----------------------------------------------------------------------
*
      Subroutine Test0(Var)
*
      Implicit None
*
*-- Common Definitions
*
*   Write file index
      Integer           IOut
      Common /IO/       IOut
*   Error Information
      Integer           IStat
      Common /ErStat/   IStat
*   Zero Tolerence (Two values refers to the error or warning situation)
      Double Precision  Thr0, ThrW0
      Common /ZeroTh/   Thr0, ThrW0
*
*-- Variables
      Double Precision  Var
*
*-----------------------------------------------------------------------

      If (Var.LT.0) then
        If (-Var.LE.ThrW0) then
          Var = 0.0D0
        Else If (-Var.LE.Thr0) then
      Write(IOut,*)'Tested variable less than 0 but larger than 1.0D-15'
          Write(IOut,*) 'Variable value: ', Var
          Var = 0.0D0
        Else
          Write(IOut,*) 'Tested variable less than 0. Error return.'
          Write(IOut,*) 'Variable value: ', Var
          IStat = 1
          Return
        End If
      End If

      Call Err
      Return
      End Subroutine Test0

