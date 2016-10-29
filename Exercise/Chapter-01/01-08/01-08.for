*----------------------------------------------------------------------_
*
      Program MyEmp1
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
      Double Precision  Gmr
      Common /CalPar/   Gmr
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
      Double Precision  Itg5P1, Itg5P2, RtSec, Func1
      External          Itg5P1, ITg5P2, RtSec, Func1
      Integer           IFlr
      External          IFlr
*
*-- Internal Variables
*
*   Maximun Quantum Number (Related with DevX)
      Integer           IQntN
*
*-----------------------------------------------------------------------

*-- Initialize Common Variables

*   IO
      IOut = 11
*   CalPar
      Gmr = 50.0D0
*   Const
      Pi = DACos(-1.0D0)
*   ErStat
      IStat = 0
*   Grid
      IGridN = 100000
*   ZeroTh
      Thr0 = 1.0D-7
      ThrW0 = 1.0D-10
*   DevThr
      ItrMax = 150
      DevX = 1.0D-10
      DevY = 1.0D-10

*-- File I/O 
      Open(IOut, File='01-08.txt')

*-- I/O Format
C1000 Format (1x,A50)
C1010 Format (1x,A50,A10,F15.7)
C1020 Format (1x,A50,A10,8x,I7)
C1030 Format (1x,A50,A10,D15.8)

*-- Find Maximun Quantum Number and Call the Subroutine Main
      IQntN = IFlr(Func1(-DevX)/Pi-0.5D0)
      If (IQntN.EQ.0) then 
        Print *, 'No acceptable quantum number! Exit Program!'
        Stop
      End If
      Call MainPg(IQntN)

      Stop
      End Program MyEmp1

*----------------------------------------------------------------------_
*
      Subroutine MainPg(IQntN)
*
*     Purpose: Main program for this problem
*
      Implicit None
*
*-- Common Definitions
*
*   Write file index
      Integer           IOut
      Common /IO/       IOut
*   Constants
      Double Precision  Pi
      Common /Const/    Pi
*   Deviation Threahold
*!!   Non-FORTRAN77 syntax: An effort to make the variables stored into
*!!   the memories of the same length
      Integer*8         ItrMax
      Double Precision  DevX, DevY
      Common /DevThr/   ItrMax, DevX, DevY
*
*-- Internal Functions
*
      Double Precision  Itg5P1, Itg5P2, RtSec, Func1, IFlr
      External          Itg5P1, ITg5P2, RtSec, Func1, IFlr
*
*-- Internal Variables
*
*   Dummy
      Integer           I
*   Maximun Quantum Number (Related with DevX)
      Integer           IQntN
*   Internal Root Guess 
      Double Precision  RtLoI, RtUpI
*   Internal Result Temporary Store
      Double Precision  RtI
*   Result
      Integer           IQntV(IQntN+1)
      Double Precision  EngV(IQntN+1), ItgLoV(IQntN+1), ItgUpV(IQntN+1)
*
*-----------------------------------------------------------------------

*-- I/O Format
 1000 Format(A47)
 1010 Format(A8,A15,2(A12))
 1020 Format(I8,F15.6,1(F12.6))
 1030 Format(I8,D15.5,1(F12.6))

*-- Set the first root search condition
      RtLoI = -1.0D0+DevX
      RtUpI = -DevX
 
      Do 10 I = 0, IQntN
*     Root Solve
        RtI = RtSec(Func1, (DBLE(I)+0.5D0)*Pi, RtLoI, RtUpI)
*     Result Write In
        IQntV(I+1)  = I
        EngV(I+1)   = RtI
        ItgLoV(I+1) = ( (DSqrt(1.0D0+RtI)-1.0D0) * 2.0D0 / RtI ) ** 
     $         (1.0D0 / 6.0D0)
        ItgUpV(I+1) = ( (-DSqrt(1.0D0+RtI)-1.0D0) * 2.0D0 / RtI ) **
     $         (1.0D0 / 6.0D0)
*     Prepare for the next solve
        RtLoI = RtI
        RtUpI = RtI * 1.0D-5
        If (RtUpI.GT.-DevX) RtUpI = DevX
  10  Continue

      Write(IOut,1000) '-----------------------------------            '
      Write(IOut,1010) '   Level', '         Energy', '        XMin',
     $                 '            '
      Write(IOut,1000) '-----------------------------------            '
      Do 20 I = 0, IQntN
        If(-EngV(I+1).GE.1.0D-2) then
        Write(IOut,1020) IQntV(I+1), EngV(I+1),
     $                   (EngV(I+1)+1.0D0) / (EngV(1)+1.0D0)
        Else
        Write(IOut,1030) IQntV(I+1), EngV(I+1),
     $                   (EngV(I+1)+1.0D0) / (EngV(1)+1.0D0)
        End If
  20  Continue
      Write(IOut,1000) '-----------------------------------            '

      Call Ex0108

      End Subroutine MainPg

*-----------------------------------------------------------------------
*
      Subroutine Ex0108
*
      Implicit None
*
*-- Common Definitions
*
*   Write file index
      Integer           IOut
      Common /IO/       IOut
*   Calculation Parameter
      Double Precision  Gmr
      Common /CalPar/   Gmr
*
*-- Functions
      Double Precision  Dv25P1, FuncP
      External          Dv25P1, FuncP
*
*-- Internal Variables
*   Lowerest Point for the Potential (If Known as 2^(1/6))
      Double Precision  PotMin
*   Temporary Store
      Double Precision  Expt, Numr
*
*-----------------------------------------------------------------------

      PotMin = (2.0D0)**(1.0D0/6.0D0)

      Expt = DSqrt( Dv25P1(FuncP, PotMin, 1.0D-6) * 2.0D0 ) / Gmr
      Numr = 6.0D0 * 2.0D0 ** (5.0D0/6.0D0) / Gmr

      Write(IOut,*)
      Write(IOut,*)'Second Derviation of Potential at Potential Minimun'
     $           , Expt
      Write(IOut,*) 'Expected Value', Numr
      Write(IOut,*) 'Numerical / Expect', Expt / Numr

      End Subroutine Ex0108

************************************************************************
*                                                                      *
*     USER DEFINE SUBROUTINES                                          *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*
      Function Func1(Eng)
*
*     Purpose: The integral of the Gmr * (e - v(x))^0.5
*     Note:    This function actually has two variables: x, Gmr. Gmr is
*              pre-defined by the common blocks.
*
      Implicit None
*
*   Calculation Parameter
      Double Precision  Gmr
      Common /CalPar/   Gmr
*   Global Integral grid number
      Integer           IGridN
      Common /Grid/     IGridN
*
*-- Variables
*
      Double Precision  Func1, Eng
*
*-- Internal Variables
*
*   Integral Limit
      Double Precision  IntgLo, IntgUp
*   Integral Core Function
      Double Precision  Func1C
      External          Func1C
*   Integral Routine
      Double Precision  Itg5P2
      External          Itg5P2
*
*-----------------------------------------------------------------------

      IntgLo = ( (DSqrt(1.0D0+Eng)-1.0D0) * 2.0D0 / Eng ) ** 
     $         (1.0D0 / 6.0D0)
      IntgUp = ( (-DSqrt(1.0D0+Eng)-1.0D0) * 2.0D0 / Eng ) **
     $         (1.0D0 / 6.0D0)

      Func1 = Gmr * Itg5P2(Func1C, Eng, IntgLo, IntgUp, IGridN)

      Call Err
      Return
      End Function Func1

*-----------------------------------------------------------------------
*
      Function Func1C(Eng, Var)
*
*     Purpose: The integral core of Func1
*
      Implicit None
*
*-- Variables
*
      Double Precision  Func1C, Eng, Var
      Double Precision  Tmp, FuncP
      External          FuncP
*
*-----------------------------------------------------------------------

      Tmp = Eng - FuncP(Var)
      Call Test0(Tmp)
      Func1C = DSqrt(Tmp)

      End Function Func1C

*-----------------------------------------------------------------------
*
      Function FuncP(Var)
*
*     Purpose: Potential Function
*
      Implicit None
*
*-- Variables
*
      Double Precision  FuncP, Var
*
*-----------------------------------------------------------------------

      FuncP = 4.0D0 * ( Var**(-12.0D0) - Var**(-6.0D0) )  

      End Function FuncP

************************************************************************
*                                                                      *
*     MATHEMATICAL BACKGROUND UTILITIES                                *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*
      Subroutine Err
*
*     Purpose: The integral core of Func1
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
      End If

      End Subroutine

************************************************************************
*                                                                      *
*     MATHEMATICAL BACKGROUND UTILITIES                                *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*
      Function RtSec(Func, Num, RtLo, RtUp)
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
      Double Precision  RtSec, Func, Num, RtLo, RtUp
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
        RtSec = ( RtLoI * FncUp - RtUpI * FncLo ) / ( FncUp - FncLo )
        FncSec = Func(RtSec) - Num
        RtLoI = RtUpI
        FncLo = FncUp
        RtUpI = RtSec
        FncUp = FncSec
*   Converge Judgement
      If ( (Itr.LT.ItrMax).And.( DAbs(RtUpI-RtLoI).GT.DevX.Or.
     $      (DAbs(FncUp).GT.DevY.And.(FncUp-FncLo).GT.DevY))) then
*     Not converged
        Go To 10
      Else If ( (Itr.GE.ItrMax).And.( DAbs(RtUpI-RtLoI).GT.DevX.Or.
     $      (DAbs(FncUp).GT.DevY.And.(FncUp-FncLo).GT.DevY))) then
*     Not converged but reached the maximun iteration
        Print *, 'Error Info in Function [RtSec]'
        Print *, 'Not converged but reached the maximun iteration'
        Print *, 'Last Iteration Information:'
        Print *, 'Upper Guess: ', RtUpI, '  Lower Guess: ', RtLoI, 
     $           'Upper Func: ', FncUp, '  Lower Guess: ', FncLo
        IStat = 1
      Else
        Return
      End If

      Call Err
      Return
      End Function RtSec

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
      H = (IntgUp - IntgLo) / (4.0D0 * DBLE(Itr))

*   Initial and final point calculation
      Itg5P2 = 7.0D0 * (Func(VM, IntgLo) - Func(VM, IntgUp))

*   Internal grid point calculation
      Do 10 I = 1, Itr
        Itg5P2 = Itg5P2 + 32.0D0 * Func(VM, IntgLo + DBLE(I*4-3) * H)
     $                  + 12.0D0 * Func(VM, IntgLo + DBLE(I*4-2) * H)
     $                  + 32.0D0 * Func(VM, IntgLo + DBLE(I*4-1) * H)
     $                  + 14.0D0 * Func(VM, IntgLo + DBLE(I*4  ) * H)
  10  Continue

*   Coefficient multiplication
      Itg5P2 = Itg5P2 * H * 2.0D0 / 45.0D0

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

