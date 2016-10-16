CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM EXMPL4
C     Example 4: Born and Eikonal Approximations to Quantum Scattering
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company, Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop/ execute once for each set of param
        CALL PARAM        !get input from screen
        CALL ARCHON       !calculate cross sections
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON
C calculates the Born, eikonal, and optical theorem total cross sections
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E4'
C Local variables:
      REAL SINHLF                !sin theta_cut/2
      REAL THTCUT,COSCUT         !theta_cut/2, cos of the same
      REAL SIGBT,SIGBTF,SIGBTB   !Born cross sections
      REAL SIGET,SIGETF,SIGETB   !eikonal cross sections
      REAL SIGOPT,IMFE0          !optical cross sect, imag part of f(0)
      REAL COSTH                 !integration variable
      INTEGER ILEG,LANG          !index quad points, angle
      INTEGER NLINES             !number of lines printed to terminal
      INTEGER SCREEN             !send to terminal
      INTEGER PAPER              !make a hardcopy
      INTEGER FILE               !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     output summary of parameters                 
      IF (TTERM) CALL PRMOUT(OUNIT,NLINES)
      IF (TFILE) CALL PRMOUT(TUNIT,NLINES)
      IF (GFILE) CALL PRMOUT(GUNIT,NLINES)
C
C     find angle at which q*rcut=2*pi to divide integral into
C       forward and backward scattering
      SINHLF=PI/(K*RCUT)
      IF (SINHLF.GT.SQHALF) SINHLF=SQHALF
      THTCUT=2*ATAN(SINHLF/SQRT(1-SINHLF*SINHLF))
      COSCUT=COS(THTCUT)                
C
      SIGBT=0.0                         !zero sums
      SIGET=0.0
      IMFE0=0.0
C 
C     integrate from theta=0 to theta=theta_cut 
C      to find forward angle cross section
      DO 100 ILEG=1,NLEG
         COSTH=COSCUT+(XLEG(ILEG)+1)*(1-COSCUT)/2.0
         LANG=ILEG
         CALL DIFFCS(COSTH,ILEG,LANG,SIGET,SIGBT,IMFE0,NLINES)
100   CONTINUE
      SIGBTF=PI*(1-COSCUT)*SIGBT
      SIGETF=PI*(1-COSCUT)*SIGET
C
      SIGBT=0                           !zero sums
      SIGET=0
C 
C     integrate from theta=theta_cut to theta=pi
C      to find backward angle cross section
      DO 200 ILEG=1,NLEG
         COSTH=-1+(XLEG(ILEG)+1)*(COSCUT+1)/2.0
         LANG=ILEG+NLEG
         CALL DIFFCS(COSTH,ILEG,LANG,SIGET,SIGBT,IMFE0,NLINES)
200   CONTINUE
      SIGBTB=PI*(COSCUT+1)*SIGBT
      SIGETB=PI*(COSCUT+1)*SIGET
C 
C     add backward and forward scattering cross sections
      SIGBT=SIGBTF+SIGBTB
      SIGET=SIGETF+SIGETB
      SIGOPT=4*PI/K*IMFE0
C
C     output results
      IF (TTERM) CALL TOTOUT(OUNIT,THTCUT,SIGBT,SIGET,SIGOPT)
      IF (TFILE) CALL TOTOUT(TUNIT,THTCUT,SIGBT,SIGET,SIGOPT)
C
      IF (TTERM) CALL PAUSE('to continue...',1)
      IF (TTERM) CALL CLEAR
C
C     graphics output
      IF (GTERM) CALL GRFOUT(SCREEN)
      IF (GFILE) CALL GRFOUT(FILE)
      IF (GHRDCP) CALL GRFOUT(PAPER)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DIFFCS(COSTH,ILEG,LANG,SIGET,SIGBT,IMFE0,NLINES)
C calculate the differential cross section at a give angle (COSTH)
C and its contribution to the total cross section
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global Variables: 
      INCLUDE 'PARAM.E4'
      INCLUDE 'IO.ALL'
C Passed Variables
      REAL SIGBT              !Born cross sections (output)
      REAL SIGET              !eikonal cross sections (output)
      REAL IMFE0              !imag part of f(0) (output)
      REAL COSTH              !integration variable (input)
      INTEGER ILEG,LANG       !index quad points, angle (input)
      INTEGER NLINES          !number of lines printed to terminal (I/O)
C Local Variables:
      REAL THETA,Q            !angle, momentum transfer
      REAL DEGREE(96)         !angle
      REAL SIGE(96)           !eikonal differential cross section
      REAL SIGB(96)           !born differential cross section
      COMMON/RESULT/DEGREE,SIGE,SIGB  !pass to graphics routine only
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      THETA=ATAN(SQRT(1-COSTH*COSTH)/COSTH)    !find the angle
      IF (THETA.LT.0.0) THETA=THETA+PI
      DEGREE(LANG)=THETA*180.0/PI
C                                                                     
      Q=2*K*SIN(THETA/2.0)            !momentum transfer at this angle
C
      CALL BORN(Q,LANG,SIGB)          !calculate differential cross sect
      SIGBT=SIGBT+SIGB(LANG)*WLEG(ILEG) !and add to total cross sect
      CALL EIKONL(Q,LANG,SIGE,IMFE0)
      SIGET=SIGET+SIGE(LANG)*WLEG(ILEG)
C 
      IF (TTERM)                        !text output
     &  CALL TXTOUT(OUNIT,DEGREE(LANG),Q,SIGB(LANG),SIGE(LANG),NLINES)
      IF (TFILE) 
     &  CALL TXTOUT(TUNIT,DEGREE(LANG),Q,SIGB(LANG),SIGE(LANG),NLINES)
C      
      RETURN
      END                                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BORN(Q,LANG,SIGB)          
C calculates the differential cross section in the Born approximation
C for one value of the momentum transfer (Q)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global Variables:
      INCLUDE 'PARAM.E4'
C Passed Variables
      REAL    Q                !momentum transfer (input)
      INTEGER LANG             !index of the angle (input)
      REAL SIGB(96)            !Born differential cross section (output)
C Local Variables:
      REAL    R                !radial variable
      REAL    FBORN            !scattering amplitude
      INTEGER J                !index of quadrature points
C Function:
      REAL V                   !potential 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FBORN=0                                !zero sum
      DO 100 J=1,NLEG
         R=RMAX*(XLEG(J)+1)/2.0              !scale the variable
         FBORN=FBORN+SIN(Q*R)*V(R)*R*WLEG(J) !Gaussian quadrature
100   CONTINUE
      FBORN=-RMAX/(Q*HBARM)*FBORN
      SIGB(LANG)=FBORN*FBORN
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE EIKONL(Q,LANG,SIGE,IMFE0)
C calculates differential cross secton in the eikonal approximation
C for one value of the momentum transfer Q
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global Variables:
      INCLUDE 'PARAM.E4'
C Input/Output Variables
      REAL Q                       !momentum transfer (input)
      INTEGER LANG                 !index of the angle (input)
      REAL SIGE(96)                !eikonal diff cross section (output)
      REAL IMFE0                   !imag part of f(0) (output)
C Local Variables:
      REAL CHI(96)                 !profile function
      REAL B                       !impact parameter
      COMPLEX FE                   !scattering amplitude
      REAL ZMAX,ZZ                 !integration variable for CHI
      REAL R                       !radius
      REAL J0                      !Bessel function
      INTEGER M,J                  !index for quadrature
      REAL X                       !argument of Bessel function
      COMPLEX SQRTM1               !square root of minus 1 
C Functions:
      REAL BESSJ0                  !zeroth Bessel function eval at X
      REAL V                       !potential
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      SQRTM1=(0.,1.)               !constant
      FE=(0.,0.)                   !zero the scatt ampl
C
      DO 100 J=1,NLEG              !integrate over B
         B=RMAX*(XLEG(J)+1)/2.0
C
         IF (LANG .EQ. 1) THEN         !calculate profile function once
             CHI(J)=0.0
             ZMAX=SQRT(RMAX*RMAX-B*B)  !integration variable z
             DO 200 M=1,NLEG
                ZZ=ZMAX*(XLEG(M)+1)/2. !scale variable
                R=SQRT(ZZ*ZZ+B*B)      !calculate radius
                CHI(J)=CHI(J)+V(R)*WLEG(M)           !Gaussian quad
200          CONTINUE
             CHI(J)=-(ZMAX/(2*K*HBARM))*CHI(J)
             IMFE0=IMFE0+B*(SIN(CHI(J)))**2*WLEG(J)  !imag f(0)
         ENDIF
C
         X=Q*B
         J0=BESSJ0(X)
         FE=FE+B*J0*WLEG(J)*(EXP(2*SQRTM1*CHI(J))-1.)  !Gaussian quad
C
100   CONTINUE
C
      FE=-SQRTM1*K*RMAX*FE/2. !factor of RMAX/2 from change in variables
      SIGE(LANG)=ABS(FE)**2
      IF (LANG .EQ. 1) IMFE0=K*RMAX*IMFE0
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION V(R)
C calculates the potential at fixed R
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E4'
C Input/Output Variables:
      REAL R                        !radius (input)
C Local Variables
      REAL RR                       !radius for Lenz-Jensen (.ge. R1S)
      REAL U                        !temp variable for Lenz-Jensen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (POT .EQ. SQUARE) THEN     !square well
         V=-VZERO
C
      ELSE IF (POT .EQ. GAUSS) THEN !Gaussian well
         V=-VZERO*EXP(-2*R*R)
C
      ELSE                       !Lenz-Jensen
         RR=R                    !cutoff radius at location of 1S shell
         IF (R .LT. R1S) RR=R1S
         U=4.5397*Z6*SQRT(RR)                             
         V=-((Z*E2)/RR)*EXP(-U)*(1+U+U*U*(.3344+U*(.0485+.002647*U)))
      ENDIF
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION BESSJ0(X)
C calculates the zeroth order Bessel function at X
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global Variables:
      INCLUDE 'PARAM.E4'
C Input Variables:
      REAL X                      !argument of J_0
C Local Variables
      REAL Y,Y2,TEMP,TEMP2        !temporary variables to save values
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Y=X/3.
      Y2=Y*Y
C
      IF ( ABS(X) .LE. 3.) THEN
          TEMP=-.0039444+.00021*Y2
          TEMP=.0444479+Y2*TEMP
          TEMP=-.3163866+Y2*TEMP
C 
          TEMP=1.2656208+Y2*TEMP
          TEMP=-2.2499997+Y2*TEMP
          BESSJ0=1.+Y2*TEMP
      ELSE
          Y=1./Y
          TEMP=-7.2805E-04+1.4476E-04*Y
          TEMP=1.37237E-03+TEMP*Y
          TEMP=-9.512E-05+TEMP*Y
          TEMP=-.0055274+TEMP*Y
          TEMP=-7.7E-07+TEMP*Y
          TEMP=.79788456+TEMP*Y
          TEMP2=-2.9333E-04+1.3558E-04*Y
          TEMP2=-5.4125E-04+TEMP2*Y
          TEMP2=2.62573E-03+TEMP2*Y
          TEMP2=-3.954E-05+TEMP2*Y
          TEMP2=-4.166397E-02+TEMP2*Y
          TEMP2=X-.78539816+TEMP2*Y
          BESSJ0=TEMP*COS(TEMP2)/SQRT(X)
      ENDIF
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE QUAD(NLBAD)   
C subroutine to establish Gauss-Legendre abscissae and weights
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E4'
C Passed Variables
      LOGICAL NLBAD                   !allowed value of nleg?(output)
C Local Variables:
      INTEGER ILEG                    !index for abscissae and weights
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      NLBAD=.FALSE.
C 
      IF (NLEG.EQ.2) THEN
         XLEG(1)=.577350269189626
         WLEG(1)=1.                   !nleg=2
      ELSE IF (NLEG.EQ.3) THEN
         XLEG(1)=.774596669241483
         WLEG(1)=.555555555555556     !NLEG=3
         XLEG(2)=0.                  
         WLEG(2)=.888888888888889
      ELSE IF (NLEG.EQ.4) THEN
         XLEG(1)=.861136311594053  
         WLEG(1)=.347854845137454     !NLEG=4
         XLEG(2)=.339981043584856  
         WLEG(2)=.652145154862546
      ELSE IF (NLEG.EQ.5) THEN
         XLEG(1)=.906179845938664  
         WLEG(1)=.236926885056189     !NLEG=5
         XLEG(2)=.538469310105683 
         WLEG(2)=.478628670499366
         XLEG(3)=0.                  
         WLEG(3)=.568888888888889
      ELSE IF (NLEG.EQ.6) THEN
         XLEG(1)=.932469514203152  
         WLEG(1)=.17132449237917      !NLEG=6
         XLEG(2)=.661209386466266  
         WLEG(2)=.360761573048139
         XLEG(3)=.238619186083197  
         WLEG(3)=.467913934572691
      ELSE IF (NLEG.EQ.8) THEN
         XLEG(1)=.960289856497536  
         WLEG(1)=.101228536290376    !NLEG=8
         XLEG(2)=.796666477413627  
         WLEG(2)=.222381034453374
         XLEG(3)=.525532409916329  
         WLEG(3)=.313706645877887
         XLEG(4)=.18343464249565   
         WLEG(4)=.362683783378362
      ELSE IF (NLEG.EQ.10) THEN
         XLEG(1)=.973906528517172  
         WLEG(1)=.066671344308688    !NLEG=10
         XLEG(2)=.865063366688985  
         WLEG(2)=.149451349150581
         XLEG(3)=.679409568299024  
         WLEG(3)=.219086362515982
         XLEG(4)=.433395394129247  
         WLEG(4)=.269266719309996
         XLEG(5)=.148874338981631  
         WLEG(5)=.295524224714753
      ELSE IF (NLEG.EQ.12) THEN
         XLEG(1)=.981560634246719  
         WLEG(1)=.047175336386512    !NLEG=12
         XLEG(2)=.904117256370475  
         WLEG(2)=.106939325995318
         XLEG(3)=.769902674194305  
         WLEG(3)=.160078328543346
         XLEG(4)=.587317954286617  
         WLEG(4)=.203167426723066
         XLEG(5)=.36783149899818   
         WLEG(5)=.233492536538355
         XLEG(6)=.125233408511469  
         WLEG(6)=.249147045813403
      ELSE IF (NLEG.EQ.16) THEN
         XLEG(1)=.98940093499165   
         WLEG(1)=.027152459411754    !NLEG=16
         XLEG(2)=.944575023073233  
         WLEG(2)=.062253523938648
         XLEG(3)=.865631202387832  
         WLEG(3)=.095158511682493
         XLEG(4)=.755404408355003  
         WLEG(4)=.124628971255534
         XLEG(5)=.617876244402644  
         WLEG(5)=.149595988816577
         XLEG(6)=.458016777657227  
         WLEG(6)=.169156519395003
         XLEG(7)=.281603550779259  
         WLEG(7)=.182603415044924
         XLEG(8)=.095012509837637  
         WLEG(8)=.189450610455069
      ELSE IF (NLEG.EQ.20) THEN
         XLEG(1)=.993128599185094  
         WLEG(1)=.017614007139152    !NLEG=20
         XLEG(2)=.963971927277913  
         WLEG(2)=.040601429800386
         XLEG(3)=.912234428251325  
         WLEG(3)=.062672048334109
         XLEG(4)=.839116971822218  
         WLEG(4)=.083276741576704
         XLEG(5)=.74633190646015   
         WLEG(5)=.10193011981724
         XLEG(6)=.636053680726515  
         WLEG(6)=.118194531961518
         XLEG(7)=.510867001950827  
         WLEG(7)=.131688638449176
         XLEG(8)=.373706088715419  
         WLEG(8)=.142096109318382
         XLEG(9)=.227785851141645  
         WLEG(9)=.149172986472603
         XLEG(10)=.076526521133497 
         WLEG(10)=.152753387130725
      ELSE IF (NLEG.EQ.24) THEN
         XLEG(1)=.995187219997021  
         WLEG(1)=.012341229799987    !NLEG=24
         XLEG(2)=.974728555971309  
         WLEG(2)=.028531388628933
         XLEG(3)=.938274552002732  
         WLEG(3)=.044277438817419
         XLEG(4)=.886415527004401  
         WLEG(4)=.059298584915436
         XLEG(5)=.820001985973902  
         WLEG(5)=.07334648141108
         XLEG(6)=.740124191578554  
         WLEG(6)=.086190161531953
         XLEG(7)=.648093651936975  
         WLEG(7)=.097618652104113
         XLEG(8)=.545421471388839  
         WLEG(8)=.107444270115965
         XLEG(9)=.433793507626045  
         WLEG(9)=.115505668053725
         XLEG(10)=.315042679696163 
         WLEG(10)=.121670472927803
         XLEG(11)=.191118867473616 
         WLEG(11)=.125837456346828
         XLEG(12)=.064056892862605 
         WLEG(12)=.127938195346752
      ELSE IF (NLEG.EQ.32) THEN
         XLEG(1)=.997263861849481  
         WLEG(1)=.00701861000947     !NLEG=32
         XLEG(2)=.985611511545268  
         WLEG(2)=.016274394730905
         XLEG(3)=.964762255587506  
         WLEG(3)=.025392065309262
         XLEG(4)=.934906075937739  
         WLEG(4)=.034273862913021
         XLEG(5)=.896321155766052  
         WLEG(5)=.042835898022226
         XLEG(6)=.849367613732569  
         WLEG(6)=.050998059262376
         XLEG(7)=.794483795967942  
         WLEG(7)=.058684093478535
         XLEG(8)=.732182118740289  
         WLEG(8)=.065822222776361
         XLEG(9)=.663044266930215  
         WLEG(9)=.072345794108848
         XLEG(10)=.587715757240762 
         WLEG(10)=.07819389578707
         XLEG(11)=.506899908932229 
         WLEG(11)=.083311924226946
         XLEG(12)=.421351276130635 
         WLEG(12)=.087652093004403
         XLEG(13)=.331868602282127 
         WLEG(13)=.091173878695763
         XLEG(14)=.239287362252137 
         WLEG(14)=.093844399080804
         XLEG(15)=.144471961582796 
         WLEG(15)=.095638720079274
         XLEG(16)=.048307665687738 
         WLEG(16)=.096540088514727
      ELSE IF (NLEG.EQ.40) THEN
         XLEG(1)=.998237709710559  
         WLEG(1)=.004521277098533    !NLEG=40
         XLEG(2)=.990726238699457  
         WLEG(2)=.010498284531152
         XLEG(3)=.977259949983774  
         WLEG(3)=.016421058381907
         XLEG(4)=.957916819213791  
         WLEG(4)=.022245849194166
         XLEG(5)=.932812808278676  
         WLEG(5)=.027937006980023
         XLEG(6)=.902098806968874  
         WLEG(6)=.033460195282547
         XLEG(7)=.865959503212259  
         WLEG(7)=.038782167974472
         XLEG(8)=.824612230833311  
         WLEG(8)=.043870908185673
         XLEG(9)=.778305651426519  
         WLEG(9)=.048695807635072
         XLEG(10)=.727318255189927 
         WLEG(10)=.053227846983936
         XLEG(11)=.671956684614179 
         WLEG(11)=.057439769099391
         XLEG(12)=.61255388966798  
         WLEG(12)=.061306242492928
         XLEG(13)=.549467125095128 
         WLEG(13)=.064804013456601
         XLEG(14)=.483075801686178 
         WLEG(14)=.067912045815233
         XLEG(15)=.413779204371605 
         WLEG(15)=.070611647391286
         XLEG(16)=.341994090825758 
         WLEG(16)=.072886582395804
         XLEG(17)=.268152185007253 
         WLEG(17)=.074723169057968
         XLEG(18)=.192697580701371 
         WLEG(18)=.076110361900626
         XLEG(19)=.116084070675255 
         WLEG(19)=.077039818164247
         XLEG(20)=.03877241750605  
         WLEG(20)=.077505947978424
      ELSE IF (NLEG.EQ.48) THEN
         XLEG(1)=.998771007252426  
         WLEG(1)=.003153346052305     !NLEG=48
         XLEG(2)=.99353017226635   
         WLEG(2)=.007327553901276
         XLEG(3)=.984124583722826  
         WLEG(3)=.011477234579234
         XLEG(4)=.970591592546247  
         WLEG(4)=.015579315722943
         XLEG(5)=.95298770316043   
         WLEG(5)=.019616160457355
         XLEG(6)=.931386690706554  
         WLEG(6)=.023570760839324
         XLEG(7)=.905879136715569  
         WLEG(7)=.027426509708356
         XLEG(8)=.876572020274247  
         WLEG(8)=.031167227832798
         XLEG(9)=.843588261624393  
         WLEG(9)=.03477722256477
         XLEG(10)=.807066204029442 
         WLEG(10)=.03824135106583
         XLEG(11)=.76715903251574  
         WLEG(11)=.041545082943464
         XLEG(12)=.724034130923814 
         WLEG(12)=.044674560856694
         XLEG(13)=.677872379632663 
         WLEG(13)=.04761665849249
         XLEG(14)=.628867396776513 
         WLEG(14)=.050359035553854
         XLEG(15)=.577224726083972 
         WLEG(15)=.052890189485193
         XLEG(16)=.523160974722233 
         WLEG(16)=.055199503699984
         XLEG(17)=.466902904750958 
         WLEG(17)=.057277292100403
         XLEG(18)=.408686481990716 
         WLEG(18)=.059114839698395
         XLEG(19)=.34875588629216  
         WLEG(19)=.060704439165893
         XLEG(20)=.287362487355455 
         WLEG(20)=.062039423159892
         XLEG(21)=.224763790394689 
         WLEG(21)=.063114192286254
         XLEG(22)=.161222356068891 
         WLEG(22)=.063924238584648
         XLEG(23)=.097004699209462 
         WLEG(23)=.06446616443595
         XLEG(24)=.032380170962869 
         WLEG(24)=.064737696812683
      ELSE
         NLBAD=.TRUE.     !if NLEG is a disallowed value
      ENDIF
C 
      IF (.NOT. NLBAD) THEN
         DO 100 ILEG=1,NLEG/2            !since the weights and abscissa
            XLEG(NLEG-ILEG+1)=-XLEG(ILEG)!are even and odd functions
            WLEG(NLEG-ILEG+1)=WLEG(ILEG) !functions respectively
100      CONTINUE
      ENDIF
C 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INIT
C initializes constants, displays header screen,
C initializes arrays for input parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'MENU.ALL'
      INCLUDE 'PARAM.E4'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL SETUP                !get environment parameters
C                                                                     
C     display header screen     
      DESCRP(1)= 'EXAMPLE 4'
      DESCRP(2)= 'Born and Eikonal Approximations'
      DESCRP(3)= 'to Quantum Scattering'
      NHEAD=3
C 
C     text output description
      DESCRP(4)= 'angle, momentum transfer, Born and Eikonal '
     +   //'differential cross sections,'
      DESCRP(5)= 'and Born, Eikonal, and optical theorem '
     +   //'total cross section'
      NTEXT=2                               
C 
C     graphics output description
      DESCRP(6)= 'differential cross section vs. theta'
      NGRAPH=1
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C 
C     set constant values
      PI=3.1415926
      SQHALF=SQRT(0.5)
      E2=14.409                      !electron charge, squared
      RMAX=2.0
      HBARM=7.6359                   !hbar**2/mass electron
C 
      CALL MENU                      !setup constant part of menu
C                
      MTYPE(13)=MTITLE
      MPRMPT(13)='Potential Function Options:'
      MLOLIM(13)=2
      MHILIM(13)=1
C 
      MTYPE(14)=MTITLE
      MPRMPT(14)='1) Lenz-Jensen Potential: electron & neutral atom'
      MLOLIM(14)=0
      MHILIM(14)=0
C 
      MTYPE(15)=MTITLE
      MPRMPT(15)='2) Square Well'
      MLOLIM(15)=0            
      MHILIM(15)=0
C 
      MTYPE(16)=MTITLE
      MPRMPT(16)='3) Gaussian Well'
      MLOLIM(16)=0
      MHILIM(16)=1
C 
      MTYPE(17)=MCHOIC
      MPRMPT(17)='Enter Choice'
      MTAG(17)='18 20 20'
      MLOLIM(17)=1
      MHILIM(17)=3
      MINTS(17)=1
      MREALS(17)=1.
C 
      MTYPE(18)=NUM
      MPRMPT(18)='Enter charge of the atomic nucleus'
      MTAG(18)='Z'
      MLOLIM(18)=1
      MHILIM(18)=108
      MINTS(18)=4
C
      MTYPE(19)=SKIP
      MREALS(19)=21
C 
      MTYPE(20)=FLOAT
      MPRMPT(20)='Enter Vzero (eV)'
      MTAG(20)='Vzero (eV)'
      MLOLIM(20)=0.0
      MHILIM(20)=5000.
      MREALS(20)=20.0
C 
      MTYPE(21)=FLOAT
      MPRMPT(21)='Enter Energy (eV)'
      MTAG(21)='Energy (eV)'
      MLOLIM(21)=0.0
      MHILIM(21)=5.0E+06
      MREALS(21)=20.0
C 
      MTYPE(22)=SKIP
      MREALS(22)=35
C 
      MTYPE(38)=NUM
      MPRMPT(38)='Number of quadrature points'
      MTAG(38)='number of quadrature points'
      MLOLIM(38)=2
      MHILIM(38)=48
      MINTS(38)=20
C 
      MTYPE(39)=SKIP
      MREALS(39)=60
C 
      MSTRNG(MINTS(75))= 'exmpl4.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C 
      MSTRNG(MINTS(86))= 'exmpl4.grf'
C                     
      MTYPE(87)=SKIP
      MREALS(87)=90.
C                
      RETURN
      END                            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PARAM
C gets parameters from screen
C ends program on request
C closes old files
C maps menu variables to program variables
C opens new files
C calculates all derivative parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E4'
C Local Variables:
      LOGICAL NLBAD                  !allowed value of nleg?
C map between menu items and parameters
      INTEGER IPOT,IZ,IVZERO,IE,INLEG
      PARAMETER (IPOT   = 17 )
      PARAMETER (IZ     = 18 )
      PARAMETER (IVZERO = 20 )
      PARAMETER (IE     = 21 )
      PARAMETER (INLEG  = 38 )
C Functions:
      INTEGER GETINT                 !gets integer value from screen
      LOGICAL LOGCVT                 !converts 1 and 0 to true and false 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get input from terminal
      CALL CLEAR
      CALL ASK(1,ISTOP)
C 
C     stop program if requested
      IF (MREALS(IMAIN) .EQ. STOP) CALL FINISH 
C
C     close files if necessary
      IF (TNAME .NE. MSTRNG(MINTS(ITNAME))) 
     +     CALL FLCLOS(TNAME,TUNIT)
      IF (GNAME .NE. MSTRNG(MINTS(IGNAME))) 
     +     CALL FLCLOS(GNAME,GUNIT)
C
C     physical and numerical parameters
      POT=MINTS(IPOT)                                       
      Z=MINTS(IZ)
      VZERO=MREALS(IVZERO)
      E=MREALS(IE)
      NLEG=MINTS(INLEG)
C
C     derived constants
      K=SQRT(2*E/HBARM)
      IF (POT .EQ. LENZ) THEN
         RCUT=1.
         Z6=Z**.166667
         R1S=.529/Z
      ELSE
         RCUT=2.
      END IF
C
C     find weights and abscissae for Gauss-Legendre integration
      CALL QUAD(NLBAD)
100   IF (NLBAD) THEN
         WRITE (OUNIT,*) ' Allowed values for quadrature points are ',
     +                   '2 3 4 5 6 8 10 12 16 20 24 32 40 48'
         NLEG=GETINT(20,2,48,'Enter new number of quadrature points')
         MINTS(INLEG)=NLEG
         CALL QUAD(NLBAD)
      GOTO 100
      END IF
      CALL CLEAR
C 
C     text output
      TTERM=LOGCVT(MINTS(ITTERM))
      TFILE=LOGCVT(MINTS(ITFILE))
      TNAME=MSTRNG(MINTS(ITNAME))
C
C     graphics output
      GTERM=LOGCVT(MINTS(IGTERM))
      GHRDCP=LOGCVT(MINTS(IGHRD))
      GFILE=LOGCVT(MINTS(IGFILE))                              
      GNAME=MSTRNG(MINTS(IGNAME))
C 
C     open files 
      IF (TFILE) CALL FLOPEN(TNAME,TUNIT)
      IF (GFILE) CALL FLOPEN(GNAME,GUNIT)
      !files may have been renamed
      MSTRNG(MINTS(ITNAME))=TNAME
      MSTRNG(MINTS(IGNAME))=GNAME
C 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRMOUT(MUNIT,NLINES)
C outputs parameters to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'IO.ALL'
       INCLUDE 'PARAM.E4'
C Passed variables:               
       INTEGER MUNIT            !unit number for output (input)
       INTEGER NLINES           !number of lines written so far (I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (MUNIT .EQ. OUNIT) CALL CLEAR
C 
       WRITE (MUNIT,2)
       WRITE (MUNIT,4)
       WRITE (MUNIT,6) E
       WRITE (MUNIT,8) NLEG
C
       IF (POT .EQ. LENZ) THEN
           WRITE (MUNIT,10) Z
       ELSE IF (POT .EQ. SQUARE) THEN
           WRITE (MUNIT,12) VZERO
       ELSE IF (POT .EQ. GAUSS) THEN
           WRITE (MUNIT,14) VZERO
       END IF
       WRITE (MUNIT,15)
C 
C      different header for text and graphics files
       IF (MUNIT .EQ. GUNIT) THEN
         WRITE (MUNIT,2)
       ELSE
         WRITE (MUNIT,2)
         WRITE (MUNIT,16)
         WRITE (MUNIT,17)
         WRITE (MUNIT,18)
       END IF
C 
       NLINES=10
C
2      FORMAT (' ')
4      FORMAT (' Output from example 4: Born and eikonal '
     +         'approximations to quantum scattering')
6      FORMAT (' Energy (eV) =', 1PE10.3)
8      FORMAT (' Number of quadrature points = ' , I2)
10     FORMAT (' Lenz Jensen potential with Z = ', I3)
12     FORMAT (' Square well potential with Vzero = ' 1PE10.3)
14     FORMAT (' Gaussian potential with Vzero = ' 1PE10.3)
15     FORMAT (' Theta is in degrees')
16     FORMAT (11X,'Theta',13X,'Q',14X,'Born',14X,'Eikonal')
17     FORMAT (10X,'degrees',11X,'A*-1',12X,'A**2',14X,'  A**2 ')
18     FORMAT (10X,'-------',11X,'----',12X,'----',14X,'-------')
C 
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,DEGREE,Q,SIGB,SIGE,NLINES)
C writes out differential cross section at each angle
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E4'
C Input variables:
      INTEGER MUNIT           !output unit specifier
      INTEGER NLINES          !number of lines printed to screen (I/O)
      REAL DEGREE,Q           !angle and momentum transfer
      REAL SIGB,SIGE          !Born and eikonal diff cross section
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     if screen is full, clear and type headings again
      IF ((MOD(NLINES,TRMLIN-5) .EQ. 0) 
     +                          .AND. (MUNIT .EQ. OUNIT)) THEN
         CALL PAUSE('to continue...',1)
         CALL CLEAR
         WRITE (MUNIT,16)
         WRITE (MUNIT,17)
         WRITE (MUNIT,18)
         WRITE (MUNIT,2)
         NLINES=NLINES+3
      END IF
C
      WRITE (MUNIT,20) DEGREE,Q,SIGB,SIGE
C
C     keep track of printed lines only for terminal output
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+1
C
20    FORMAT(9X,F7.2,9X,F7.2,9X,1PE10.3,9X,1PE10.3)
2     FORMAT (' ')
16    FORMAT (11X,'Theta',13X,'Q',14X,'Born',14X,'Eikonal')
17    FORMAT (10X,'degrees',11X,'A*-1',12X,'A**2',14X,'  A**2 ')
18    FORMAT (10X,'-------',11X,'----',12X,'----',14X,'-------')
C
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TOTOUT(MUNIT,THTCUT,SIGBT,SIGET,SIGOPT)
C writes out total cross section to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E4'
C Input variables:
      INTEGER MUNIT               !output unit specifier
      REAL THTCUT                 !theta cut
      REAL SIGBT,SIGET,SIGOPT     !total cross section with diff methods
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE (MUNIT,*) ' '
      WRITE (MUNIT,10) THTCUT*180/PI
      WRITE (MUNIT,20) SIGBT
      WRITE (MUNIT,30) SIGET
      WRITE (MUNIT,40) SIGOPT
C
10    FORMAT ('                   theta cut = ', F7.2,' degrees')
20    FORMAT ('            Born sigma total =',1PE10.3,' Angstroms**2')
30    FORMAT ('         Eikonal sigma total =',1PE10.3,' Angstroms**2')
40    FORMAT (' Optical theorem sigma total =',1PE10.3,' Angstroms**2')
C
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFOUT(DEVICE)
C graphs differential cross section vs. angle on a semi-log scale
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables                                                   
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E4'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE         !which device is being used?
C Local variables
      INTEGER ILEG           !indexes angle              
      INTEGER EXPMAX,EXPMIN  !min and max exp for diff cross section
      CHARACTER*9 CE         !E as a character string
      REAL DEGREE(96)        !angle
      REAL SIGE(96)          !eikonal differential cross section
      REAL SIGB(96)          !Born differential cross section
      INTEGER LEN            !string length
      INTEGER SCREEN         !send to terminal
      INTEGER PAPER          !make a hardcopy
      INTEGER FILE           !send to a file
      COMMON/RESULT/DEGREE,SIGE,SIGB
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     calculate parameters for graphing
      IF (DEVICE .NE. FILE) THEN
          NPLOT=1                        !how many plots?
          IPLOT=1
C 
          YMAX=0.                        !find limits on data points
          YMIN=SIGB(1)
          DO 20 ILEG=1,2*NLEG
             IF (SIGB(ILEG) .GT. YMAX) YMAX=SIGB(ILEG)
             IF (SIGE(ILEG) .GT. YMAX) YMAX=SIGE(ILEG)
             IF (SIGB(ILEG) .LT. YMIN) YMIN=SIGB(ILEG)
             IF (SIGE(ILEG) .LT. YMIN) YMIN=SIGE(ILEG)
20        CONTINUE
C         find integer limits on exponent
          EXPMAX=INT(LOG10(YMAX))
          IF (YMAX .GT. 1.) EXPMAX =EXPMAX+1
          EXPMIN=INT(LOG10(YMIN))
          IF (YMIN .LT. 1.) EXPMIN=EXPMIN-1
          YMAX=10.**EXPMAX
          YMIN=10.**EXPMIN
C 
          XMIN=DEGREE(1)
          XMAX=DEGREE(2*NLEG)
          Y0VAL=XMIN
          X0VAL=YMIN
C          
          NPOINT=2*NLEG
C 
          ILINE=1                      !line and symbol styles
          ISYM=4
          IFREQ=1
          NXTICK=5
          NYTICK=EXPMAX-EXPMIN
          IF (NYTICK .GT. 8) THEN      !keep number of ticks small
             IF (MOD(NYTICK,2) .EQ. 0) THEN
                NYTICK=NYTICK/2
             ELSE
                NYTICK=8
             END IF
          END IF                           
C 
          CALL CONVRT(E,CE,LEN)                      !titles and labels
          INFO = 'Born(X) and Eikonal(0) Approximations'
          LABEL(1)= 'Angle (Degrees)'
          IF (POT .EQ. LENZ) THEN
            TITLE = ' Lenz-Jensen Potential, Energy (eV)='//CE
            LABEL(2)= 'Differential Cross Section (Angstroms**2)'
          ELSE IF (POT .EQ. SQUARE) THEN                       
            TITLE = 'Square Well Potential, Energy (eV)='//CE
            LABEL(2)= 'Differential Cross Section (Angstroms**2)'
          ELSE IF (POT .EQ. GAUSS) THEN
            TITLE = 'Gaussian Potential, Energy (eV)='//CE
            LABEL(2)= 'Differential Cross Section (Angstroms**2)'
          END IF
C 
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
          CALL LGLNAX                          !draw axes
      END IF
C                                                      
C     output results
      IF (DEVICE .EQ. FILE) THEN
          WRITE (GUNIT,16)
          WRITE (GUNIT,18)
          WRITE (GUNIT,70) (DEGREE(ILEG),SIGB(ILEG),SIGE(ILEG),
     +        ILEG=1,2*NLEG)
      ELSE
         CALL XYPLOT (DEGREE,SIGB)              !plot born
         ISYM=1
         CALL XYPLOT (DEGREE,SIGE)              !plot eikonal
      END IF
C 
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)  !end graphics package
      IF (DEVICE .EQ. SCREEN) CALL TMODE        !switch to text mode
C
16    FORMAT (10X,'Theta',13X,'Born',14X,'Eikonal')
18    FORMAT (10X,'-----',13X,'----',14X,'-------')
70    FORMAT (9X,F7.2,9X,1PE10.3,9X,1PE10.3)
100   FORMAT (/,' Patience, please; output going to a file.')
C
      RETURN
      END
