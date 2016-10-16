CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM PROJ4
C     PROJECT 4: Partial-wave solution of quantum scattering
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT          !display header screen, setup parameters
5     CONTINUE           !main loop/ execute once for each set of param
        CALL PARAM       !get input from screen
        CALL ARCHON      !calculate total and differential cross section
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON
C calculates differential and total cross section using 
C partial-wave expansion; scattering is from a central potential 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P4'
      INCLUDE 'IO.ALL'
C Local variables:
      INTEGER  L                  !which partial wave
      REAL DELTA(0:MAXL)          !phase shift
      REAL DELDEG(0:MAXL)         !phase shift in degrees
      REAL SIGMA(0:MAXL),SIGTOT   !total cross section
      REAL DSIGMA(0:MAXANG)       !differential cross section
      COMPLEX F(0:MAXANG)         !scattering amplitude
      REAL CD,SD                  !sins and cosines
      INTEGER ITHETA              !index the angle
      INTEGER NLINES              !number of lines written to terminal
      INTEGER DEVICE              !current graphics device
      INTEGER SCREEN              !send to terminal
      INTEGER PAPER               !make a hardcopy
      INTEGER FILE                !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     output summary of parameters
      IF ((TTERM) .AND. (.NOT. GTERM)) CALL PRMOUT(OUNIT,NLINES)
      IF (TFILE) CALL PRMOUT(TUNIT,NLINES)
      IF (GFILE) CALL PRMOUT(GUNIT,NLINES)
C
C     find the phase shifts by integrating the time-indp Schrodinger eq.
      SIGTOT=0.
      DO 100 L=LSTART,LSTOP
        CALL NUMERV(L,DELTA)
        DELDEG(L)=180.0*DELTA(L)/PI
        SIGMA(L)=4*PI/K2*(2*L+1)*SIN(DELTA(L))**2
        SIGTOT=SIGTOT+SIGMA(L)
C
        IF (GTERM) THEN             !output results for this L value
            CALL GRFOUT(SCREEN,L,DELDEG(L),SIGMA(L))
        ELSE IF (TTERM) THEN
            CALL TXTOUT(OUNIT,NLINES,L,DELDEG(L),SIGMA(L))
        END IF
        IF (TFILE) CALL TXTOUT(TUNIT,NLINES,L,DELDEG(L),SIGMA(L))
100   CONTINUE                                                            
C
C     output data for all L 
      IF ((TTERM) .AND. (GTERM)) CALL DELOUT(DELDEG,SIGMA)
C
C     calculate differential cross section by summing over partial waves
      DO 300 ITHETA=0,NANG
         F(ITHETA)=(0.,0.)
         DO 200 L=LSTART,LSTOP
            SD=SIN(DELTA(L))
            CD=COS(DELTA(L))     
            F(ITHETA)=F(ITHETA)+(2*L+1)*SD*PL(L,ITHETA)*CMPLX(CD,SD)
200      CONTINUE
         F(ITHETA)=F(ITHETA)/K
         DSIGMA(ITHETA)=(ABS(F(ITHETA)))**2
300   CONTINUE           
C
      IF (TTERM) CALL PAUSE('to see total cross sections...',1)
      IF (TTERM) CALL CLEAR
C
      IF (TTERM) CALL SIGOUT(OUNIT,F,DSIGMA,SIGTOT)
      IF (TFILE) CALL SIGOUT(TUNIT,F,DSIGMA,SIGTOT)
C                         
      IF (TTERM) CALL PAUSE('to continue...',1)
      IF (TTERM) CALL CLEAR
C              
C     graphics output
      IF (GTERM) CALL GRFSIG(SCREEN,DSIGMA,SIGTOT)
      IF (GFILE) CALL GRFSIG(FILE,DSIGMA,SIGTOT)
      IF (GHRDCP) CALL GRFSIG(PAPER,DSIGMA,SIGTOT)
C                     
      RETURN                                               
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine NUMERV(L,DELTA)
C find the phase shifts DELTA for fixed L by integrating 
C the radial wave equation (for both the scattered and free wave)
C using the Numerov algorithm
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P4'
C Input/Output Variables:
      INTEGER  L               !which partial wave (input)
      REAL DELTA(0:MAXL)       !phase shift (output)
C Local Variables:                    
      REAL PSI(MAXN+MAXE),PSIF(MAXN+MAXE)   !wave functions
      REAL PSIMAX,PSFMAX       !maxima of wave functions
      REAL VEFF(MAXN+MAXE)     !effective pot for fixed L
      REAL CONST               !factor in Numerov algorithm
      REAL VFREE(MAXN+MAXE)    !centrifugal potential
      REAL KI,KIM1,KIP1        !Numerov coefficients for scatt wave
      REAL KIF,KIM1F,KIP1F     !Numerov coefficients for free wave
      REAL G,NUMER,DENOM       !factors to calculate phase shift
      REAL LL                  !useful constant for centrifugal pot
      INTEGER IR,MR            !index for radius 
      COMMON/RESULT/PSI,PSIF,PSIMAX,PSFMAX,VEFF  !pass results to GRFOUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CONST=DR*DR/HBARM/6         !useful constant
C             
C     calculate the centrifugal barrier and effective potential
      LL=L*(L+1)*HBARM/2.0
      DO 200 IR=1,NPTS+NXTRA
         VFREE(IR)=LL/(R(IR)*R(IR))     
         VEFF(IR)=V(IR)+VFREE(IR)
200   CONTINUE
C
C     beginning values for Numerov algorithm
C       we are finding both scattered and free wave functions
      PSI(1)=1E-25                   !scattered wave initial conditions
      KIM1=CONST*(E-VEFF(1))         !k of Numerov algorithm
      PSI(2)=(2.-12.*KIM1)*PSI(1)    !use information from the second
      KI=CONST*(E-VEFF(2))           ! derivative to get started
C
      PSIF(1)=1E-25                  !free wave initial conditions
      KIM1F=CONST*(E-VFREE(1))
      PSIF(2)=(2.-12.*KIM1F)*PSIF(1)
      KIF=CONST*(E-VFREE(2))
C                       
      DO 300 IR=2,NPTS+NXTRA-1
C
        KIP1=CONST*(E-VEFF(IR+1))        !scattered wave K
        PSI(IR+1)=((2-10*KI)*PSI(IR)-(1+KIM1)*PSI(IR-1))/(1+KIP1)
        IF (PSI(IR+1).GT.1E+15) THEN     !prevent overflows
            DO 400 MR=1,IR+1
400         PSI(MR)=PSI(MR)*0.00001
        ENDIF
        KIM1=KI                          !roll values
        KI=KIP1
C
        KIP1F=CONST*(E-VFREE(IR+1))      !free wave  K
        PSIF(IR+1)=((2-10*KIF)*PSIF(IR)-(1+KIM1F)*PSIF(IR-1))/(1+KIP1F)
        IF (PSIF(IR+1).GT.1E+15) THEN    !prevent overflows
           DO 500 MR=1,IR+1
500        PSIF(MR)=PSIF(MR)*0.00001
        ENDIF
        KIM1F=KIF                        !roll values
        KIF=KIP1F
C
300   CONTINUE
C 
C     calculate delta (see Eq. IV.8)
      G=RMAX*PSI(NPTS+NXTRA)/((RMAX+RXTRA)*PSI(NPTS))
      NUMER=G*JL(1,L)-JL(2,L)
      DENOM=G*NL(1,L)-NL(2,L)
      DELTA(L)=ATAN(NUMER/DENOM)
      CALL FIXDEL(DELTA(L),PSI,PSIF,PSIMAX,PSFMAX) !fix delta absolutely
C
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FIXDEL(DELTA,PSI,PSIF,PSIMAX,PSFMAX)  
C resolves the ambiguity in DELTA using info from the wave functions
C PSI and PSIF; finds maxima in wave functions PSIMAX, PSIFMAX
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P4'
C Passed Variables:
      REAL DELTA                          !phase shift (I/O)
      REAL PSI(MAXN+MAXE),PSIF(MAXN+MAXE) !wave functions (input)
      REAL PSIMAX,PSFMAX                  !max of wave functions(output)
C Local Variables:
      INTEGER BUMP,NODES,BUMPF,NODESF  !number of nodes and antinodes
      INTEGER I                        !radial index
      LOGICAL DPOS                     !is the phase shift positive?
      INTEGER N                        !how many factors of pi to add
C Function:
      REAL SGN                         !returns sign of argument
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     find number of maxima, nodes, and psimax for the scattered wave
      BUMP=0
      NODES=0
      PSIMAX=0.0
      DO 100 I=2,NPTS+NXTRA-1
         IF (SGN(PSI(I)) .NE. SGN(PSI(I-1))) NODES=NODES+1
         IF (ABS(PSI(I)).GT. PSIMAX) PSIMAX=ABS(PSI(I))                 
         IF ( SGN(PSI(I+1)-PSI(I)) .NE. SGN(PSI(I)-PSI(I-1)) )
     +        BUMP=BUMP+1
100   CONTINUE
C             
C     find number of maxima, nodes, and psimax for the free wave
      BUMPF=0
      NODESF=0
      PSFMAX=0.0
      DO 200 I=2,NPTS+NXTRA-1
         IF (SGN(PSIF(I)) .NE. SGN(PSIF(I-1))) NODESF=NODESF+1
         IF (ABS(PSIF(I)).GT.PSFMAX) PSFMAX=ABS(PSIF(I))
         IF ( SGN(PSIF(I+1)-PSIF(I)) .NE. SGN(PSIF(I)-PSIF(I-1)) )
     +              BUMPF=BUMPF+1
200   CONTINUE
C
      DPOS=DELTA .GT. 0.0            !is the phase shift positive?
C
C     make sign of phase shift agree with convention
C     viz., an attractive potential has positive delta
      N=(NODES+BUMP-NODESF-BUMPF)
      IF (ATTRCT .AND. .NOT. DPOS) THEN
         N=(N+1)/2
      ELSE IF (.NOT. ATTRCT .AND. DPOS) THEN
         N=(N-1)/2
      ELSE          
         N=N/2
      END IF
      DELTA=DELTA+N*PI
C
      RETURN 
      END                                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION SGN(X)
C returns the sign of x
      REAL X
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (X.EQ.0) THEN
         SGN=0.0
      ELSE IF (X.LT.0) THEN
         SGN=-1.0
      ELSE IF (X.GT.0) THEN
         SGN=1.0
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE POTNTL
C fill the array V for choice of potential
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P4'
C Local variables:
      INTEGER IR                    !index for radius
      REAL U                        !temp value for Lenz-Jensen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     save values of the radius
      DO 100 IR=1,NPTS+NXTRA
         R(IR)=IR*DR
100   CONTINUE
C
      IF (POT .EQ. LENZ) THEN       !Lenz-Jensen
        DO 200 IR=1,NPTS
          U=4.5397*Z6*SQRT(R(IR))
          V(IR)=-((Z*E2)/R(IR))*EXP(-U)*
     +        (1+U+U*U*(0.3344+U*(0.0485+0.002647*U)))
200     CONTINUE
        ATTRCT=.TRUE.
C
      ELSE IF (POT .EQ. SQUARE) THEN      !square well
         DO 300 IR=1,NPTS*3/4
           V(IR)=-VZERO
300      CONTINUE
         DO 400 IR=NPTS*3/4,NPTS
           V(IR)=0.0
400      CONTINUE
         IF (VZERO .LE. 0) ATTRCT=.FALSE.
         IF (VZERO .GT. 0) ATTRCT=.TRUE.
C
      ELSE IF (POT .EQ. GAUSS) THEN       !Gaussian well
         DO 500 IR=1,NPTS
            V(IR)=-VZERO*EXP(-R(IR)*R(IR))
500      CONTINUE
         IF (VZERO .LE. 0) ATTRCT=.FALSE.
         IF (VZERO .GT. 0) ATTRCT=.TRUE.
C
      ENDIF
C
C     for all potentials, V is zero outside RMAX
      DO 600 IR=NPTS+1,NPTS+NXTRA
         V(IR)=0.
600   CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LPINIT
C calculates the Legendre Polynomials P0 and P1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global Variables:
      INCLUDE 'PARAM.P4'
C Local Variables:
      INTEGER    ITHETA                 !indexes the angle
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 100 ITHETA=0,NANG
         THETA(ITHETA)=ITHETA*PI/NANG
         DEGREE(ITHETA)=ITHETA*180/NANG
         CTHETA(ITHETA)=COS(THETA(ITHETA))
         PL(0,ITHETA)=1
         PL(1,ITHETA)=CTHETA(ITHETA)
100   CONTINUE
      LTABLE=1
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LEGPOL
C calculates the Legendre Polynomials from LTABLE to LSTOP
C  using recursion relation
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE    'PARAM.P4'
C Local Variables:
      INTEGER   ITHETA,L               !index angle and L
      REAL X                           !cos(theta)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 100 ITHETA=0,NANG
       X=CTHETA(ITHETA)
       DO 200 L=LTABLE,LSTOP-1
         PL(L+1,ITHETA)=((2*L+1)*X*PL(L,ITHETA)-L*PL(L-1,ITHETA))/(L+1)
200    CONTINUE
100   CONTINUE
      LTABLE=LSTOP
C
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SBSSEL
C calculates the Spherical Bessel functions at RMAX and RXTRA+RMAX
C  for L values from 0 to LSTOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE    'PARAM.P4'
C Local Variables:
      REAL       J,JM1,JP1            !temp values for JL backwd recursn
      REAL       NORM                 !normalizing factor for JL's
      INTEGER    L,LUPPER             !index of JL and NL
      INTEGER    IR                   !indexes RMAX or RMAX+RXTRA
      REAL       X                    !function argument
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 1000 IR=1,2                  !needed at RMAX and RMAX+RXTRA
C       X values for which we are calculating spherical Bessel functions
        IF (IR .EQ. 1) X=K*RMAX
        IF (IR .EQ. 2) X=K*(RMAX+RXTRA)
C         
C       obtain NL's by stable forward recursion
        NL(IR,0)=-COS(X)/X
        NL(IR,1)=-COS(X)/(X*X)-SIN(X)/X
        DO 100 L=1,LSTOP-1
           NL(IR,L+1)=(2*L+1)*NL(IR,L)/X-NL(IR,L-1)
100     CONTINUE
C 
C       obtain JL's by stable backward recursion
        LUPPER=X+10.             !start at L such that JL is neglgbl
        JP1=0
        J=9.9999999E-21          !arbitrary beginning value
        DO 200 L=LUPPER,1,-1
           JM1=(2*L+1)*J/X-JP1
           IF ((L-1).LE. LSTOP) JL(IR,L-1)=JM1   !save needed L values
           JP1=J                                 !roll values
           J=JM1
200     CONTINUE
        NORM=SIN(X)/X/JL(IR,0)   !normalize using analytic form of JO
        DO 300 L=0,LSTOP
          JL(IR,L)=NORM*JL(IR,L)
300     CONTINUE
1000  CONTINUE
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
      INCLUDE 'PARAM.P4'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL SETUP                !get environment parameters
C 
C     display header screen     
      DESCRP(1)= 'PROJECT 4'
      DESCRP(2)= 'Partial-wave solution to quantum scattering'
      NHEAD=2
C
C     text output description
      DESCRP(3)= 'phase shift and total cross section for each L;'
      DESCRP(4)= 'differential cross section at several angles'
      NTEXT=2
C
C     graphics output description
      DESCRP(5)='Veffective, scattered wave, and free wave vs. radius;'
      DESCRP(6)='differential cross section vs. angle' 
      NGRAPH=2
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C 
C     set constant values
      PI=4.0*ATAN(1.)
      E2=14.409
      RMAX=2.0
      HBARM=7.6359
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
      MPRMPT(20)='Enter depth of potential well (eV)'
      MTAG(20)='Vzero (eV)'
      MLOLIM(20)=-5000.
      MHILIM(20)=5000.
      MREALS(20)= 50.0
C 
      MTYPE(21)=FLOAT
      MPRMPT(21)='Enter Energy (eV)'
      MTAG(21)='Energy (eV)'
      MLOLIM(21)=0.0001
      MHILIM(21)=1000.    
      MREALS(21)=20.0
C 
      MTYPE(22)=SKIP
      MREALS(22)=35
C 
      MTYPE(38)=NUM
      MPRMPT(38)='Number of integration points to r1'
      MTAG(38)='number of integration points to r1'
      MLOLIM(38)=100
      MHILIM(38)=MAXN
      MINTS(38)=400
C               
      MTYPE(39)=NUM
      MPRMPT(39)='Number of integration points between r1 and r2'
      MTAG(39)='number of extra integration points'
      MLOLIM(39)=20
      MHILIM(39)=MAXE
      MINTS(39)=60
C               
      MTYPE(40)=NUM
      MPRMPT(40)='Number of angles'
      MTAG(40)='number of angles'
      MLOLIM(40)=1
      MHILIM(40)=MAXANG
      MINTS(40)=36
C 
      MTYPE(41)=SKIP
      MREALS(41)=60
C 
      MSTRNG(MINTS(75))= 'proj4.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C                      
      MSTRNG(MINTS(86))= 'proj4.grf'
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
C Global Variables:
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P4'
C Local Variables:
      INTEGER GETINT             !get integer input from terminal
      INTEGER LDEF               !default value for L stop
C map between menu items and parameters 
      INTEGER IPOT,IZ,IVZERO,IE,INPTS,INXTRA,INANG
      PARAMETER (IPOT   = 17 )
      PARAMETER (IZ     = 18 )
      PARAMETER (IVZERO = 20 )
      PARAMETER (IE     = 21 )
      PARAMETER (INPTS  = 38 )
      PARAMETER (INXTRA = 39 )
      PARAMETER (INANG  = 40 )
C Function:
      LOGICAL LOGCVT             !converts 1 and 0 to true and false
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
C     set new parameter values   
C     physical and numerical 
      POT=MINTS(IPOT)
      Z=MINTS(IZ)
      VZERO=MREALS(IVZERO)
      E=MREALS(IE)
      NPTS=MINTS(INPTS)
      NXTRA=MINTS(INXTRA)
C
C     calculate derivative parameters:
      Z6=Z**0.166667
      K2=2*E/HBARM
      K=SQRT(K2)
      LMAX=K*RMAX/2+4
      DR=RMAX/NPTS
      RXTRA=DR*NXTRA
C
      CALL POTNTL           !fill the potential array
C
C     prompt for range of partial waves after displaying LMAX
      WRITE (OUNIT,10) LMAX
10    FORMAT (' The value for LMAX at this energy is = ',I3)
      WRITE (OUNIT,*) ' '
      LSTART=GETINT(0,0,MAXL,' Enter a value for L start')
      LDEF=MAX(LSTART,LMAX)
      LSTOP=GETINT(LDEF,LSTART,MAXL,' Enter a value for L stop')
      CALL CLEAR
C
C     calculate Legendre Polynomials if not already done
      IF (NANG .NE. MINTS(INANG)) THEN       !angles have changed
         NANG=MINTS(INANG)
         CALL LPINIT
      END IF
C     need to have Legendre polynomials all the way to LSTOP
      IF (LTABLE .LT. LSTOP) CALL LEGPOL
C
C     calculate spherical bessel functions
      CALL SBSSEL
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
C outputs parameters to MUNIT, keeping track of number of lines printed
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'IO.ALL'
       INCLUDE 'PARAM.P4'
C Passed variables:               
       INTEGER MUNIT            !unit number for output (input)
       INTEGER NLINES           !number of lines written so far (I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (MUNIT .EQ. OUNIT) CALL CLEAR
C 
       WRITE (MUNIT,2)
       WRITE (MUNIT,4)
       WRITE (MUNIT,6) E
       WRITE (MUNIT,8) LSTART,LSTOP
       WRITE (MUNIT,20) NPTS
       WRITE (MUNIT,22) NXTRA
       WRITE (MUNIT,24) NANG
C
       IF (POT .EQ. LENZ) THEN
           WRITE (MUNIT,10) Z
       ELSE IF (POT .EQ. SQUARE) THEN
           WRITE (MUNIT,12) VZERO
       ELSE IF (POT .EQ. GAUSS) THEN
           WRITE (MUNIT,14) VZERO      
       END IF
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
       NLINES=11
C
2      FORMAT (' ')
4      FORMAT (' Output from project 4: Partial-wave solution'
     +         ' to quantum scattering')
6      FORMAT (' Energy (eV) =', 1PE10.3)
8      FORMAT (' L start = ' , I3, 10X,' L stop = ', I3)
10     FORMAT (' Lenz Jensen potential with Z = ', I3)
12     FORMAT (' Square well potential with Vzero = ' 1PE10.3)
14     FORMAT (' Gaussian potential with Vzero = ' 1PE10.3)
20     FORMAT (' Number of integration points = ', I4)
22     FORMAT (' Number of points between r1 and r2 = ', I4)
24     FORMAT (' Number of angles = ', I3)
16     FORMAT (15X,'L',15X,'Delta(L)',17X,'Sigma(L)')
17     FORMAT (15X,' ',15X,'degrees ',15X,'Angstroms**2')
18     FORMAT (15X,'-',15X,'--------',15X,'------------')
C 
       RETURN
       END        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,NLINES,L,DELDEG,SIGMA)
C writes out phase shift and total cross section for each partial wave 
C to MUNIT, keeping track of number of lines printed
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P4'
C Passed variables:
      INTEGER MUNIT              !output unit specifier (input)
      INTEGER NLINES             !number of lines printed to screen(I/O)
      INTEGER L                  !partial wave number(input)
      REAL DELDEG,SIGMA          !phase shift and cross section(input)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     if screen is full, clear and type headings again
      IF ((MOD(NLINES,TRMLIN-4) .EQ. 0) 
     +                          .AND. (MUNIT .EQ. OUNIT)) THEN
         CALL PAUSE('to continue...',1)
         CALL CLEAR
         WRITE (MUNIT,16)
         WRITE (MUNIT,17)
         WRITE (MUNIT,18)
         WRITE (MUNIT,2)
         NLINES=NLINES+4
      END IF
C
      WRITE (MUNIT,20) L,DELDEG,SIGMA
C     keep track of printed lines only for terminal output
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+1
C
20    FORMAT(13X,I3,13X,1PE10.3,13X,1PE15.8)
2     FORMAT (' ')
16    FORMAT (15X,'L',15X,'Delta(L)',17X,'Sigma(L)')
17    FORMAT (15X,' ',15X,'degrees ',15X,'Angstroms**2')
18    FORMAT (15X,'-',15X,'--------',15X,'------------')
C
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DELOUT(DELDEG,SIGMA)
C outputs delta and partial wave cross sections to screen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P4'
C Passed variables:               
      REAL DELDEG(0:MAXL)         !phase shift in degrees (input)
      REAL SIGMA(0:MAXL),SIGTOT   !total cross section (input)
      INTEGER L                   !partial wave index (input)
C Local variables:
      INTEGER NLINES              !number of lines written so far (I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       CALL PRMOUT(OUNIT,NLINES)
       DO 100 L=LSTART,LSTOP
          CALL TXTOUT(OUNIT,NLINES,L,DELDEG(L),SIGMA(L))
100    CONTINUE
       RETURN
       END    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SIGOUT(MUNIT,F,DSIGMA,SIGTOT)
C writes out scattering angle, scattering amplitude, 
C differential cross section and total cross section to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P4'
C Input variables:
      INTEGER MUNIT                 !output unit specifier 
      REAL DSIGMA(0:MAXANG)         !differential cross section
      COMPLEX F(0:MAXANG)           !scattering amplitude
      INTEGER ITHETA                !index of theta
      REAL SIGTOT                   !total cross section
C Local variables:
      INTEGER NLINES                !keep track of lines printed out
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE (MUNIT,2)
      WRITE (MUNIT,2)
      WRITE (MUNIT,4) SIGTOT
      WRITE (MUNIT,2)
      WRITE (MUNIT,2)
      WRITE (MUNIT,10)
      WRITE (MUNIT,11)
      WRITE (MUNIT,12)
C
      NLINES=8
C
C     write out data, allowing time for user to examine full pages
      DO 100 ITHETA=0,NANG
        WRITE (MUNIT,14) DEGREE(ITHETA),F(ITHETA),DSIGMA(ITHETA)
        NLINES=NLINES+1
        IF ((MOD(NLINES,TRMLIN-4) .EQ. 0) .AND. 
     +               (MUNIT .EQ. OUNIT)) THEN
           CALL PAUSE('to continue...',1)
           CALL CLEAR
        END IF
100   CONTINUE
C
2      FORMAT (' ')
4      FORMAT (15X,'TOTAL CROSS SECTION = ',1PE15.8,' Angstroms**2')
10     FORMAT (10X,'Theta' ,14X,'Amplitude (Re,Im)',15X,'dSigma/dTheta')
11     FORMAT (9X,'degrees',13X,'  Angstroms**2   ',15X,' Angstroms**2')
12     FORMAT (9X,'-------',13X,'-----------------',15X,'-------------')
14     FORMAT (7X,F8.2,6X,'(',1PE15.8,','1PE15.8,')',6X,1PE15.8)
C
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFOUT(DEVICE,L,DELDEG,SIGMA)
C outputs potential and both free and scattered wave functions for one L:
C does not allow for hardcopy or file output
C in order to avoid the creation of voluminous files
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables                                                   
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P4'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE               !which device is being used?
      INTEGER L                    !partial wave
      REAL DELDEG,SIGMA            !phase shift and cross section
      REAL PSI(MAXN+MAXE),PSIF(MAXN+MAXE)   !wave functions
      REAL PSIMAX,PSFMAX           !maxima of wave functions
      REAL VEFF(MAXN+MAXE)         !effective pot for fixed L
      COMMON/RESULT/PSI,PSIF,PSIMAX,PSFMAX,VEFF!pass results from Numerv
C Local variables
      INTEGER IR                   !index of lattice
      REAL VGRF(MAXE+MAXN)         !V scaled for graphing
      REAL PSIG(MAXE+MAXN),PSIFG(MAXE+MAXN)  !rescaled wave functions
      REAL ENERGY(2),RADIUS(2)     !array for plotting energy  
      CHARACTER*9 CL,CDEL,CSIG,CE  !character version of data
      INTEGER LEN,LCDEL,LCSIG      !length of strings
      INTEGER SCREEN               !send to terminal
      INTEGER PAPER                !make a hardcopy
      INTEGER FILE                 !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     setup graphing paramters
      IF (DEVICE .NE. FILE) THEN
          NPLOT=1                        !how many plots?
          IPLOT=1
C
          XMIN=0.0                       !axis parameters
          XMAX=RMAX+RXTRA            
          Y0VAL=XMIN
          YMIN=-2*E                      !energy axis varies between
          YMAX=+4*E                      !-2E and 4E
          X0VAL=0.
C
          NPOINT=NPTS+NXTRA              
C
          LABEL(1)= 'radius (Angstroms)' !titles and labels
          LABEL(2)= 'Veff(eV), scattered wave (solid),'
     +              //' and free (ooo) wave'
          CALL CONVRT(E,CE,LEN)
          IF (POT .EQ. LENZ) THEN
            TITLE = ' Lenz-Jensen Potential, Energy (eV)='//CE
          ELSE IF (POT .EQ. SQUARE) THEN
            TITLE = 'Square Well Potential, Energy (eV)='//CE
          ELSE IF (POT .EQ. GAUSS) THEN
            TITLE = 'Gaussian Potential, Energy (eV)='//CE
          END IF
C
C         arrays for plotting energy
          RADIUS(1)=0.
          RADIUS(2)=XMAX
          ENERGY(1)=E
          ENERGY(2)=E
C 
C         plot Veffective between -2E and 4E
          DO 30 IR=1,NPTS+NXTRA         !normalize VEFF
             VGRF(IR)=VEFF(IR)
             IF (VGRF(IR) .GT. YMAX) VGRF(IR)=YMAX
             IF (VGRF(IR) .LT. YMIN) VGRF(IR)=YMIN
30        CONTINUE
C 
C         rescale wave functions so that it fits on the energy scale
          DO 10 IR=1,NPTS+NXTRA         
            PSIG(IR)=PSI(IR)*3*E/PSIMAX+E
            PSIFG(IR)=PSIF(IR)*3*E/PSFMAX+E
10        CONTINUE
C                            
          ILINE=1                      !line and symbol styles
          ISYM=1
          IFREQ=0
          NXTICK=5
          NYTICK=5
C
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
C             
C         create the legend
          CALL ICNVRT(L,CL,LEN)
          CALL CONVRT(DELDEG,CDEL,LCDEL)
          CALL CONVRT(SIGMA,CSIG,LCSIG)
          INFO='Delta(L=' //CL(1:LEN) //')=' //CDEL(1:LCDEL)//
     +        ' deg'//
     +        '  Sig(L=' //CL(1:LEN) //')='//CSIG(1:LCSIG)//' A**2'
C
          CALL LNLNAX                          
C
          CALL XYPLOT (R,VGRF)                   !plot potential
          NPOINT=2
          CALL XYPLOT (RADIUS,ENERGY)            !plot energy
C
          NPOINT=NPTS+NXTRA
          CALL XYPLOT (R,PSIG)                   !plot wave function     
          IFREQ=4             
          CALL XYPLOT (R,PSIFG)                  !plot free wave function
C
      END IF                 
C 
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)  !end graphics package
      IF (DEVICE .EQ. SCREEN) CALL TMODE        !switch to text mode
C
100   FORMAT (/,' Patience, please; output going to a file.')
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFSIG(DEVICE,DSIGMA,SIGTOT)
C graphs differential cross section vs. angle for all partial waves
C on a semi-log scale
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables                                                   
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P4'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE         !which device is being used?
      REAL SIGTOT            !total cross section
      REAL DSIGMA(0:MAXANG)  !differential cross section
C Local variables
      INTEGER ITHETA         !indexes angle                   
      INTEGER EXPMAX,EXPMIN  !min and max exp for diff cross section
      CHARACTER*9 CE,CSIG    !energy,sigma as a character string
      INTEGER LEN            !string length
      INTEGER SCREEN         !send to terminal
      INTEGER PAPER          !make a hardcopy
      INTEGER FILE           !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     calculate parameters for graphing
      IF (DEVICE .NE. FILE) THEN
C
          NPLOT=1                        !how many plots?
          IPLOT=1
C 
          YMAX=0.                        !find limits on data points
          YMIN=DSIGMA(1)
          DO 20 ITHETA=0,NANG
             IF (DSIGMA(ITHETA) .GT. YMAX) YMAX=DSIGMA(ITHETA)
             IF (DSIGMA(ITHETA) .LT. YMIN) YMIN=DSIGMA(ITHETA)
20        CONTINUE
          !find integer limits on exponent
          EXPMAX=INT(LOG10(YMAX))
          IF (YMAX .GT. 1.) EXPMAX =EXPMAX+1
          EXPMIN=INT(LOG10(YMIN))
          IF (YMIN .LT. 1.) EXPMIN=EXPMIN-1
          YMAX=10.**EXPMAX
          YMIN=10.**EXPMIN
C 
          XMIN=DEGREE(0)
          XMAX=DEGREE(NANG)
          Y0VAL=XMIN                       
          X0VAL=YMIN
C          
          NPOINT=NANG
C 
          ILINE=1                      !line and symbol styles
          ISYM=4
          IFREQ=1
          NXTICK=4
          NYTICK=EXPMAX-EXPMIN
          IF (NYTICK .GT. 8) THEN      !keep number of ticks small
             IF (MOD(NYTICK,2) .EQ. 0) THEN
                NYTICK=NYTICK/2
             ELSE
                NYTICK=8
             END IF
          END IF
C 
          CALL CONVRT(E,CE,LEN)             !titles and labels
          IF (POT .EQ. LENZ) THEN
            TITLE = ' Lenz-Jensen Potential, Energy (eV)='//CE
          ELSE IF (POT .EQ. SQUARE) THEN
            TITLE = 'Square Well Potential, Energy (eV)='//CE
          ELSE IF (POT .EQ. GAUSS) THEN
            TITLE = 'Gaussian Potential, Energy (eV)='//CE
          END IF
          CALL CONVRT(SIGTOT,CSIG,LEN)
          INFO='Total cross section in A**2 = '//CSIG
          LABEL(1)= 'Angle (degrees)'
          LABEL(2)= 'Differential Cross Section (Angstroms**2)'
C 
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
          CALL LGLNAX                          !draw axes
      END IF
C                                                      
C     output results
      IF (DEVICE .EQ. FILE) THEN
         WRITE (GUNIT,10)
         WRITE (GUNIT,11)
         WRITE (GUNIT,12)
         WRITE (GUNIT,14) (DEGREE(ITHETA),DSIGMA(ITHETA),ITHETA=0,NANG)
      ELSE          
         CALL XYPLOT (DEGREE,DSIGMA)             
      END IF
C                    
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE) !end graphics package
      IF (DEVICE .EQ. SCREEN) CALL TMODE       !switch to text mode
C
10    FORMAT (21X,'Theta',  21X,'dSigma/dTheta')
11    FORMAT (20X,'degrees',20X,' Angstroms**2')             
12    FORMAT (20X,'-------',20X,'-------------')             
14    FORMAT (18X,1PE10.3,18X,1PE15.8)
100   FORMAT (/,' Patience, please; output going to a file.')
C
      RETURN
      END                   
                     
