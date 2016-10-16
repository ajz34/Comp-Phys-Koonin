CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM EXMPL2
C     Example 2:  Trajectories in the Henon-Heiles potential  
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company, Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop/ execute once for each set of param
        CALL PARAM        !get input from screen
        CALL ARCHON       !calculates surface of section for one energy
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON
C calculates surface of section (SOS) for one energy
C allows for several sets of initial conditions
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E2'
C Local variables:
      REAL SSY(MAXSS),SSPY(MAXSS)  !surface of section points
      INTEGER NCROSS               !total number of SOS crossings     
      INTEGER ANOTHR               !start with another set of init cond?
      INTEGER DEF                  !default number of SOS points
      INTEGER DEFMIN               !minimum number of SOS points
      INTEGER MORE                 !how many more SOS points?
      INTEGER SCREEN               !send to terminal
      INTEGER PAPER                !make a hardcopy
      INTEGER FILE                 !send to a file
C Functions:
      INTEGER YESNO                !yes or no input 
      INTEGER GETINT               !integer input 
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NCROSS=0                     !no surface of section points yet
50    CONTINUE                     !loop over initial conditions
         CALL SOS(NCROSS,SSPY,SSY) !calculate surface of section
      CALL CLEAR                   !allow for another set of init cond
      ANOTHR=YESNO(0,'Do you want another set of initial conditions?')
      IF (ANOTHR .EQ. 1) THEN
         DEF=MIN(100,MAXSS-NSOS)   !don't exceed storage space
         DEFMIN=MIN(10,MAXSS-NSOS)
         MORE=GETINT(DEF,DEFMIN,MAXSS-NSOS,
     +        'How many surface of section points?')
         NSOS=NSOS+MORE
      GOTO 50                          
      END IF
C 
      IF (GHRDCP) CALL GRFOUT(PAPER,SSY,SSPY,1)  !hrdcpy once per energy
C      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOS(NCROSS,SSPY,SSY)
C finds the surface of section points (SSPY,SSY) 
C for a single set of initial conditions;
C NCROSS keeps track of the total number of SOS points found
C for all initial conditions
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E2'
C Passed variables:
      REAL SSY(MAXSS),SSPY(MAXSS) !surface of section points (output)
      INTEGER NCROSS              !total number of SOS crossings (I/O)    
C Local variables:
      REAL T                      !time
      REAL VAR(4)                 !coordinates and momenta
      REAL OLDVAR(4)              !previous values of VAR       
      INTEGER NLINES              !number of lines written to terminal
      INTEGER NBEGIN     !index of first SOS point for current init cond
      LOGICAL CROSS               !have we just crossed a SOS?
      REAL PYMAX                  !limits on PY for this E
      INTEGER MORE                !how many more SOS points?
      INTEGER IVAR                !dependent variable index
      REAL TSTOP                  !time at last stop in integration
      INTEGER SCREEN              !send to terminal
      INTEGER PAPER               !make a hardcopy
      INTEGER FILE                !send to a file
C Functions:
      INTEGER GETINT,YESNO        !get integer,boolean input from screen
      REAL V                      !potential function
      REAL GETFLT                 !get real input from screen
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     prompt for initial Y and PY; obtain PX from Energy conservation
      WRITE (OUNIT,*) ' '
      WRITE (OUNIT,*) ' Input initial conditions (recall X=0):'
      WRITE (OUNIT,*) ' '
      VAR(1)=0.
      VAR(2)=GETFLT(0.,CMIN,CMAX,' Enter initial value for y')
      PYMAX=SQRT(2*(EINIT-V(VAR(1),VAR(2))))
      VAR(4)=GETFLT(0.,-PYMAX,PYMAX,' Enter initial value for py')
      VAR(3)=SQRT(2*(EINIT-V(VAR(1),VAR(2)))-VAR(4)**2 )
C
C     output summary of parameters
      IF ((TTERM) .AND. (.NOT. GTERM)) CALL PRMOUT(OUNIT,NLINES,VAR)
      IF (TFILE) CALL PRMOUT(TUNIT,NLINES,VAR)
      IF (GFILE) CALL PRMOUT(GUNIT,NLINES,VAR)
C                                                         
C     setup initial values
      DO 500 IVAR=1,4
         OLDVAR(IVAR)=VAR(IVAR)
500   CONTINUE
      T=0.                          
      TSTOP=T
      NBEGIN=NCROSS+1
      IF (TRAJCT) CALL TRJINT(TDEV)   !prepare for graphing traj
C          
C     integrate until we get requested number of SOS crossings
10    CONTINUE                               
         CALL STEP(T,VAR,OLDVAR,NCROSS,SSY,SSPY,NLINES,TSTOP)
      IF (NCROSS .LT. NSOS) GOTO 10 
C
      IF (.NOT. GTERM) THEN
         IF (TTERM) CALL PAUSE('to continue...',1)
         IF (TTERM) CALL CLEAR
      ELSE
          CALL GRFOUT(SCREEN,SSY,SSPY,NBEGIN)   !graph SOS points
      END IF
C      
C     see if user wishes more SOS points
      CALL CLEAR
      MORE=GETINT(0,0,MAXSS-NSOS,'How many more SOS points?')
      NSOS=NSOS+MORE
      IF (MORE .GT. 0) THEN
        NLINES=0
        IF (TRAJCT) CALL TRJINT(TDEV)
        TSTOP=T
      GOTO 10
      END IF
C
C     graphics output to a file, once per initial condition
      IF (GFILE) CALL GRFOUT(FILE,SSY,SSPY,NBEGIN)
C      
      RETURN                                        
      END          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE STEP(T,VAR,OLDVAR,NCROSS,SSY,SSPY,NLINES,TSTOP)
C integrates one time step and checks for SOS crossing
C finds SSY and SSPY on surface of section
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'                                          
      INCLUDE 'PARAM.E2'
C Passed variables:
      REAL SSY(MAXSS),SSPY(MAXSS) !surface of section points (output)
      INTEGER NCROSS              !total number of SOS crossings (I/O)    
      REAL T                      !time (I/O)
      REAL VAR(4)                 !coordinates and momenta (I/O)
      REAL OLDVAR(4)              !previous values of VAR (I/O)
      INTEGER NLINES              !num of lines written to terminal(I/O)
      REAL TSTOP                  !time at last stop in integ (input)
C Local variables:
      LOGICAL CROSS               !have we just crossed a SOS?             
      REAL E,EPOT,EKIN            !energies
      INTEGER ITIME               !number of time steps
      INTEGER IVAR                !dependent function index
C Function:
      REAL V                      !value of the potential 
      INTEGER SCREEN              !send to terminal
      INTEGER PAPER               !make a hardcopy
      INTEGER FILE                !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CROSS=.FALSE.
      T=T+TSTEP
      CALL RUNGE(CROSS,VAR,OLDVAR)             !take a Runge-Kutta step
C                                        
      EPOT=V(VAR(1),VAR(2))                    !calculate energies    
      EKIN=(VAR(3)**2+VAR(4)**2)/2.0
      E=EKIN+EPOT
C
      IF (TRAJCT) THEN                 !output trajectories, but only
        ITIME=(T-TSTOP)/TSTEP          !for a limited time
        IF ((ITIME .EQ. NTRJ+1) .AND. (TDEV .EQ. SCREEN)) THEN
            CALL NOTICE                !let user know we're still integ
        ELSE IF (ITIME .LE. NTRJ) THEN
            CALL TRJOUT(VAR,OLDVAR)    !output trajectories
        END IF
      END IF
C 
      CROSS=(VAR(1)*OLDVAR(1)) .LT. 0.0        !check for sos crossing
C                                      
      DO 100 IVAR=1,4                          !update variables
         OLDVAR(IVAR)=VAR(IVAR)
100   CONTINUE
C
      IF (CROSS) THEN              !find SOS crossing and output results
          NCROSS=NCROSS+1
          CALL RUNGE(CROSS,VAR,OLDVAR)         !find SSY and SSPY 
          SSY(NCROSS)=VAR(2)                   !store SOS intersections
          SSPY(NCROSS)=VAR(4)
C
          DO 200 IVAR=1,4                      !reset old values
             VAR(IVAR)=OLDVAR(IVAR)
200       CONTINUE
C
          IF ((TDEV .NE. SCREEN) .AND. (TTERM)) 
     +        CALL TXTOUT(OUNIT,T,E,EPOT,EKIN,NLINES)  !text output
          IF (TFILE) CALL TXTOUT(TUNIT,T,E,EPOT,EKIN,NLINES)
      ENDIF
C             
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RUNGE(CROSS,VAR,OLDVAR)
C uses 4th order Runge-Kutta algorithm to integrate the 
C four coupled ODE's; 
C the independent variable is determined by the value of CROSS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E2'
C Passed variables:
      LOGICAL CROSS                !is this to find sos crossing?(input)
      REAL VAR(4)                  !coordinates and momenta (output)
      REAL OLDVAR(4)               !previous values of VAR(input)
C Local variables:
      REAL F(4)                    !derivatives
      REAL K1(4),K2(4),K3(4),K4(4) !increments 
      REAL H                       !step size           
      INTEGER I                    !dependent variable index
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      H=TSTEP
      IF (CROSS) H=-VAR(1)
C      
      CALL EVAL(F,VAR,CROSS)
      DO 100 I=1,4
         K1(I)=H*F(I)
         VAR(I)=OLDVAR(I)+K1(I)/2.0
100   CONTINUE
C
      CALL EVAL(F,VAR,CROSS)
      DO 200 I=1,4
         K2(I)=H*F(I)
         VAR(I)=OLDVAR(I)+K2(I)/2.0
200   CONTINUE
C
      CALL EVAL(F,VAR,CROSS)
      DO 300 I=1,4
         K3(I)=H*F(I)
         VAR(I)=OLDVAR(I)+K3(I)
300   CONTINUE
C
      CALL EVAL(F,VAR,CROSS)
      DO 400 I=1,4
         K4(I)=H*F(I)
         VAR(I)=OLDVAR(I)+(K1(I)+2*K2(I)+2*K3(I)+K4(I))/6
400   CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE EVAL(F,VAR,CROSS)
C Evaluate derivatives F evaluated at VAR
C The independent variable is determined by CROSS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
      REAL VAR(4)              !variables (input)
      REAL F(4)                !derivatives of variables (output)
      LOGICAL CROSS            !is this for a surface of section?(input)
C Local variables:               
      REAL DENOM               !factor to calc SOS
C Functions:                  
      REAL XDERIV,YDERIV       !x and y derivatives of the potential
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     to find surface of section, all derivatives are divided by PX
      DENOM=1.0
      IF (CROSS) DENOM=VAR(3)
C
      F(1)=VAR(3)/DENOM
      F(2)=VAR(4)/DENOM                  
      F(3)=-1.0*XDERIV(VAR(1),VAR(2))/DENOM
      F(4)=-1.0*YDERIV(VAR(1),VAR(2))/DENOM
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION V(X,Y)
C Calculates the potential and forces
C
C If you change the potential, you may also need to change DY0 and TOLY
C for Y limit searches, as well as limits for X and Y in TRJINT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
      REAL X,Y                     !coordinates
C Functions:
      REAL XDERIV,YDERIV           !x and y derivatives of the potential
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      V=(X**2+Y**2)/2+X**2*Y-Y**3/3
      RETURN
C
      ENTRY XDERIV(X,Y)
      XDERIV=X+2*X*Y
      RETURN
C
      ENTRY YDERIV(X,Y)
      YDERIV=Y+X**2-Y**2
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LIMITS
C Find limits on Y from energy conservation using a simple search
C This limit is on the surface of section where X=0.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E2'
C Local variables:
      REAL DY                      !step in simple search
      REAL V                       !potential (function)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DY=DY0                       !search for CMAX starting at 0
      CMAX=0.0
10    CONTINUE
         CMAX=CMAX+DY              
         IF (V(0.0,CMAX) .GT. EINIT) THEN      !CMAX is where all the
            CMAX=CMAX-DY                       !energy is potntl energy
            DY=DY/2.                           !recall that X=0.0 on SOS
         ENDIF
      IF (DY.GE.TOLY) GOTO 10
C
      DY=DY0                       !search for CMIN starting at 0
      CMIN=0.0
20    CONTINUE
         CMIN=CMIN-DY
         IF ( V(0.0,CMIN) .GT. EINIT) THEN     !CMIN is where all the 
            CMIN=CMIN+DY                       !energy is potntl energy
            DY=DY/2.                           !recall that X=0.0 on SOS
         ENDIF
      IF (DY.GE.TOLY) GOTO 20
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
      INCLUDE 'PARAM.E2'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL SETUP                !get environment parameters
C 
C     display header screen     
      DESCRP(1)= 'EXAMPLE 2'
      DESCRP(2)= 'Trajectories in the Henon-Heiles potential'
      NHEAD=2
C
C     text output description
      DESCRP(3)= 'kinetic, potential, and total energy;'
      DESCRP(4)= 
     + 'and percent change in energy at surface of section crossing'
      NTEXT=2
C
C     graphics output description
      DESCRP(5)= 'trajectory and surface of section'
      NGRAPH=1
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C       
C     calculate constants  (used in search for limits on y)
      DY0=0.1
      TOLY=.0005 
C
      CALL MENU                        !constant part of menu          
C
      MTYPE(13)=FLOAT
      MPRMPT(13)= 'Enter Energy'
      MTAG(13)= 'Energy'
      MLOLIM(13)=0.0
      MHILIM(13)=1.0/6.
      MREALS(13)=.1
C                      
      MTYPE(14)=SKIP
      MREALS(14)=35.
C                                      
      MTYPE(38)=FLOAT
      MPRMPT(38)= 'Enter time step'
      MTAG(38)= 'Time step'
      MLOLIM(38)=.00001
      MHILIM(38)=.2
      MREALS(38)=.12
C 
      MTYPE(39)=NUM
      MPRMPT(39)= 'Enter number of surface of section points'
      MTAG(39)= 'Number of surface of section points'
      MLOLIM(39)=10
      MHILIM(39)=MAXSS
      MINTS(39)=100
C            
      MTYPE(40)=SKIP
      MREALS(40)=60.
C 
      MSTRNG(MINTS(75))= 'exmpl2.txt'
C
      MTYPE(76)=SKIP
      MREALS(76)=80.
C 
      MSTRNG(MINTS(86))= 'exmpl2.grf'
C 
      MTYPE(87)=NUM
      MPRMPT(87)= 'Enter number of points in trajectory to plot'
      MTAG(87)= 'Number of points in trajectory'
      MLOLIM(87)=0
      MHILIM(87)=10000
      MINTS(87)=4000
C 
      MTYPE(88)=SKIP
      MREALS(88)=90.
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
      INCLUDE 'PARAM.E2'
C Local variables:
      INTEGER SCREEN              !send to terminal
      INTEGER PAPER               !make a hardcopy
      INTEGER FILE                !send to a file
C map between menu items and parameters
      INTEGER IE,ITSTEP,INSOS,INTRJ
      PARAMETER (IE = 13 )
      PARAMETER (ITSTEP =38 )
      PARAMETER (INSOS =39 )
      PARAMETER (INTRJ  = 87 )
C Function:
      LOGICAL LOGCVT              !converts 1 and 0 to true and false
      DATA SCREEN,PAPER,FILE/1,2,3/
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
      EINIT=MREALS(IE)
      TSTEP=MREALS(ITSTEP)
      NSOS=MINTS(INSOS)
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
      NTRJ=MINTS(INTRJ)
C             
C     trajectories are output ONLY if graphics are available
C      (it's too much data to output to a file)
C     since traj can only go to one device, these are sent to the screen
      TRAJCT=.FALSE.
      TDEV=0
      IF (GTERM) THEN
         TDEV=SCREEN
         TRAJCT=.TRUE.
      ELSE IF (GHRDCP) THEN
         TDEV=PAPER
         TRAJCT=.TRUE.
      END IF 
C
C     open files 
      IF (TFILE) CALL FLOPEN(TNAME,TUNIT)
      IF (GFILE) CALL FLOPEN(GNAME,GUNIT)
C     files may have been renamed
      MSTRNG(MINTS(ITNAME))=TNAME
      MSTRNG(MINTS(IGNAME))=GNAME
C
C     calculate CMIN and CMAX (limits on Y on the surface of section)
      CALL LIMITS
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE PRMOUT(MUNIT,NLINES,VAR)
C outputs parameter summary to the specified unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:   
       INCLUDE 'IO.ALL'
       INCLUDE 'PARAM.E2'
C Passed variables:               
       INTEGER MUNIT            !unit number for output (input)
       INTEGER NLINES           !number of lines written so far (I/O)
       REAL VAR(4)              !initial values of coord and momenta (I)
C Local variables: 
       INTEGER I                !independent variable index
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (MUNIT .EQ. OUNIT) CALL CLEAR
C 
       WRITE (MUNIT,2)
       WRITE (MUNIT,4)
       WRITE (MUNIT,6) TSTEP
       WRITE (MUNIT,10) EINIT
       WRITE (MUNIT,8) CMIN,CMAX,SQRT(2*EINIT)
       WRITE (MUNIT,12) (VAR(I),I=1,4)
       WRITE (MUNIT,2)
C
C      different header for text and graphics files
       IF (MUNIT .EQ. GUNIT) THEN
         WRITE (MUNIT,20)
         WRITE (MUNIT,25)
       ELSE
         WRITE (MUNIT,30)
         WRITE (MUNIT,35)
         WRITE (MUNIT,40)
         WRITE (MUNIT,2)
       END IF
C 
       NLINES=11
C
2      FORMAT (' ')
4      FORMAT (' Output from example 2: Trajectories in the Henon-',
     +         'Heiles Potential')
6      FORMAT (' Time step =', E12.5)
10     FORMAT (' Energy =',F6.3) 
8      FORMAT (' Ymin =',F6.3, 5X,' Ymax =',F6.3,5X,'Pymax =',F6.3)
12     FORMAT (' Xinit=',F6.3, 5X,' Yinit=',F6.3,5X,'PXinit=',F6.3,
     +         5X,'PYinit=',F6.3)
20     FORMAT (7X,'Y on SOS',7X, 'PY on SOS')
25     FORMAT (7X,'--------',7X, '---------')
30     FORMAT (9X,'Time',9X,'Kinetic',4X,'Potential',5X,'Total',
     +         7X,'Percent') 
35     FORMAT (23X,'Energy',6x,'Energy',6X,'Energy',7X,'Change')
40     FORMAT (9X,'----',9X,'-------',4X,'---------',6X,'-----',
     +         6X,'-------') 
C 
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,T,E,EPOT,EKIN,NLINES)
C writes out energy data at each sos crossing
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E2'
C Passed variables:
      INTEGER MUNIT              !output unit specifier (I)
      INTEGER NLINES             !number of lines printed to screen (I/O)
      REAL E,EPOT,EKIN           !energies (I)
      REAL T                     !time (I)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     if this is a new page, retype headings
      IF ((NLINES .EQ. 0) .AND. (MUNIT .EQ. OUNIT)) THEN
         CALL CLEAR
         WRITE (MUNIT,30)
         WRITE (MUNIT,35)
         WRITE (MUNIT,40)
         WRITE (MUNIT,2)
         NLINES=NLINES+4
C     else if screen is full, clear screen and retype headings          
      ELSE IF ((MOD(NLINES,TRMLIN-4) .EQ. 0) 
     +                          .AND. (MUNIT .EQ. OUNIT)) THEN
         CALL PAUSE('to continue...',1)
         CALL CLEAR
         WRITE (MUNIT,30)
         WRITE (MUNIT,35)
         WRITE (MUNIT,40)
         WRITE (MUNIT,2)
         NLINES=NLINES+4
      END IF
C
      WRITE (MUNIT,20)T,EKIN,EPOT,E,ABS((E-EINIT)/EINIT)
C                                   
C     keep track of printed lines only for terminal output
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+1
C
20     FORMAT (5X,1PE12.5,3(5X,0PF7.5),5X,1PE10.3)
2      FORMAT (' ')
30     FORMAT (9X,'Time',9X,'Kinetic',4X,'Potential',5X,'Total',
     +         7X,'Percent') 
35     FORMAT (23X,'Energy',6x,'Energy',6X,'Energy',7X,'Change')
40     FORMAT (9X,'----',9X,'-------',4X,'---------',6X,'-----',
     +         6X,'-------') 
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE GRFOUT(DEVICE,SSY,SSPY,NBEGIN)
C outputs surface of section from NBEGIN to NSOS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global parameters:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E2'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE                  !which device is being used?
      REAL SSY(MAXSS),SSPY(MAXSS)     !Y and PY values on SOS
      INTEGER NBEGIN                  !beginning SOS for these init cond
C Local parameters:                         
      INTEGER I                       !indexes SOS points
      INTEGER SCREEN                  !send to terminal
      INTEGER PAPER                   !make a hardcopy
      INTEGER FILE                    !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     calculate parameters for graphing
      IF (DEVICE .NE. FILE) THEN
C 
          NPLOT=2                       !how many plots
          IPLOT=2
C                                         
C         if both screen and paper are used , only SOS are sent to paper
          IF ((GTERM) .AND. (DEVICE .EQ. PAPER)) THEN
             NPLOT=1
             IPLOT=1
             CALL GTDEV(DEVICE)                   !device nomination
          END IF
C 
          YMAX=SQRT(2*EINIT)                      !horiz axis is y
          YMIN=-YMAX                              !vert axis is py
          XMAX=CMAX 
          XMIN=CMIN 
          Y0VAL=CMIN
          X0VAL=0.
C 
          NPOINT=NSOS
C 
          ILINE=5                   !line and symbol styles
          ISYM=5                    !this choice gives unconnected dots
          IFREQ=1
          NXTICK=7
          NYTICK=7
C 
          INFO=' '
          LABEL(1)='Y'
          LABEL(2)='PY'
C 
          CALL LNLNAX                          !draw axes
      END IF                                   
C         
C     output results
      IF (DEVICE .EQ. FILE) THEN
          WRITE (GUNIT,70) (SSY(I),SSPY(I),I=NBEGIN,NSOS) 
      ELSE  
          CALL XYPLOT(SSY,SSPY)
      END IF
C 
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE) !close graphing package
      IF (DEVICE .EQ. SCREEN) CALL TMODE     !switch back to text mode
C              
70    FORMAT (2(5X,E11.3))
100   FORMAT (/,' Patience, please; output going to a file.')     
C 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE TRJINT(DEVICE)
C prepares to graph trajectories
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global parameters:
      INCLUDE 'IO.ALL'
      INCLUDE 'GRFDAT.ALL'
      INCLUDE 'PARAM.E2'
C Input variables:
      INTEGER DEVICE                  !which device is being used?
C Local parameters:                         
      CHARACTER *9 CE                 !energy as a character string
      INTEGER LEN                     !length of string
      REAL X(4),Y(4)                  !corners of the triangle
      INTEGER SCREEN                  !send to terminal
      INTEGER PAPER                   !make a hardcopy
      INTEGER FILE                    !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NPLOT=2                       !how many plots
      IPLOT=1
C 
      YMAX=1.                       !corners of bounding triangle
      YMIN=-.5
      XMAX=SQRT(3.)/2. 
      XMIN=-XMAX 
      Y0VAL=0. 
      X0VAL=0.
C 
      NPOINT=4
C 
      ILINE=1                       !line and symbol styles
      ISYM=1 
      IFREQ=0
      NXTICK=4
      NYTICK=4
C
      CALL CONVRT(EINIT,CE,LEN)
      TITLE='Henon-Heiles Potential, Energy='//CE
      LABEL(1)='X'
      LABEL(2)='Y'
C 
      CALL GTDEV(DEVICE)                   !device nomination
      IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
      CALL LNLNAX                          !draw axes
C             
      X(1)=XMAX                            !draw the bounding triangle
      Y(1)=YMIN
      X(2)=XMIN
      Y(2)=YMIN
      X(3)=0.
      Y(3)=YMAX
      X(4)=X(1)
      Y(4)=Y(1)
      CALL XYPLOT(X,Y)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE TRJOUT(VAR,OLDVAR)
C outputs trajectory, one line segment per call
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global parameters:
      INCLUDE 'PARAM.E2'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      REAL VAR(4)                     !coordinates and momenta
      REAL OLDVAR(4)                  !previous values of VAR
C Local parameters:                         
      REAL X(2),Y(2)                  !coordinates 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NPOINT=2
      X(1)=OLDVAR(1)
      Y(1)=OLDVAR(2)
      X(2)=VAR(1)
      Y(2)=VAR(2)
      CALL XYPLOT(X,Y)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NOTICE
C let user know we're still computing, even though trajectory isn't
C being plotted
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'GRFDAT.ALL'
      INFO=' Still computing, but no longer plotting ...'
      CALL LEGEND
      RETURN
      END
