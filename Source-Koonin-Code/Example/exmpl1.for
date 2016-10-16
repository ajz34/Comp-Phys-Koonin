CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM EXMPL1
C     Example 1: Bohr-Sommerfeld quantization for bound states of the
C               Lennard-Jones Potential
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company, Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop/ execute once for each set of param
        CALL PARAM        !get input from screen
        CALL ARCHON       !search for bound states
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON  
C finds the bound states of the Lennard-Jones potential
C from the Bohr-Sommerfeld quantization rule
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E1'
C Local variables:
      REAL S                        !current value of action
      REAL E1                       !current value of energy
      REAL X1,X2                    !current turning points
      REAL F1                       !f=action/2 - pi/2 -ilevel*pi
      INTEGER ILEVEL                !current level 
      REAL ENERGY(0:MAXLVL)         !energy of bound state
      REAL XIN(0:MAXLVL)            !inner turning point
      REAL XOUT(0:MAXLVL)           !outer turning point
      INTEGER NLINES                !number of lines printed to terminal
      INTEGER SCREEN                !send to terminal
      INTEGER PAPER                 !make a hardcopy
      INTEGER FILE                  !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     output summary of parameters
      IF (TTERM) CALL PRMOUT(OUNIT,NLINES)
      IF (TFILE) CALL PRMOUT(TUNIT,NLINES)
      IF (GFILE) CALL PRMOUT(GUNIT,NLINES)
C      
C     search for bound states
      E1=-1.                      !begin at the well bottom
      F1=-PI/2                    !the action is zero there     
C
C     find the NLEVEL bound states
      DO 100 ILEVEL=0,NLEVEL-1
C
         CALL SEARCH(ILEVEL,E1,F1,X1,X2) !search for eigenvalue 
         ENERGY(ILEVEL)=E1        !store values for this state
         XIN(ILEVEL)=X1
         XOUT(ILEVEL)=X2
C
C        text output                           
         IF (TTERM) CALL TXTOUT(OUNIT,ILEVEL,E1,X1,X2,NLINES)
         IF (TFILE) CALL TXTOUT(TUNIT,ILEVEL,E1,X1,X2,NLINES)
C       
         F1=F1-PI                 !guess to begin search for next level
C
100   CONTINUE                             
C                                                                   
      IF (TTERM) CALL PAUSE('to continue...',1) 
      IF (TTERM) CALL CLEAR
C
C     graphics output
      IF (GTERM) CALL GRFOUT(SCREEN,ENERGY,XIN,XOUT)
      IF (GFILE) CALL GRFOUT(FILE,ENERGY,XIN,XOUT)
      IF (GHRDCP) CALL GRFOUT(PAPER,ENERGY,XIN,XOUT)
C      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE SEARCH(N,E1,F1,X1,X2)
C finds the N'th bound state 
C E1 is passed in as initial guess for the bound state energy
C    and returned as the true bound state energy with turning points
C    X1 and X2
C F1 is the function which goes to zero at a bound state
C    F1 = action/2-(n+1/2)*pi         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E1'
C Input/Output variables:
      INTEGER N                !current level (input)
      REAL E1,E2               !trial energies (I/O)
      REAL F1,F2               !f=action/2-pi*(n+1/2) (I/O)
      REAL S                   !action (output)
      REAL X1,X2               !turning points (output)
C Local variables:      
      REAL DE                  !increment in energy search
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     guess the next energy in order to begin search
      E2=E1+ABS(E1)/4.
      DE=2*ETOL
C                                                           
C     use secant search to find the bound state
50    IF (ABS(DE) .GT. ETOL) THEN
         CALL ACTION(E2,X1,X2,S)          !S at new energy
         F2=S-(N+.5)*PI                   !F at new energy
         IF (F1 .NE. F2) THEN             !calculate new DE
            DE=-F2*(E2-E1)/(F2-F1)          
         ELSE 
            DE=0.
         END IF
C
         E1=E2                             !roll values
         F1=F2
         E2=E1+DE                          !increment energy
         IF (E2 .GE. 0) E2=-ETOL           !keep energy negative
       GOTO 50
       END IF
C
       RETURN     
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE ACTION(E,X1,X2,S)
C calculates the (action integral)/2 (S) and the classical turning 
C points (X1,X2) for a given energy (E)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E1'
C Input/Output variables:     
      REAL E                   !energy (input)
      REAL S                   !action (output)
      REAL X1,X2               !turning points (output)
C Local variables:   
      REAL DX                  !increment in turning point search
      REAL H                   !quadrature step size
      REAL SUM                 !sum for integral
      INTEGER IFAC             !coefficient for Simpson's rule
      INTEGER IX               !index on X
      REAL X                   !current X value in sum
      REAL POT                 !potential as a function of X (function)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     find inner turning point; begin search at the well bottom
      X1=POTMIN  
      DX=.1
50    IF (DX .GT. XTOL) THEN
         X1=X1-DX              !use simple search, going inward
         IF (POT(X1) .GE. E) THEN
            X1=X1+DX
            DX=DX/2
         END IF
      GOTO 50                                   
      END IF
C
C     find the outer turning point; begin search at the well bottom
      X2=POTMIN
      DX=.1
120   IF (DX .GT. XTOL) THEN
         X2=X2+DX              !use simple search going outward
         IF (POT(X2) .GE. E) THEN
            X2=X2-DX    
            DX=DX/2
         END IF
      GOTO 120
      END IF
C 
C     Simpson's rule from X1+H to X2-H 
      IF (MOD(NPTS,2) .EQ. 1) NPTS=NPTS+1    !NPTS must be even
      H=(X2-X1)/NPTS                         !step size
      SUM=SQRT(E-POT(X1+H))                  
      IFAC=2
      DO 200 IX=2,NPTS-2
          X=X1+IX*H
          IF (IFAC .EQ. 2) THEN              !alternate factors
              IFAC=4
          ELSE
              IFAC=2
          END IF
          SUM=SUM+IFAC*SQRT(E-POT(X))
200   CONTINUE
      SUM=SUM+SQRT(E-POT(X2-H))
      SUM=SUM*H/3
C
C     special handling for sqrt behavior of first and last intervals
      SUM=SUM+SQRT(E-POT(X1+H))*2*H/3
      SUM=SUM+SQRT(E-POT(X2-H))*2*H/3
      S=SUM*GAMMA
C
      RETURN
      END     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       REAL FUNCTION POT(X)
C evaluates the Lennard-Jones potential at X
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
       REAL X                     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C If you change the potential, normalize to a minimum of -1
C and change the value of POTMIN in subroutine INIT to the
C new equilibrium position (i.e. the X value at which the force is zero)
C      Lennard-Jones potential in scaled variables
       POT=4*(X**(-12)-X**(-6))
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INIT
C initializes constants, displays header screen,
C initializes menu arrays for input parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'MENU.ALL'
      INCLUDE 'PARAM.E1'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)              
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get environment parameters
      CALL SETUP                
C 
C     display header screen     
      DESCRP(1)= 'EXAMPLE 1'
      DESCRP(2)= 'Bohr-Sommerfeld quantization for bound state'
      DESCRP(3)= 'energies of the 6-12 potential'
      NHEAD=3
C 
C     text output description
      DESCRP(4)= 'energy and classical turning points for each state'
      NTEXT=1
C 
C     graphics output description
      DESCRP(5)= 'phase space (wavenumber vs. position) portrait'
      DESCRP(6)= 'of classical trajectories'
      NGRAPH=2
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C 
C     calculate constants                
      PI=4*ATAN(1.0)
      POTMIN=2**(1.0/6.0)
C 
C     setup menu arrays, beginning with constant part
      CALL MENU                      
C                
      MTYPE(13)=FLOAT
      MPRMPT(13)= 'Enter gamma=sqrt(2*m*a**2*V/hbar**2) (dimensionless)'
      MTAG(13)= 'Gamma (dimensionless)'
      MLOLIM(13)=1.
      MHILIM(13)=500.
      MREALS(13)=50.
C                
      MTYPE(14)=SKIP
      MREALS(14)=35.
C 
      MTYPE(38)=FLOAT
      MPRMPT(38)= 'Enter tolerance for energy search (scaled units)'       
      MTAG(38)= 'Energy search tolerance (scaled units)'
      MLOLIM(38)=.00001
      MHILIM(38)=.01
      MREALS(38)=.0005
C 
      MTYPE(39)=FLOAT
      MPRMPT(39)= 
     +'Enter tolerance for turning point search (scaled units)'
      MTAG(39)= 'Turning point search tolerance (scaled units)'
      MLOLIM(39)=.00001
      MHILIM(39)=.01
      MREALS(39)=.0005
C             
      MTYPE(40)=NUM
      MPRMPT(40)= 'Enter number of points for action integral'
      MTAG(40)= 'Number of quadrature points for action integral'
      MLOLIM(40)=20.
      MHILIM(40)=5000.
      MINTS(40)=100
C 
      MTYPE(41)=SKIP
      MREALS(41)=60.
C 
      MSTRNG(MINTS(75))= 'exmpl1.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C 
      MSTRNG(MINTS(86))= 'exmpl1.grf'
C                     
      MTYPE(87)=NUM
      MPRMPT(87)= 'Enter number of points to be used in graphing'
      MTAG(87)= 'Number of graphing points'
      MLOLIM(87)= 10.
      MHILIM(87)= MAXGRF-2
      MINTS(87)= 80
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
C performs checks on parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'                                            
      INCLUDE 'PARAM.E1'
      INCLUDE 'MAP.E1'
C Local variables:
      REAL S                  !current value of action
      REAL E1                 !current value of energy
      REAL X1,X2              !current turning points
C Function:
      LOGICAL LOGCVT          !converts 1 and 0 to true and false 
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
      GAMMA=MREALS(IGAMMA)
      ETOL=MREALS(IETOL)
      XTOL=MREALS(IXTOL)
      NPTS=MINTS(INPTS)
C
C     text output parameters
      TTERM=LOGCVT(MINTS(ITTERM))
      TFILE=LOGCVT(MINTS(ITFILE))
      TNAME=MSTRNG(MINTS(ITNAME))
C
C     graphics output parameters
      GTERM=LOGCVT(MINTS(IGTERM))
      GHRDCP=LOGCVT(MINTS(IGHRD))
      GFILE=LOGCVT(MINTS(IGFILE))
      GNAME=MSTRNG(MINTS(IGNAME))
      NGRF=MINTS(INGRF)
C 
C     open files 
      IF (TFILE) CALL FLOPEN(TNAME,TUNIT)
      IF (GFILE) CALL FLOPEN(GNAME,GUNIT)
C     files may have been renamed
      MSTRNG(MINTS(ITNAME))=TNAME
      MSTRNG(MINTS(IGNAME))=GNAME
C
C     calculate total number of levels
      E1=-ETOL
      CALL ACTION(E1,X1,X2,S)
      NLEVEL=INT(S/PI-.5)+1
C     check value of GAMMA
      CALL PCHECK        
C                                                                         
      CALL CLEAR
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE PCHECK
C ensure that the number of states is not greater than the size of
C the data arrays; if so prompt for smaller GAMMA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global parameters:
      INCLUDE 'PARAM.E1'
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'
      INCLUDE 'MAP.E1'
C Local parameters:
      REAL S                   !action
      REAL E                   !small negative energy
      REAL X1,X2               !classical turning points
C Function:
      REAL GETFLT              !returns a floating point variable 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
10    IF ((NLEVEL-1) .GT. MAXLVL) THEN
          WRITE (OUNIT,15) NLEVEL,MAXLVL
          MHILIM(IGAMMA)=GAMMA            
          MREALS(IGAMMA) =  GETFLT(MREALS(IGAMMA)/2,MLOLIM(IGAMMA),
     +                  MHILIM(IGAMMA), 'Enter a smaller gamma') 
          GAMMA=MREALS(IGAMMA)
C 
          E=-ETOL
          CALL ACTION(E,X1,X2,S)
          NLEVEL=INT(S/PI+.5)+1
      GOTO 10
      END IF
C 
15    FORMAT (' Total number of levels (=',i5, 
     +        ') is larger than maximum allowable (=',i3,')')
C
      RETURN  
      END 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE PRMOUT(MUNIT,NLINES)
C outputs parameter summary to the specified unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'IO.ALL'
       INCLUDE 'PARAM.E1'
C Input/Output variables:               
       INTEGER MUNIT            !unit number for output (input)
       INTEGER NLINES           !number of lines written so far (output)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (MUNIT .EQ. OUNIT) CALL CLEAR
C 
       WRITE (MUNIT,2)
       WRITE (MUNIT,4)
       WRITE (MUNIT,6) ETOL,XTOL
       WRITE (MUNIT,8) NPTS
       WRITE (MUNIT,10) GAMMA,NLEVEL
       WRITE (MUNIT,12)                                                 
       WRITE (MUNIT,2)
C                                          
       IF (MUNIT .NE. GUNIT) THEN
         WRITE (MUNIT,20)
         WRITE (MUNIT,25)
       END IF
C 
       NLINES=7
C
2      FORMAT (' ')
4      FORMAT (' Output from example 1: Bohr Sommerfeld Quantization')
6      FORMAT (' Energy tolerance =', E12.5,
     +         '   position tolerance =', E12.5)
8      FORMAT (' number of quadrature points =', I4)
10     FORMAT (' For gamma =', F8.2,' there are ', I4, ' levels:')
12     FORMAT (' (all quantities are expressed in scaled units)')
20     FORMAT (8X, 'Level', 8X, 'Energy', 12X, 'Xmin', 12X, 'Xmax')
25     FORMAT (8X, '-----', 8X, '------', 12X, '----', 12X ,'----')
C 
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,ILEVEL,E,X1,X2,NLINES)       
C writes results for one state to the requested unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:
      INTEGER MUNIT                 !output unit specifier
      INTEGER ILEVEL                !current level
      REAL E                        !eigen energy
      REAL X1,X2                    !classical turning points
      INTEGER NLINES                !number of lines printed so far
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     if screen is full, clear screen and retype headings          
      IF ((MOD(NLINES,TRMLIN-6) .EQ. 0) 
     +                          .AND. (MUNIT .EQ. OUNIT)) THEN
         CALL PAUSE('to continue...',1)
         CALL CLEAR
         WRITE (MUNIT,20)
         WRITE (MUNIT,25)
      END IF
C
      WRITE (MUNIT,30) ILEVEL,E,X1,X2
C 
C     keep track of printed lines only for terminal output
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+1
C
20    FORMAT (8X, 'Level', 8X, 'Energy', 12X, 'Xmin', 12X, 'Xmax')
25    FORMAT (8X, '-----', 8X, '------', 12X, '----', 12X ,'----')
30    FORMAT(8X,I4,3(8X,F8.5))
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE GRFOUT(DEVICE,ENERGY,XIN,XOUT)
C outputs phase space portraits of the bound states to the terminal
C and/or a file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global parameters:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E1'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE               !which device is being used?
      REAL ENERGY(0:MAXLVL)        !energy of bound state
      REAL XIN(0:MAXLVL)           !inner turning point
      REAL XOUT(0:MAXLVL)          !outer turning point
C Local parameters:                         
      INTEGER ILEVEL               !level index
      REAL H                       !step size for x
      INTEGER IX                   !x index 
      REAL E                       !current energy    
      REAL K(MAXGRF)               !current wavenumber
      REAL X(MAXGRF)               !current position
      CHARACTER*9 CGAMMA,CN        !Gamma,nlevel as a character string
      REAL POT                     !potential (function)
      INTEGER LEN,NLEN             !length of character data
      INTEGER SCREEN               !send to terminal
      INTEGER PAPER                !make a hardcopy
      INTEGER FILE                 !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     calculate parameters for graphing
      IF (DEVICE .NE. FILE) THEN
          NPLOT=1                       !how many plots
          IPLOT=1
C 
          YMAX=GAMMA*SQRT(1.+ENERGY(NLEVEL-1))   !limits on data points
          YMIN=-YMAX  
          XMIN=XIN(NLEVEL-1)
          XMAX=XOUT(NLEVEL-1)
          Y0VAL=XMIN
          X0VAL=0.
C 
          IF (MOD(NGRF,2) .EQ. 0) NGRF=NGRF+1
          NPOINT=NGRF                   !keep number of points odd
C                                                            
          ILINE=1                       !line and symbol styles
          ISYM=1
          IFREQ=0
          NXTICK=5
          NYTICK=5
C 
          CALL CONVRT(GAMMA,CGAMMA,LEN) !titles and labels
          CALL ICNVRT(NLEVEL,CN,NLEN)
          INFO=' NLEVEL = '//CN(1:NLEN)
          TITLE='Semiclassically Quantized Trajectories, Gamma='//CGAMMA
          LABEL(1)='scaled position'
          LABEL(2)='scaled wave number'
C 
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
          CALL LNLNAX                          !draw axes
      END IF                                   
C 
C     calculate classical phase space trajectory for each bound state 
C     by finding the scaled wavenumber as a function of X and Energy
      DO 50 ILEVEL=0,NLEVEL-1
         E=ENERGY(ILEVEL)
         H=(XOUT(ILEVEL)-XIN(ILEVEL))/((NGRF-1)/2)  !step size
         X(1)=XIN(ILEVEL)
         K(1)=0.
C 
         DO 20 IX=1,(NGRF-1)/2
            X(IX+1)=XIN(ILEVEL)+(IX)*H
            K(IX+1)=(E-POT(X(IX+1)))                !scaled wave number
            IF (K(IX) .LE. 0) THEN
               K(IX)=0.
            ELSE
               K(IX)=GAMMA*SQRT(K(IX))
            END IF                                                  
20       CONTINUE
C 
         DO 30 IX=(NGRF+1)/2,NGRF-1   !graph is symmetric about x-axis
            X(IX+1)=X(NGRF-IX)
            K(IX+1)=-K(NGRF-IX)
30       CONTINUE
C         
C        output results
         IF (DEVICE .EQ. FILE) THEN
             WRITE (GUNIT,75) E
             WRITE (GUNIT,70) (X(IX),K(IX),IX=1,NGRF)
         ELSE  
             CALL XYPLOT(X,K)
         END IF
50    CONTINUE
C 
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)   !close graphing package
      IF (DEVICE .EQ. SCREEN) CALL TMODE         !switch to text mode
C 
70    FORMAT (2(5X,E11.3))
75    FORMAT (/,' Position vs. wave number for Energy =',1PE11.3)
100   FORMAT (/,' Patience, please; output going to a file.')
C 
      RETURN
      END
                            
