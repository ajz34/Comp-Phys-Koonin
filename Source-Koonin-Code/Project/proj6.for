CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM PROJ6
C   Project 6: 2-D viscous incompressible flow about a rectangular block
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop/ execute once for each set of param
        CALL LATTIC       !get lattice size, allow for ending
        CALL ARCHON       !get detailed bound cond and relax lattice
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON
C subroutine to get parameters, relax the lattice, and output results
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P6'
      INCLUDE 'IO.ALL'
C Local variables:
      REAL P(NXMAX,NYMAX)            !stream function
      REAL XSI(NXMAX,NYMAX)          !vorticity
      REAL VFORCE,PFORCE             !viscous force and pressure force
      REAL XSIMAX,PMAX,XSIMIN,PMIN   !min,max values of XSI, P
      INTEGER XMIN,XMAX,YMIN,YMAX    !corners of relaxing lattice
      LOGICAL ITER                   !continue iterating?
      LOGICAL PRM                    !print out parameters?
      INTEGER CHOICE                 !write out field now? 
      INTEGER NITER                  !number of iterations
      LOGICAL END                    !end this run?   
      INTEGER IX,IY                  !lattice indices
      INTEGER NLINES                 !number of lines sent to terminal
      INTEGER SCREEN                 !send to terminal
      INTEGER PAPER                  !make a hardcopy
      INTEGER FILE                   !send to a file
      REAL FGRF(NXMAX,NYMAX)         !field for graphing
C Functions:
      REAL GETFLT
      INTEGER GETINT,YESNO                          
      LOGICAL LOGCVT
      DATA SCREEN,PAPER,FILE/1,2,3/
      DATA VRTCTY,STREAM /1,2/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END=.FALSE.                 !initialize values
      DO 5 IX=1,NX
       DO 5 IY=1,NY
          P(IX,IY)=0.
          XSI(IX,IY)=0.
5     CONTINUE
C            
200   CONTINUE                    !allow for many runs with same lattice
C
        CALL PARAM(P,XSI,END)     !get new parameters
        IF (END) RETURN           !start fresh or end altogether
C 
        IF (SUBLAT) THEN
          XMIN=NXLL               !set limits for relaxation
          XMAX=NXUR
          YMIN=NYLL
          YMAX=NYUR
        ELSE
          XMIN=1
          YMIN=1
          XMAX=NX
          YMAX=NY
        END IF
C
        IF (TFILE) CALL PRMOUT(TUNIT,NLINES)
        NITER=0
        PRM=.TRUE.
C
99      CONTINUE                   !begin iterations
          IF ((TTERM) .AND. (PRM)) CALL PRMOUT(OUNIT,NLINES)
          NITER=NITER+1
C
C         relax stream function, then vorticity; calculate forces
          CALL PRELAX(P,XSI,XMIN,XMAX,YMIN,YMAX,PMAX,PMIN)
          CALL VRELAX(P,XSI,XMIN,XMAX,YMIN,YMAX,XSIMAX,XSIMIN)
          CALL FORCES(XSI,VFORCE,PFORCE)
C          
          IF (TFILE) CALL TXTOUT(TUNIT,NITER,
     +               VFORCE,PFORCE,XSIMAX,PMAX,XSIMIN,PMIN,NLINES)
          IF (TTERM) CALL TXTOUT(OUNIT,NITER,
     +               VFORCE,PFORCE,XSIMAX,PMAX,XSIMIN,PMIN,NLINES)
C
          ITER=.FALSE.
          PRM=.TRUE.
          IF (MOD(NITER,NFREQ) .NE. 0) THEN         
             ITER=.TRUE.        ! continue iterating without asking
             PRM=.FALSE.        ! don't print out header
          ELSE                  ! otherwise display fields 
             IF (GTERM) THEN
                CALL PAUSE('to see stream function and vorticity...',1)
                CALL CLEAR                                             
                CALL GRFOUT(SCREEN,STREAM,P,PMIN,PMAX,FGRF,NX,NY)
                CALL GRFOUT(SCREEN,VRTCTY,XSI,XSIMIN,XSIMAX,FGRF,NX,NY)
             ELSE IF (TTERM) THEN
                CALL PAUSE('to see stream function...',1)
                CALL CLEAR                                             
                CALL DISPLY(P,OUNIT,STREAM,PMAX,PMIN)
                CALL PAUSE('to see vorticity...',0)
                CALL CLEAR
                CALL DISPLY(XSI,OUNIT,VRTCTY,XSIMAX,XSIMIN)
             END IF
             IF (SUBLAT) THEN   !and provide options for continuing
               SUBLAT=LOGCVT(YESNO(1,'Continue relaxing sublattice?'))
               IF (SUBLAT) ITER=.TRUE.
             ELSE
               ITER=LOGCVT(YESNO(1,'Continue relaxing lattice?'))
             END IF
          END IF
C
       IF (ITER) THEN       !continue iterating
          GOTO 99
C
       ELSE
C       prompt for writing out of field values
        IF ((GFILE) .OR. (GHRDCP) .OR. (TFILE)) THEN
          CHOICE=YESNO(1,'Do you want to write out field now?')
          IF (CHOICE .EQ. 1) THEN
           IF (GFILE) CALL WRTOUT(P,XSI)   !write out field if requested
           IF (GHRDCP) THEN
              CALL GRFOUT(PAPER,STREAM,P,PMIN,PMAX,FGRF,NX,NY)
              CALL GRFOUT(PAPER,VRTCTY,XSI,XSIMIN,XSIMAX,FGRF,NX,NY)
           END IF
           IF (TFILE) CALL DISPLY(P,TUNIT,STREAM,PMAX,PMIN)
           IF (TFILE) CALL DISPLY(XSI,TUNIT,VRTCTY,XSIMAX,XSIMIN)
          END IF
        END IF
        GOTO 200       !display menu
       END IF
C     
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FORCES(XSI,VFORCE,PFORCE)
C calculate pressure force (PFORCE) and viscous forces (VFORCE)
C from the vorticity (XSI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global Variables:
      INCLUDE 'PARAM.P6'
C Passed Variables:
      REAL XSI(NXMAX,NYMAX)    !vorticity (input)
      REAL VFORCE,PFORCE       !viscous force and pressure force (output)
C Local variables:
      INTEGER IX,IY            !lattice indices
      REAL P                   !pressure
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     calculate viscous forces by integrating vorticity 
C     along the top of the plate
      VFORCE=XSI(FRONT,HFWID)/2
      DO 10 IX=FRONT+1,BACK
         VFORCE=VFORCE+XSI(IX,HFWID)
10    CONTINUE
      VFORCE=VFORCE-XSI(BACK,HFWID)/2
      VFORCE=VFORCE/REYNLD/(HFWID-1)
C
C     calculate the pressure force by integrating vorticity along the
C     edges of the plate to obtain the pressure, then integrating
C     the pressure to obtain the pressure force
      P=0
      PFORCE=0
      DO 20 IY=2,HFWID         !up front
         P=P-(XSI(FRONT,IY)-XSI(FRONT-1,IY)
     +        +XSI(FRONT,IY-1)-XSI(FRONT-1,IY-1))/2/REYNLD
         PFORCE=PFORCE+P
20    CONTINUE
      PFORCE=PFORCE-P/2
C
      DO 30 IX=FRONT+1,BACK    !across top
         P=P+(XSI(IX,HFWID+1)-XSI(IX,HFWID)  
     +        +XSI(IX-1,HFWID+1)-XSI(IX-1,HFWID))/2/REYNLD
30    CONTINUE
      PFORCE=PFORCE-P/2
C
      DO 40 IY=HFWID-1,1,-1    !down back
         P=P+(XSI(BACK+1,IY)-XSI(BACK,IY)
     +        +XSI(BACK+1,IY+1)-XSI(BACK,IY+1))/2/REYNLD
         PFORCE=PFORCE-P
40    CONTINUE
      PFORCE=PFORCE/(HFWID-1)
C
      RETURN
      END 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRELAX(P,XSI,XMIN,XMAX,YMIN,YMAX,PMAX,PMIN)
C relaxes the stream function (P) given previous P and vorticity (XSI)
C on a lattice with corners at XMIN,XMAX,YMIN,YMAX; 
C returns min, max value of P in PMIN, PMAX
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P6'
C Passed variables:
      REAL P(NXMAX,NYMAX)            !stream function (I/O)
      REAL XSI(NXMAX,NYMAX)          !vorticity (output)
      INTEGER XMIN,XMAX,YMIN,YMAX    !edges of relaxing lattice (input)
      REAL PMAX,PMIN                 !min,max value of P (output)
C Local variables:                        
      INTEGER IX,IY                  !lattice indices
      REAL TEMP                      !temporary storage
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PMAX=0.
      PMIN=0.
      DO 100 IX=XMIN,XMAX
        DO 200 IY=YMIN,YMAX
          IF (BNDCND(IX,IY) .EQ. 0) THEN  !relax only non-bound points
C
            IF (IY .EQ. NY-1) THEN        !just below top edge
C
              IF (IX .EQ. 2) THEN         !left corner
                 TEMP=P(IX,IY-1)+P(IX+1,IY)-XSI(IX,IY)+1 
                 P(IX,IY)=SOMEGA/2*TEMP+MSOMEG*P(IX,IY)
                 P(IX-1,IY)=P(IX,IY)
                 P(IX-1,NY)=P(IX-1,IY)+1
C
              ELSE IF (IX .EQ. NX-1) THEN !right corner
                 TEMP=P(IX,IY-1)+P(IX-1,IY)-XSI(IX,IY)+1
                 P(IX,IY)=SOMEGA/2*TEMP+MSOMEG*P(IX,IY)
                 P(NX,IY)=P(IX,IY)
                 P(NX,NY)=P(NX,IY)+1                                    
C
               ELSE                       !not a corner
                 TEMP=P(IX,IY-1)+P(IX-1,IY)+P(IX+1,IY)-XSI(IX,IY)+1.
                 P(IX,IY)=SOMEGA*TEMP/3+MSOMEG*P(IX,IY)
               END IF
C
               P(IX,NY)=P(IX,IY)+1        !top edge given by bound cond
C
           ELSE                           !not just below top edge
C
              IF (IX .EQ. 2) THEN         !front edge
                 TEMP=P(IX,IY-1)+P(IX,IY+1)+P(IX+1,IY)-XSI(IX,IY)
                 P(IX,IY)=SOMEGA*TEMP/3+MSOMEG*P(IX,IY)
                 P(IX-1,IY)=P(IX,IY)
C
              ELSE IF (IX .EQ. NX-1) THEN !back edge
                 TEMP=P(IX,IY-1)+P(IX-1,IY)+P(IX,IY+1)-XSI(IX,IY)
                 P(IX,IY)=SOMEGA/3*TEMP+MSOMEG*P(IX,IY)
                 P(NX,IY)=P(IX,IY)
C
               ELSE                       !interior point
                 TEMP=P(IX,IY-1)+P(IX,IY+1)+P(IX-1,IY)+P(IX+1,IY)
     +                -XSI(IX,IY)
                 P(IX,IY)=SOMEGA*TEMP/4+MSOMEG*P(IX,IY)
               END IF
C
            END IF
            IF (P(IX,IY) .GT. PMAX) PMAX=P(IX,IY)
            IF (P(IX,NY) .GT. PMAX) PMAX=P(IX,NY)
            IF (P(IX,IY) .LT. PMIN) PMIN=P(IX,IY)
            IF (P(IX,NY) .LT. PMIN) PMIN=P(IX,NY)
C
          END IF
200     CONTINUE
100   CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE VRELAX(P,XSI,XMIN,XMAX,YMIN,YMAX,XSIMAX,XSIMIN)
C relaxes the vorticity (XSI) given previous XSI and stream function (P)
C on a lattice defined by XMIN,XMAX,YMIN,YMAX            
C returns min, max value of XSI in XSIMIN, XSIMAX
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P6'
C Passed variables:
      REAL P(NXMAX,NYMAX)            !stream function (input)
      REAL XSI(NXMAX,NYMAX)          !vorticity (I/O)
      INTEGER XMIN,XMAX,YMIN,YMAX    !edges of relaxing lattice (input)
      REAL XSIMAX,XSIMIN             !max and min value of XSI (output)         
C Local variables:
      INTEGER IX,IY                  !lattice indices
      REAL TEMP,TEMP2,TEMP3          !temporary storage
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      XSIMAX=0.
      XSIMIN=0.
C     impose Dirichlet boundary conditons along the plate
      DO 10 IY=1,HFWID           
         XSI(FRONT,IY)=2*P(FRONT-1,IY)
         XSI(BACK,IY)=2*P(BACK+1,IY)
         IF (XSI(FRONT,IY) .GT. XSIMAX) XSIMAX=XSI(FRONT,IY)
         IF (XSI(BACK,IY) .GT. XSIMAX) XSIMAX=XSI(BACK,IY)
         IF (XSI(FRONT,IY) .LT. XSIMIN) XSIMIN=XSI(FRONT,IY)
         IF (XSI(BACK,IY) .LT. XSIMIN) XSIMIN=XSI(BACK,IY)
10    CONTINUE
      DO 20 IX=FRONT+1,BACK-1
         XSI(IX,HFWID)=2*P(IX,HFWID+1)
         IF (XSI(IX,HFWID) .GT. XSIMAX) XSIMAX=XSI(IX,HFWID)
         IF (XSI(IX,HFWID) .LT. XSIMIN) XSIMIN=XSI(IX,HFWID)
20    CONTINUE
C
      DO 30 IX=XMIN,XMAX
        DO 40 IY=YMIN,YMAX
         IF (BNDCND(IX,IY) .EQ. 0) THEN   !don't relax bound points
C            
          IF (IX .EQ. NX-1) THEN          !right edge, Neumann BC
           TEMP=XSI(IX,IY+1)+XSI(IX,IY-1)+XSI(IX-1,IY)
           TEMP2=REYND4*(P(IX,IY+1)-P(IX,IY-1))*XSI(IX-1,IY)
           TEMP3=REYND4*(P(IX+1,IY)-P(IX-1,IY))
     +           *(XSI(IX,IY+1)-XSI(IX,IY-1))
           XSI(IX,IY)=MVOMEG*XSI(IX,IY)+
     +      VOMEGA*(TEMP+TEMP2+TEMP3)/(3+REYND4*(P(IX,IY+1)-P(IX,IY-1)))
           XSI(NX,IY)=XSI(IX,IY)
C
          ELSE                             !interior point
           TEMP=XSI(IX,IY+1)+XSI(IX,IY-1)+XSI(IX+1,IY)+XSI(IX-1,IY)
           TEMP2=-REYND4*(P(IX,IY+1)-P(IX,IY-1))
     +           *(XSI(IX+1,IY)-XSI(IX-1,IY))
           TEMP3=REYND4*(P(IX+1,IY)-P(IX-1,IY))
     +           *(XSI(IX,IY+1)-XSI(IX,IY-1))
           XSI(IX,IY)=VOMEGA/4*(TEMP+TEMP2+TEMP3)+MVOMEG*XSI(IX,IY)
C
          END IF
          IF (XSI(IX,IY) .GT. XSIMAX) XSIMAX=XSI(IX,IY)
          IF (XSI(IX,IY) .LT. XSIMIN) XSIMIN=XSI(IX,IY)
C
         END IF
40      CONTINUE
30    CONTINUE
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
      INCLUDE 'PARAM.P6'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)         
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get environment parameters
      CALL SETUP                
C 
C     display header screen     
      DESCRP(1)= 'PROJECT 6'
      DESCRP(2)= '2-D viscous incompressible flow about a '
     +           //'rectangular block'
      NHEAD=2
C      
C     text output description
      DESCRP(3)= 'maximum vorticity and stream function,'
      DESCRP(4)= 'forces at each iteration'
      DESCRP(5)= '(all values are in scaled units)'
      NTEXT=3
C 
C     graphics output description
      DESCRP(6)= 'stream function and vorticity'
      NGRAPH=1
C
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C
      FIRST=.TRUE.             !is this the first time through the menu?
C
C     set up constant part of menu
      CALL MENU
C
C     item 9 has a different use in this program
      MTYPE(9)=MTITLE
      MPRMPT(9)='7) Change lattice size or end altogether'
      MLOLIM(9)=0
      MHILIM(9)=1
C
C     default values and limits of next 3 items are set in subroutine 
C     LATTIC since they depend on the size of lattice chosen
      MTYPE(13)=NUM
      MPRMPT(13)= 'Enter half-width of plate in lattice spacings'
      MTAG(13)= 'Half-width of plate in lattice spacings'
C                
      MTYPE(14)=NUM
      MPRMPT(14)= 'Enter length of plate in lattice spacings'
      MTAG(14)= 'Length of plate in lattice spacings'
C                
      MTYPE(15)=NUM
      MPRMPT(15)= 
     + 'Enter location of plate''s front edge in lattice spacings'
      MTAG(15)= 'Location of plate''s front edge in lattice spacings'
C                
      MTYPE(16)=FLOAT
      MPRMPT(16)= 'Enter lattice Reynold''s number'
      MTAG(16)= 'Lattice Reynold''s number'
      MLOLIM(16)=0.01
      MHILIM(16)=20.
      MREALS(16)=1.                               
C                
      MTYPE(17)=SKIP
      MREALS(17)=35.
C
      MTYPE(38)=FLOAT
      MPRMPT(38)= 'Enter vorticity relaxation parameter'
      MTAG(38)= 'Vorticity relaxation parameter'        
      MLOLIM(38)=0.
      MHILIM(38)=1.9
      MREALS(38)=.3
C            
      MTYPE(39)=FLOAT
      MPRMPT(39)= 'Enter stream function relaxation parameter'
      MTAG(39)= 'Stream function  relaxation parameter'
      MLOLIM(39)=0.
      MHILIM(39)=1.9                      
      MREALS(39)=.3
C            
      MTYPE(40)=NOSKIP
      MPRMPT(40)='Do you want to specify a sublattice to relax first?'
      MTAG(40)='Relax a sublattice first?'
      MINTS(40)=0
      MREALS(40)=45.
C                     
      MTYPE(41)=NUM
      MPRMPT(41)=
     + 'Enter lower left X (in lattice spacings) for sublattice'
      MTAG(41)='Lower left X sublattice value'
      MLOLIM(41)=1
      MINTS(41)=1
C
      MTYPE(42)=NUM
      MPRMPT(42)=
     + 'Enter lower left Y (in lattice spacings) for sublattice'
      MTAG(42)='Lower left Y sublattice value'
      MLOLIM(42)=1
      MINTS(42)=1
C               
      MTYPE(43)=NUM
      MPRMPT(43)='Enter upper right X value for sublattice'
      MTAG(43)='Upper right X sublattice value'
      MLOLIM(43)=1
      MINTS(43)=1
C               
      MTYPE(44)=NUM
      MPRMPT(44)='Enter upper right Y value for sublattice'
      MTAG(44)='Upper right Y sublattice value'
      MLOLIM(44)=1
      MINTS(44)=1
C
      MTYPE(45)=MTITLE
      MPRMPT(45)='Starting Values (for stream function and Vorticity)'
      MLOLIM(45)=2
      MHILIM(45)=1
C
      MTYPE(46)=MTITLE
      MPRMPT(46)='1) Read in starting values from a file'
      MLOLIM(46)=0
      MHILIM(46)=0
C                
      MTYPE(47)=MTITLE
      MPRMPT(47)='2) Set values to free flow values'
      MLOLIM(47)=0            
      MHILIM(47)=0
C 
      MTYPE(48)=MTITLE
      MPRMPT(48)='3) Leave field values unchanged'
      MLOLIM(48)=0
      MHILIM(48)=0
C                   
      MTYPE(49)=MCHOIC
      MPRMPT(49)='Enter Choice'
      MTAG(49)='50 51 51'
      MLOLIM(49)=1
      MHILIM(49)=3
      MINTS(49)=2
      MREALS(49)=2
C                        
      MTYPE(50)=CHSTR
      MPRMPT(50)= 'Enter name of data file'
      MTAG(50)= 'File with initial values for the fields'
      MHILIM(50)=12.
      MINTS(50)=3.
      MSTRNG(MINTS(50))= 'proj6.in'
C
      MTYPE(51)=SKIP
      MREALS(51)=60.
C                              
      MSTRNG(MINTS(75))= 'proj6.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C 
      MSTRNG(MINTS(86))= 'proj6.grf'
C                     
      MTYPE(87)=NUM
      MPRMPT(87)= 'Enter the display frequency for the fields'
      MTAG(87)= 'Field display frequency'
      MLOLIM(87)= 1.
      MHILIM(87)= 100.
      MINTS(87)= 10
C                   
      MTYPE(88)=NUM
      MPRMPT(88)= 'Enter number of contour levels'
      MTAG(88)= 'Number of contour levels'
      MLOLIM(88)= 1.
      MHILIM(88)= MAXLEV
      MINTS(88)= 10
C                   
      MTYPE(89)=SKIP
      MREALS(89)=90.
C
      RETURN
      END          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LATTIC
C gets lattice size from screen and calculates best way to display the
C field as ascii charcters based on lattice size and terminal size;
C resets all boundary conditions and default menu values
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P6'
      INCLUDE 'IO.ALL'             
      INCLUDE 'MENU.ALL'
      INCLUDE 'MAP.P6'
C Local variables:
      INTEGER END                   !end program
      INTEGER IX,IY,IBC             !lattice indices, BC index
      INTEGER NXHI,NYHI,NXLO,NYLO   !limits on lattice size
      INTEGER NXDEF,NYDEF           !default lattice sizes
      LOGICAL RESET                 !reset parameters?
C Functions:
      INTEGER YESNO,GETINT          !user input functions
      LOGICAL LOGCVT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     allow user to end the program
      CALL CLEAR
      IF (.NOT. FIRST) THEN
         END=YESNO(0,' Do you want to end the program?')
         IF (END .EQ. 1) CALL FINISH
      ELSE
C        the lattice size is determined by array size and terminal size;
C        if you're using graphics, terminal size won't matter
C        set NXHI=NXMAX and NYHI=NYMAX
         NXHI=MIN(TRMWID-2,NXMAX)          
         NYHI=MIN(TRMLIN-3,NYMAX)
         NXDEF=MIN(NXHI,20)
         NYDEF=MIN(NYHI,20)
         NXLO=MIN(5,NXDEF)
         NYLO=MIN(5,NYDEF)
      END IF
C
C     get lattice parameters from the terminal
      NX=GETINT(NXDEF,NXLO,NXHI,' Enter number of X lattice points')
      NY=GETINT(NYDEF,NYLO,NYHI,' Enter number of Y lattice points')
      NXDEF=NX
      NYDEF=NY
C
C     calculate parameters for best looking display
      IF (2*NX .LE. TRMWID) THEN
           XSKIP=.TRUE.              !skip spaces in x
           XCNTR=(TRMWID-2*NX)/2     !how to center display
      ELSE
           XSKIP=.FALSE.
           XCNTR=(TRMWID-NX)/2
      END IF
      IF (XCNTR .LT. 1) XCNTR=1
C
      IF (2*NY .LE. TRMLIN-3) THEN
          YSKIP=.TRUE.               !skip lines in y
          YCNTR=(TRMLIN-2*NY)/2-2    !how to center display
      ELSE
          YSKIP=.FALSE.
          YCNTR=(TRMLIN-NY)/2-2
      END IF
      IF (YCNTR .LT. 0) YCNTR=0
C
C     set up default and limits for the plate location and geometry
      MHILIM(IWID)=NY
      MLOLIM(IWID)=2
      MINTS(IWID)=MAX(2,NY/3)
      MHILIM(ILEN)=NX
      MLOLIM(ILEN)=0
      MINTS(ILEN)=NX/4
      MHILIM(IFRNT)=NX
      MLOLIM(IFRNT)=1
      MINTS(IFRNT)=MAX(1,NX/3)
C
C     set up limits on sublattice
      MHILIM(INXLL)=NX          
      MHILIM(INYLL)=NY
      MHILIM(INXUR)=NX
      MHILIM(INYUR)=NY
      MINTS(41)=1                    !sublattice  x, lower left
      MINTS(42)=1                    !sublattice  y, lower left
      MINTS(43)=1                    !sublattice  x, upper right
      MINTS(44)=1                    !sublattice  y, lower left
C                                                  
C     allow for resetting of defaults
      IF (FIRST) THEN
         FIRST=.FALSE.
      ELSE
         RESET=LOGCVT(YESNO(0,' Do you want to reset default values?'))
         IF (RESET) THEN
            MREALS(16)=1.                  !Reynolds number
            MREALS(38)=0.3                 !vorticity relaxation
            MREALS(39)=0.3                 !stream relaxation
            MINTS(40)=0                    !no sublattice
            MSTRNG(MINTS(50))= 'proj6.in'  !input file
            MSTRNG(MINTS(75))= 'proj6.txt' !text file
            MSTRNG(MINTS(86))= 'proj6.grf' !graphics file
            MINTS(87)= 10                  !graphing frequency
            MINTS(88)= 10                  !number of contours
            MINTS(49)=2                    !start from free 
            MREALS(49)=2                   !   streaming values
            MINTS(73)=TXTTRM               !default i/o
            MINTS(74)=TXTFIL
            MINTS(83)=GRFTRM
            MINTS(84)=GRFHRD
            MINTS(85)=GRFFIL
         END IF
      END IF
C
      RETURN 
      END                                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PARAM(P,XSI,END)
C gets parameters from screen                                      
C ends program on request
C closes old files
C maps menu variables to program variables              
C opens new files
C calculates all derivative parameters
C performs checks on sublattice parameters
C set the field to its initial values (P,XSI)
C and controls ending of program (END)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'                                            
      INCLUDE 'PARAM.P6'
      INCLUDE 'MAP.P6'
C Input/output variables:
      REAL P(NXMAX,NYMAX)     !stream function
      REAL XSI(NXMAX,NYMAX)   !vorticity
      LOGICAL END             !end program?
C Local variables:
      INTEGER IX,IY           !lattice indices
C Functions:
      LOGICAL LOGCVT          !converts 1 and 0 to true and false 
      INTEGER GETINT          !get integer from screen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get input from terminal
      CALL CLEAR
      CALL ASK(1,ISTOP)
C
C     start fresh or end altoghether, if so requested
      IF (MREALS(IMAIN) .EQ. STOP)  THEN
         END=.TRUE.
         RETURN
      END IF
C
C     close files if necessary
      IF (TNAME .NE. MSTRNG(MINTS(ITNAME))) 
     +     CALL FLCLOS(TNAME,TUNIT)
C                                           
C     physical and numerical parameters
      HFWID=MINTS(IWID)
      LENGTH=MINTS(ILEN)
      FRONT=MINTS(IFRNT)
      REYNLD=MREALS(IRNLDS)
      VOMEGA=MREALS(IVOMEG)
      SOMEGA=MREALS(ISOMEG)
      SUBLAT=LOGCVT(MINTS(ISUB))
      NXLL=MINTS(INXLL)
      NYLL=MINTS(INYLL)
      NXUR=MINTS(INXUR)
      NYUR=MINTS(INYUR)
      PTYPE=ABS(MREALS(IPTYPE))
      PFILE=MSTRNG(MINTS(IPFILE))
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
      NFREQ=MINTS(INFREQ)
      NLEV=MINTS(INLEV)     
C             
C     open files 
      IF (TFILE) CALL FLOPEN(TNAME,TUNIT)
      !files may have been renamed
      MSTRNG(MINTS(ITNAME))=TNAME
C
C     check sublattice parameters
39    IF ((NXLL .GT. NXUR) .OR. (NYLL .GT. NYUR)) THEN
         WRITE (OUNIT,40)
         CALL ASK(41,44)
         NXLL=MINTS(INXLL)
         NYLL=MINTS(INYLL)
         NXUR=MINTS(INXUR)
         NYUR=MINTS(INYUR)
         GOTO 39
      END IF
40    FORMAT (' Sublattice parameters must have upper right values '
     +        'greater than lower left')
      CALL CLEAR
C 
C     calculate derivative quantities
      BACK=FRONT+LENGTH-1
      IF (BACK .GT. NX) BACK=NX  !make sure block doesn't extend too far
      MVOMEG=1.-VOMEGA
      MSOMEG=1.-SOMEGA
      REYND4=REYNLD/4
C
      CALL INTCND(P,XSI)         !get starting values
C
C     reset PTYPE menu parameters so that P remains unchanged during
C     runs unless the user explicitly requests otherwise 
      MINTS(IPTYPE)=3
      MREALS(IPTYPE)=3              
C                                   
      RETURN                                                           
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTCND(P,XSI)
C get starting values for stream function (P) and vorticity (XSI)
C and set boundary conditions in BNDCND, depending on plate location
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P6'
C Output variables:
      REAL P(NXMAX,NYMAX)           !stream function
      REAL XSI(NXMAX,NYMAX)         !vorticity
C Local variables:
      INTEGER IX,IY                 !lattice indices
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     set P to it's initial value
      IF (PTYPE .EQ. 1) CALL READIN(P,XSI)         !read data in
C     sometimes in READIN, PTYPE is changed
      IF (PTYPE .EQ. 2) THEN
        DO 60 IX=1,NX
           DO 70 IY=1,NY               !set to free stream values
              P(IX,IY)=REAL(IY-1)    
              XSI(IX,IY)=0.
70         CONTINUE
60      CONTINUE
      ELSE IF (PTYPE .EQ. 3) THEN      !keep P and XSI the same
        CONTINUE      
      END IF 
C
C     release previous boundary conditions
      DO 61 IX=1,NX
         DO 71 IY=1,NY             
            BNDCND(IX,IY)=0.
71       CONTINUE
61    CONTINUE
C
C     set up boundary conditions on plate and lattice edges
      DO 100 IY=1,NY               !left and right edges are boundaries
         BNDCND(1,IY)=1
         BNDCND(NX,IY)=1
100   CONTINUE
      DO 110 IX=1,NX               !upper and lower edges are boundaries
         BNDCND(IX,1)=1       
         BNDCND(IX,NY)=1
110   CONTINUE
      DO 120 IX=FRONT,BACK         !top of plate is boundary
         P(IX,HFWID)=0
         BNDCND(IX,HFWID)=1
120   CONTINUE
      DO 130 IY=1,HFWID            !front and back of plate are bound
         P(FRONT,IY)=0.  
         P(BACK,IY)=0.
         BNDCND(FRONT,IY)=1
         BNDCND(BACK,IY)=1
130   CONTINUE
      DO 80 IX=FRONT+1,BACK-1      !plate interior is a boundary
         DO 90 IY=1,HFWID-1
            P(IX,IY)=0.
            XSI(IX,IY)=0.
            BNDCND(IX,IY)=2        !2 indicates actual plate
90       CONTINUE
80    CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE READIN(P,XSI)
C read in stream function P and vorticity XSI from file written
C by subroutine WRTOUT; files must have same lattice parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P6'
      INCLUDE 'MENU.ALL'
      INCLUDE 'MAP.P6'
C Output variables:
      REAL P(NXMAX,NYMAX)           !stream function
      REAL XSI(NXMAX,NYMAX)         !vorticity
C Local variables:
      INTEGER IX,IY                 !lattice indices
      CHARACTER*80 JUNK             !first line if of no interest
      LOGICAL SUCESS                !did we find a file to open?
      INTEGER CHOICE,YESNO          !get yesno input from user
      INTEGER MX,MY,MOMEGA          !parameters from PFILE
C Function:
      CHARACTER*40 CHARAC           !returns character input
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
10    CALL FLOPN2(PFILE,DUNIT,SUCESS)        !open the file for input
      MSTRNG(MINTS(IPFILE))=PFILE            !file may have been renamed
C
      IF (.NOT. SUCESS) THEN
         CALL REASK                          !prompt again for init P
         RETURN                              !if no file was found
      ELSE
C
        READ (DUNIT,5) JUNK                  !skip over title
        READ (DUNIT,5) JUNK                 
5       FORMAT (A)
C
C       lattice sizes must match; if they don't, allow for other options
        READ (DUNIT,*) MX,MY
        IF ((MX .NE. NX) .OR. (MY .NE. NY)) THEN
            CALL FLCLOS(PFILE,DUNIT)          !close it up
            WRITE (OUNIT,15) MX,MY
15          FORMAT (' Input file has does not have the correct'
     +      ' lattice size, it is ',I2,' by ',I2)
            CHOICE=YESNO(1,' Do you want to try another file?')
            IF (CHOICE .EQ. 0) THEN
                 CALL REASK                   !prompt again for init P
                 RETURN
            ELSE
               PFILE=CHARAC(PFILE,12, 'Enter another filename')
               MSTRNG(MINTS(IPFILE))=PFILE
               GOTO 10      !try to open this one
            END IF
        END IF            
C    
C       if we've gotten this far, we've opened a file with data
C       from a lattice of the same size; finally read in field
        READ (DUNIT,5) JUNK                  !skip over parameter values
        READ (DUNIT,5) JUNK                 
        READ (DUNIT,5) JUNK                 
C
        DO 100 IX=1,NX                            
          READ (DUNIT,110) (P(IX,IY),IY=1,NY)
100     CONTINUE
        DO 200 IX=1,NX                            
          READ (DUNIT,110) (XSI(IX,IY),IY=1,NY)
200     CONTINUE
110     FORMAT (5(2X,1PE14.7))
        CALL FLCLOS(PFILE,DUNIT)            !close file                
      END IF
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE REASK
C prompt again for initial P,XSI configuration; 
C called only if reading in from a file is request and failed
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'
      INCLUDE 'MAP.P6'
      INCLUDE 'PARAM.P6'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     redisplay initial Field Value Menu, but disallow choice 1
      MPRMPT(46)='1) (not allowed)'
      MLOLIM(49)=2
      MINTS(49)=2
      MREALS(49)=2
      CALL ASK(45,51)                  
C
      PTYPE=ABS(MREALS(IPTYPE))          !set parameter choices
C
      MPRMPT(46)='1) Read in starting values from a file'   !reset menu
      MLOLIM(49)=1
C            
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRMOUT(MUNIT,NLINES)
C write out parameter summary of length NLINES to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P6'
C Input variables:
      INTEGER MUNIT              !fortran unit number
      INTEGER NLINES             !number of lines sent to terminal (I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) CALL CLEAR
C
      WRITE (MUNIT,5)
      WRITE (MUNIT,6)
      IF (MUNIT .NE. GUNIT) THEN
          WRITE (MUNIT,10) NX,NY
      ELSE
          WRITE (MUNIT,*) NX,NY
      END IF
      WRITE (MUNIT,20) VOMEGA,SOMEGA
      WRITE (MUNIT,30) HFWID,LENGTH,FRONT
      WRITE (MUNIT,40) REYNLD
      IF (SUBLAT) WRITE (MUNIT,45) NXLL,NYLL,NXUR,NYUR
C
      IF (MUNIT .NE. GUNIT) THEN
         WRITE (MUNIT,*) ' '
         WRITE (MUNIT,50)
         WRITE (MUNIT,60)
         WRITE (MUNIT,70)
      END IF
C
      NLINES=8
      IF (SUBLAT) NLINES=NLINES+1
C            
5     FORMAT (' Output from project 6:')
6     FORMAT (' 2-D viscous incompressible'
     +        ' flow around a rectangular block')
10    FORMAT (' NX =',I3,5X,' NY =',I3)
20    FORMAT (' Vorticity relaxation= ',F5.3,5X,'Stream relaxation= ',
     +        F5.3)
30    FORMAT (' Plate half width=',I3,5X,' length='I3,5X,
     +        'and front edge='I3)
40    FORMAT (' Lattice Reynold''s number =',F7.3)
45    FORMAT (' Sublattice defined by  (',I3,',',I3,') and (',
     +        I3,',',I3,')')
50    FORMAT (3X,'Iter',4X,'Pressure',6X,'Viscous',11X,'Vorticity',
     +   10X,'Stream Function')
60    FORMAT (12X,'Force',9X,'Force',13X,'min,max',15X,'min,max')
70    FORMAT (3X,'----',4X,'--------',6X,'-------',11X,'---------',
     +   10X,'---------------')
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE 
     +  TXTOUT(MUNIT,NITER,VFORCE,PFORCE,XSIMAX,PMAX,XSIMIN,PMIN,NLINES)
C output forces and min, max field values to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:          
      INTEGER NITER               !number of iterations
      REAL VFORCE,PFORCE          !viscous force and pressure
      REAL XSIMAX,PMAX,XSIMIN,PMIN!min,max values of XSI, P
      INTEGER MUNIT               !fortran unit number
      INTEGER NLINES              !number of lines sent to terminal(I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     if screen is full, clear screen and retype headings          
      IF ((MOD(NLINES,TRMLIN-6) .EQ. 0) 
     +                          .AND. (MUNIT .EQ. OUNIT)) THEN
         CALL PAUSE('to continue...',1)
         CALL CLEAR
         WRITE (MUNIT,50)
         WRITE (MUNIT,60)
         WRITE (MUNIT,70)
      END IF
C             
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+1
      WRITE (MUNIT,40) NITER,PFORCE,VFORCE,XSIMIN,XSIMAX,PMIN,PMAX
C
40    FORMAT (2X,I5,2(2X,1PE12.5),2X,4(1PE10.3,1X))
50    FORMAT (3X,'Iter',4X,'Pressure',6X,'Viscous',11X,'Vorticity',
     +   10X,'Stream Function')
60    FORMAT (12X,'Force',9X,'Force',13X,'min,max',15X,'min,max')
70    FORMAT (3X,'----',4X,'--------',6X,'-------',11X,'---------',
     +   10X,'---------------')
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DISPLY(F,MUNIT,ITYPE,FMAX,FMIN)
C display stream function (P,ITYPE=STREAM) or vorticity
C (XSI,ITYPE=VRTCTY) as letters; positive values are capitals or 
C numbers, negative values are small letters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P6'
      INCLUDE 'IO.ALL'
C Input variables:
      REAL F(NXMAX,NYMAX)                !stream function or vorticity
      INTEGER MUNIT                      !unit we're writing to 
      INTEGER ITYPE                      !which field are we displaying?
      REAL FMAX,FMIN                     !min and max field values
C Local variables:
      INTEGER IX,IY                      !lattice indices
      INTEGER TEMP                       !field at current lattice site
      CHARACTER*1 FIELD(NXMAX)           !field as character data
      CHARACTER*80 BLNK                  !blanks for centering in X
      CHARACTER*1 ASKII(0:35),NEGASK(1:26)!charac data for display
      DATA BLNK /' '/
      DATA ASKII/'0','1','2','3','4','5','6','7','8','9',
     +           'A','B','C','D','E','F','G','H','I','J','K','L','M',
     +           'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA NEGASK/'a','b','c','d','e','f','g','h','i','j','k','l','m',
     +           'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA VRTCTY,STREAM /1,2/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) THEN
        DO 20 IY=1,YCNTR                   !center output
           WRITE (OUNIT,*) ' '
20      CONTINUE
      ELSE                                 !or display which field
        IF (ITYPE .EQ. STREAM) WRITE (MUNIT,*) ' stream function:'
        IF (ITYPE .EQ. VRTCTY) WRITE (MUNIT,*) ' vorticity:'
      END IF
C
      DO 100 IY=NY,1,-1
         DO 50 IX=1,NX
            IF (F(IX,IY) .GE. 0) THEN                      
               IF (FMAX .NE. 0) THEN
                    TEMP=NINT(F(IX,IY)*35/FMAX)   !scale field
               ELSE 
                    TEMP=0.
               END IF
               FIELD(IX)=ASKII(TEMP)              !convert to ascii
            ELSE
               IF (FMIN .NE. 0) THEN
                   TEMP=NINT(F(IX,IY)*26/FMIN)    !scale field
               ELSE              
                   TEMP=0
               END IF
               IF (TEMP .NE. 0) THEN
                   FIELD(IX)=NEGASK(TEMP)         !convert to ascii
               ELSE
                   FIELD(IX)=ASKII(TEMP)
               END IF
            END IF
C           leave blanks to indicate the plate
            IF (BNDCND(IX,IY) .EQ. 2) FIELD(IX)=BLNK(1:1)
50       CONTINUE    
C
C        write out a line at a time (no centering done for TUNIT)
         IF (MUNIT .EQ. TUNIT) THEN
                  WRITE (TUNIT,16) (FIELD(IX),IX=1,NX)
         ELSE
           IF (XSKIP) THEN
             WRITE (OUNIT,10) BLNK(1:XCNTR),(FIELD(IX),IX=1,NX)
           ELSE 
             WRITE (OUNIT,15) BLNK(1:XCNTR),(FIELD(IX),IX=1,NX)
           END IF
           IF (YSKIP) WRITE (OUNIT,*) ' '
         END IF
10       FORMAT (1X,A,100(A1,1X))                    
15       FORMAT (1X,A,100(A1))
16       FORMAT (1X,100(A1))
100   CONTINUE                                
C             
      RETURN 
      END  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFOUT(DEVICE,FIELD,F,FMIN,FMAX,FGRF,MX,MY)
C display contours of the sclaed stream function P and vorticity PSI
C to DEVICE 
C     the field values must be in an array that is exactly NX by NY;
C     this can be accomplished with implicit dimensioning which 
C     requires that FGRF and its dimensions be passed to this routine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'   
      INCLUDE 'PARAM.P6'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE                !which device
      REAL F(NXMAX,NYMAX)           !field
      INTEGER MX,MY                 !NX and NY in disguise
      REAL FGRF(MX,MY)              !field for graphing
      REAL FMIN,FMAX                !min,max field values
      INTEGER FIELD                 !which field
C Local variables:
      INTEGER SCREEN                !send to terminal
      INTEGER PAPER                 !make a hardcopy
      INTEGER FILE                  !send to a file
      INTEGER IX,IY                 !level index, lattice indices
      REAL PX(4),PY(4)              !edges of plate for graphing
      CHARACTER*9 CMIN,CMAX,CREY    !data as characters
      INTEGER LMIN,LMAX,RLEN        !length of string
      DATA SCREEN,PAPER,FILE/1,2,3/
      DATA VRTCTY,STREAM /1,2/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 10 IX=1,MX
         DO 15 IY=1,MY
           FGRF(IX,IY)=F(IX,IY)       !load field into FGRF
15       CONTINUE
10    CONTINUE
C
C     messages for the impatient
      IF ((DEVICE .NE. SCREEN) .AND. (FIELD .EQ. STREAM))
     +                  WRITE (OUNIT,100) 
C                           
C     calculate parameters for graphing
      NPLOT=2                       !how many plots
C 
      YMAX=MY
      YMIN=1.
      XMIN=1.
      XMAX=MX
      Y0VAL=XMIN
      X0VAL=YMIN
      NXTICK=5
      NYTICK=5
      NPOINT=5
      LABEL(1)='NX'
      LABEL(2)='NY'
      CALL CONVRT(FMIN,CMIN,LMIN)
      CALL CONVRT(FMAX,CMAX,LMAX)
      CALL CONVRT(REYNLD,CREY,RLEN)
C
      IF (FIELD .EQ. STREAM) THEN
         IPLOT=1
         INFO='Pmin='//CMIN(1:LMIN)//' Pmax='//CMAX(1:LMAX)
      ELSE IF (FIELD .EQ. VRTCTY) THEN
         IPLOT=2
         INFO='XSImin='//CMIN(1:LMIN)//' XSImax='//CMAX(1:LMAX)
      END IF
      TITLE='Stream Function and Vorticity, Reynold''s number ='
     +             //CREY(1:RLEN)
C 
      IF (FIELD .EQ. STREAM) THEN
        CALL GTDEV(DEVICE)                   !device nomination
        IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
      END IF
      CALL LNLNAX                          !draw axes
C 
C     draw in plate if it's big enough
      IF (LENGTH .GT. 2) THEN
        PX(1)=REAL(FRONT+1)
        PY(1)=1.
        PX(2)=PX(1)
        PY(2)=REAL(HFWID-1)
        PX(3)=REAL(BACK-1)
        PY(3)=PY(2)
        PX(4)=PX(3)
        PY(4)=1.
        NPOINT=4
        ILINE=2
        CALL XYPLOT(PX,PY)
      END IF
C
      CALL CONTOR(FGRF,MX,MY,FMIN,FMAX,NLEV)
C 
C     end graphing session
      IF (FIELD .EQ. VRTCTY) THEN
       IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)   !end graphics package
       IF (DEVICE .EQ. SCREEN) CALL TMODE         !switch to text mode
      END IF
C 
100   FORMAT (/,' Patience, please; output going to a file.')
C
      RETURN
      END                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE WRTOUT(P,XSI)
C write out stream function (P) and vorticity (XSI) to GUNIT for reading
C back in as initial conditions or for graphing with an external package
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P6'
C Input variables:
      REAL P(NXMAX,NYMAX)           !stream function
      REAL XSI(NXMAX,NYMAX)         !vorticity
C Local variables:
      INTEGER IX,IY                 !lattice indices
      INTEGER NLINES                !number of lines written to file
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL FLOPEN(GNAME,GUNIT)                  !open file
      MSTRNG(MINTS(IGNAME))=GNAME               !name may have changed
      CALL PRMOUT(GUNIT,NLINES)                 !write out header
C
      DO 100 IX=1,NX                            
        WRITE (GUNIT,10) (P(IX,IY),IY=1,NY)
100   CONTINUE
      DO 200 IX=1,NX                            
        WRITE (GUNIT,10) (XSI(IX,IY),IY=1,NY)
200   CONTINUE
10    FORMAT (5(2X,1PE14.7))
      CALL FLCLOS(GNAME,GUNIT)                  !close it up
C
      RETURN
      END
