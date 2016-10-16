CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM EXMPL6
C     Example 6: Solving Laplace's equation in two dimensions
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company, Inc.
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
      INCLUDE 'PARAM.E6'
      INCLUDE 'IO.ALL'
C Local variables:
      REAL P(NXMAX,NYMAX)            !field values
      REAL PGRF(NXMAX,NYMAX)         !field values for graphing
      INTEGER XMIN,XMAX,YMIN,YMAX    !corners of relaxing lattice
      REAL DELMAX,OLDMAX             !keep track of changes in P
      REAL ENERGY,DELE               !current and delta energy
      INTEGER NITER                  !number of iterations
      LOGICAL ITER                   !continue iterating?
      LOGICAL PRM                    !print out parameters?
      LOGICAL END                    !end this run?   
      INTEGER IX,IY                  !lattice indices
      INTEGER CHOICE                 !write out field now? 
      INTEGER NLINES                 !number of lines sent to terminal
      INTEGER SCREEN                 !send to terminal
      INTEGER PAPER                  !make a hardcopy
      INTEGER FILE                   !send to a file
C Functions:
      REAL GETFLT                    !user input function 
      INTEGER GETINT,YESNO           !user input functions
      LOGICAL LOGCVT                 !convert 0,1 to false,true
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END=.FALSE.                 !initial values
      DO 5 IX=1,NX
       DO 5 IY=1,NY
          P(IX,IY)=0
5     CONTINUE
C            
200   CONTINUE                    !allow for many runs with same lattice
        CALL PARAM(P,END)         !get new parameters
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
        PRM=.TRUE.
        NITER=0
C
99      CONTINUE                   !begin iterations
          IF ((TTERM) .AND. (PRM)) CALL PRMOUT(OUNIT,NLINES)
          NITER=NITER+1
C                                  !relax lattice
          CALL RELAX(P,XMIN,XMAX,YMIN,YMAX,ENERGY,DELE,DELMAX,OLDMAX)
          IF (NITER .EQ. 1) THEN   !no changes if this is the first time
              DELE=0
              OLDMAX=0
          END IF
C          
          IF (TFILE) CALL TXTOUT(TUNIT,NITER,ENERGY,DELE,DELMAX,NLINES)
          IF (TTERM) CALL TXTOUT(OUNIT,NITER,ENERGY,DELE,DELMAX,NLINES)
C
          ITER=.FALSE.
          PRM=.TRUE.
          IF (MOD(NITER,NFREQ) .NE. 0) THEN         
             ITER=.TRUE.        !continue iterating without asking
             PRM=.FALSE.        !don't print out header
          ELSE                  !display fields every NFREQ iteration
             IF (GTERM) THEN
                CALL PAUSE('to see the field ...',1)
                CALL CLEAR                                             
                CALL GRFOUT(SCREEN,P,PGRF,NX,NY,ENERGY)
             ELSE IF (TTERM) THEN
                CALL PAUSE('to see the field ...',1)
                CALL CLEAR                                             
                CALL DISPLY(P,OUNIT)
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
              IF (GFILE) CALL WRTOUT(P)    !write out field if requested
              IF (GHRDCP) CALL GRFOUT(PAPER,P,PGRF,NX,NY,ENERGY)
              IF (TFILE) CALL DISPLY(P,TUNIT)
           END IF
        END IF
        GOTO 200       !display menu
       END IF
C     
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RELAX(P,XMIN,XMAX,YMIN,YMAX,ENERGY,DELE,DELMAX,OLDMAX)
C relaxes the lattice (specified by XMIN, etc), calculates the new 
C energy (ENERGY) and maximum change in the field (DELMAX)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E6'
      INCLUDE 'IO.ALL'
C Input/output variables:
      REAL P(NXMAX,NYMAX)            !field values (I/O)
      INTEGER XMIN,XMAX,YMIN,YMAX    !corners of relaxing lattice(input)
      REAL DELMAX,OLDMAX             !keep track of change in P (output)
      REAL ENERGY,DELE               !current and delta energy (output)
C Local variables:
      INTEGER IX,IY                  !lattice indices
      REAL A,B,C,D                   !values of P at neighboring points
      REAL PNEW,POLD,DELP            !temp value of P at current point
      REAL OLDE                      !old energy
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      OLDMAX=DELMAX                  !roll values
      OLDE=ENERGY                                  
      DELMAX=0.                      !initialize values
      ENERGY=0.
C
      DO 200 IY=YMIN,YMAX
         DO 100 IX=XMIN,XMAX
            A=0.
            B=0.
            C=0.
            D=0.
            IF (IX .LT. NX) A=P(IX+1,IY)  !field values at neighboring
            IF (IX .GT. 1)  B=P(IX-1,IY)  !lattice points
            IF (IY .LT. NY) C=P(IX,IY+1)
            IF (IY .GT. 1)  D=P(IX,IY-1)
C
            POLD=P(IX,IY)
            IF (BNDCND(IX,IY) .EQ. NON) THEN
              PNEW=(1.-OMEGA)*POLD+OMEGA/4*(A+B+C+D)    !relax
              P(IX,IY)=PNEW
              IF (POLD .NE. 0.) THEN
                  DELP=ABS((POLD-PNEW)/POLD)
              ELSE
                  DELP=ABS(PNEW)
              END IF
              IF (DELP .GT. DELMAX) DELMAX=DELP         !find max change
            ELSE
              PNEW=POLD                  !don't relax DIRCHLT bound cond
            END IF
C
            IF (IX .GT. 1) ENERGY=ENERGY+(PNEW-B)**2   
            IF (IY .GT. 1) ENERGY=ENERGY+(PNEW-D)**2
C
100     CONTINUE
200   CONTINUE
      ENERGY=ENERGY/2/(NX-1)/(NY-1)
      DELE=OLDE-ENERGY                
C
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
      INCLUDE 'PARAM.E6'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)         
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get environment parameters
      CALL SETUP                
C 
C     display header screen     
      DESCRP(1)= 'EXAMPLE 6'
      DESCRP(2)= 'Solving Laplace''s Equation in Two Dimensions'
      NHEAD=2
C      
C     text output description
      DESCRP(3)= 'iteration, energy, change in energy,'
      DESCRP(4)= 'and the maximum change in field' 
      NTEXT=2
C 
C     graphics output description
      DESCRP(5)= 'field values'
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
      MTYPE(13)=FLOAT
      MPRMPT(13)= 'Enter fixed value of the field for UPPER boundary'
      MTAG(13)= 'Field on UPPER boundary'
      MLOLIM(13)=FMIN
      MHILIM(13)=FMAX
      MREALS(13)=FMIN
C                
      MTYPE(14)=FLOAT
      MPRMPT(14)= 'Enter fixed value of the field for LOWER boundary'
      MTAG(14)= 'Field on LOWER boundary'
      MLOLIM(14)=FMIN
      MHILIM(14)=FMAX
      MREALS(14)=FMAX
C                
      MTYPE(15)=FLOAT
      MPRMPT(15)= 'Enter fixed value of the field for LEFT boundary'
      MTAG(15)= 'Field on LEFT boundary'
      MLOLIM(15)=FMIN
      MHILIM(15)=FMAX
      MREALS(15)=FMAX
C                
      MTYPE(16)=FLOAT
      MPRMPT(16)= 'Enter fixed value of the field for RIGHT boundary'
      MTAG(16)= 'Field on RIGHT boundary'
      MLOLIM(16)=FMIN
      MHILIM(16)=FMAX
      MREALS(16)=FMIN                               
C                
      MTYPE(17)=NOSKIP
      MPRMPT(17)=
     +  'Do you wish to enter/review/display interior bound cond?'
      MTAG(17)='Enter/review/display interior boundary conditions'
      MINTS(17)=0
      MREALS(17)=35
C
      MTYPE(18)=QUIT
C
C     This section of the menu (from 25-33) is called from BNDRY
      MTYPE(25)=MTITLE
      MPRMPT(25)='Interior Boundary Conditions Menu'
      MLOLIM(25)=2
      MHILIM(25)=1
C 
      MTYPE(26)=MTITLE
      MPRMPT(26)='1) input rectangular boundary conditions'
      MLOLIM(26)=0
      MHILIM(26)=0
C 
      MTYPE(27)=MTITLE
      MPRMPT(27)='2) input 45 degree line boundary conditions'
      MLOLIM(27)=0            
      MHILIM(27)=0
C 
      MTYPE(28)=MTITLE
      MPRMPT(28)='3) review boundary conditions '
      MLOLIM(28)=0
      MHILIM(28)=0
C 
      MTYPE(29)=MTITLE
      MPRMPT(29)='4) display boundary conditions '
      MLOLIM(29)=0
      MHILIM(29)=0
C 
      MTYPE(30)=MTITLE
      MPRMPT(30)='5) return to main menu' 
      MLOLIM(30)=0
      MHILIM(30)=0
C 
      MTYPE(31)=MCHOIC
      MPRMPT(31)='Enter Choice'
      MTAG(31)='32 32 32 32 32 32'
      MLOLIM(31)=1
      MHILIM(31)=5
      MINTS(31)=5
      MREALS(31)=-5
C 
      MTYPE(38)=FLOAT
      MPRMPT(38)= 'Enter relaxation parameter'
      MTAG(38)= 'Relaxation parameter'
      MLOLIM(38)=0.
      MHILIM(38)=2.                       
      MREALS(38)=1.
C            
      MTYPE(39)=NOSKIP
      MPRMPT(39)='Do you want to specify a sublattice to relax first?'
      MTAG(39)='Relax a sublattice first?'
      MINTS(39)=0
      MREALS(39)=44.
C                     
      MTYPE(40)=NUM
      MPRMPT(40)='Enter lower left X value for sublattice'
      MTAG(40)='Lower left X sublattice value'
      MLOLIM(40)=1
      MHILIM(40)=NXMAX
      MINTS(40)=1
C
      MTYPE(41)=NUM
      MPRMPT(41)='Enter lower left Y value for sublattice'
      MTAG(41)='Lower left Y sublattice value'
      MLOLIM(41)=1
      MHILIM(41)=NYMAX
      MINTS(41)=1
C               
      MTYPE(42)=NUM
      MPRMPT(42)='Enter upper right X value for sublattice'
      MTAG(42)='Upper right X sublattice value'
      MLOLIM(42)=1
      MHILIM(42)=NXMAX
      MINTS(42)=1
C               
      MTYPE(43)=NUM
      MPRMPT(43)='Enter upper right Y value for sublattice'
      MTAG(43)='Upper right Y sublattice value'
      MLOLIM(43)=1
      MHILIM(43)=NXMAX
      MINTS(43)=1
C
      MTYPE(44)=MTITLE
      MPRMPT(44)='Initial Field Value Menu'
      MLOLIM(44)=1
      MHILIM(44)=1
C
      MTYPE(45)=MTITLE
      MPRMPT(45)='1) Read in initial values from a file'
      MLOLIM(45)=0
      MHILIM(45)=0
C                
      MTYPE(46)=MTITLE
      MPRMPT(46)='2) Set field at all lattice points to one value'
      MLOLIM(46)=0            
      MHILIM(46)=0
C 
      MTYPE(47)=MTITLE
      MPRMPT(47)='3) Leave field values unchanged'
      MLOLIM(47)=0
      MHILIM(47)=0
C                   
      MTYPE(48)=MCHOIC
      MPRMPT(48)='Enter Choice'
      MTAG(48)='49 51 52'
      MLOLIM(48)=1
      MHILIM(48)=3
      MINTS(48)=2
      MREALS(48)=2
C                        
      MTYPE(49)=CHSTR
      MPRMPT(49)= 'Enter name of data file'
      MTAG(49)= 'File with initial values for the field'
      MHILIM(49)=12.
      MINTS(49)=3.
      MSTRNG(MINTS(49))= 'exmpl6.in'
C
      MTYPE(50)=SKIP
      MREALS(50)=60.
C                     
      MTYPE(51)=FLOAT
      MPRMPT(51)= 'Enter initial value for the field'
      MTAG(51)= 'Initial Field Value'
      MLOLIM(51)=FMIN
      MHILIM(51)=FMAX
      MREALS(51)=(FMAX-FMIN)/2
C
      MTYPE(52)=SKIP
      MREALS(52)=60.
C
      MSTRNG(MINTS(75))= 'exmpl6.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C 
      MSTRNG(MINTS(86))= 'exmpl6.grf'
C                     
      MTYPE(87)=NUM
      MPRMPT(87)= 'Enter the display frequency for the field'
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
C field as ascii characters based on lattice size and terminal size;
C resets all boundary conditions and default menu values
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E6'
      INCLUDE 'IO.ALL'             
      INCLUDE 'MENU.ALL'
      INCLUDE 'MAP.E6'
C Local variables:
      INTEGER END                   !end program
      INTEGER IX,IY,IBC             !lattice indices, BC index
      INTEGER NXHI,NYHI,NXLO,NYLO   !limits on lattice size
      INTEGER NXDEF,NYDEF           !default lattice sizes
      LOGICAL RESET                 !reset parameters? 
C Functions:
      INTEGER YESNO,GETINT          !user input functions
      LOGICAL LOGCVT                !change 1,0 to true and false
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     allow user to end the program
      CALL CLEAR
      IF (.NOT. FIRST) THEN
         END=YESNO(0,' Do you want to end the program?')
         IF (END .EQ. 1) CALL FINISH
      ELSE
C        the lattice size is determined by array size and terminal size;
C        if you're using graphics, terminal size won't matter, 
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
C     set all bound cond parameters to zero
      NBC=0
      DO 10 IBC=1,BCMAX
          XLL(IBC)=0
          YLL(IBC)=0
          XUR(IBC)=0
          YUR(IBC)=0
          BCGEOM(IBC)=0
          BCVAL(IBC)=0.
10    CONTINUE
      DO 20,IX=1,NXMAX
         DO 30 IY=1,NYMAX
            BNDCND(IX,IY)=0      
30       CONTINUE
20    CONTINUE
C
C     set up limits on sublattice
      MHILIM(INXLL)=NX          
      MHILIM(INYLL)=NY
      MHILIM(INXUR)=NX
      MHILIM(INYUR)=NY
      MINTS(40)=1        !lower left  x sublattice
      MINTS(41)=1        !upper right x sublattice
      MINTS(42)=1        !lower left  y sublattice 
      MINTS(43)=1        !upper right y sublattice 
C
C     allow for resetting of defaults
      IF (FIRST) THEN
         FIRST=.FALSE.
      ELSE
         RESET=LOGCVT(YESNO(0,' Do you want to reset default values?'))
         IF (RESET) THEN
           MREALS(13)=FMIN    !upper bound. cond.
           MREALS(14)=FMAX    !lower bound. cond.
           MREALS(15)=FMAX    !left  bound. cond.
           MREALS(16)=FMIN    !right bound. cond.
           MREALS(38)=1.      !omega
           MINTS(39)=0        !sublattice?
           MINTS(48)=2        !PTYPE
           MREALS(48)=2
           MSTRNG(MINTS(49))= 'exmpl6.in'
           MREALS(51)=(FMAX-FMIN)/2 !initial P value
           MINTS(73)=TXTTRM
           MSTRNG(MINTS(75))= 'exmpl6.txt'
           MINTS(74)=TXTFIL
           MINTS(83)=GRFTRM
           MINTS(84)=GRFHRD
           MINTS(85)=GRFFIL
           MSTRNG(MINTS(86))= 'exmpl6.grf'
           MINTS(87)= 10       !graphing frequency
           MINTS(88)= 10       !number of contour lines
         END IF
      END IF
C
      RETURN 
      END                                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PARAM(P,END)
C gets parameters from screen                                      
C ends program on request
C closes old files
C maps menu variables to program variables              
C opens new files
C calculates all derivative parameters
C performs checks on sublattice parameters
C set the field to its initial values (P)
C and controls ending of program (END)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'                                            
      INCLUDE 'PARAM.E6'
      INCLUDE 'MAP.E6'
C Input/output variables:
      REAL P(NXMAX,NYMAX)     !field values
      LOGICAL END             !end program?
C Local variables:
      INTEGER TYPE,GEOM       !type and geom of bound cond
      INTEGER IX,IY           !lattice indices
C Functions:
      LOGICAL LOGCVT          !converts 1 and 0 to true and false 
      INTEGER GETINT          !get integer from screen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get input from terminal
10    CONTINUE
        MINTS(INTERR)=0       !reset value
C       get all parameters except interior bound cond
        CALL CLEAR                
        CALL ASK(1,ISTOP)
C
C       start fresh or end altogether, if so requested
        IF (MREALS(IMAIN) .EQ. STOP)  THEN
           END=.TRUE.
           RETURN
        END IF
C
C       set basic boundary conditions before doing interior bc
        UPPER=MREALS(IUPPER)
        LOWER=MREALS(ILOWER)
        LEFT=MREALS(ILEFT)
        RIGHT=MREALS(IRIGHT)
        GEOM=RECTGL
        TYPE=DRCHLT
        CALL SETBC(P,GEOM,TYPE,UPPER,1,NX,NY,NY)
        CALL SETBC(P,GEOM,TYPE,LOWER,1,NX,1,1)
        CALL SETBC(P,GEOM,TYPE,LEFT,1,1,1,NY)
        CALL SETBC(P,GEOM,TYPE,RIGHT,NX,NX,1,NY)
C       need to know if graphics are available
        GTERM=LOGCVT(MINTS(IGTERM))
C
C       allow for input of interior boundary conditions 
        IF (MINTS(INTERR) .EQ. 1) CALL BNDRY(P)
      IF (MINTS(INTERR) .EQ. 1) GOTO 10
C
C     close files if necessary
      IF (TNAME .NE. MSTRNG(MINTS(ITNAME))) 
     +     CALL FLCLOS(TNAME,TUNIT)
C                                           
C     physical and numerical parameters
      OMEGA=MREALS(IOMEGA)                    
      SUBLAT=LOGCVT(MINTS(ISUB))
      NXLL=MINTS(INXLL)
      NYLL=MINTS(INYLL)
      NXUR=MINTS(INXUR)
      NYUR=MINTS(INYUR)
      PTYPE=ABS(MREALS(IPTYPE))
      PFILE=MSTRNG(MINTS(IPFILE))
      PINIT=MREALS(IPINIT)
C
C     text output
      TTERM=LOGCVT(MINTS(ITTERM))
      TFILE=LOGCVT(MINTS(ITFILE))
      TNAME=MSTRNG(MINTS(ITNAME))
C
C     graphics output
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
         CALL ASK(INXLL,INYUR)
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
C     set P to its initial value
      IF (PTYPE .EQ. 1) CALL READIN(P)         !read data in
C     sometimes in READIN PTYPE is changed
      IF (PTYPE .EQ. 2) THEN
        DO 60 IX=1,NX
           DO 70 IY=1,NY         !set all non bound cond points to PINIT
              IF (BNDCND(IX,IY) .EQ. 0) P(IX,IY)=PINIT
70         CONTINUE
60      CONTINUE
      ELSE IF (PTYPE .EQ. 3) THEN              !keep P the same
        CONTINUE      
      END IF                                                     
C            
C     reset PTYPE menu parameters so that P remains unchanged during
C     runs unless the user explicitly requests otherwise 
      MINTS(IPTYPE)=3
      MREALS(IPTYPE)=3              
C                                   
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BNDRY(P)
C obtains interior boundary conditions from user (P)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E6'
      INCLUDE 'MAP.E6'
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'
C Input/output variables:
      REAL P(NXMAX,NYMAX)            !field values  
C Local variables:
      INTEGER TYPE                   !type of bound cond
      INTEGER XLDIST,XRDIST          !dist of diag to left and right bc
      INTEGER YUDIST,YLDIST          !dist of diag to upper and lower bc
      INTEGER UPLIM                  !upper limit to diag length
C Functions:
      INTEGER GETINT                 !get integer from screen
      REAL GETFLT                    !get float from screen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      TYPE=DRCHLT               !currently only dirichlet bc implemented
C
1     CALL CLEAR                  
      CALL ASK(IBCMNU,IBNDRY)                  !display bound cond menu
C                                                             
C       get detailed geometrical information for rectangular b.c.
        IF (MREALS(IBNDRY) .EQ. -1)  THEN    
          IF (NBC .LT. BCMAX) THEN
            NBC=NBC+1                          !update number of b.c.
            BCGEOM(NBC)=RECTGL                 !set geometry
            XLL(NBC)= GETINT(2,2,NX-1,' Enter lower left X value')
            YLL(NBC)= GETINT(2,2,NY-1,' Enter lower left Y value')
            XUR(NBC)= 
     +      GETINT(XLL(NBC),XLL(NBC),NX-1,' Enter upper right X value')
            YUR(NBC)= 
     +      GETINT(YLL(NBC),YLL(NBC),NY-1,' Enter upper right Y value')
            BCVAL(NBC)=GETFLT(FMIN,FMIN,FMAX,' Enter field value ')
            CALL SETBC(P,BCGEOM(NBC),TYPE,BCVAL(NBC),XLL(NBC),XUR(NBC),
     +                 YLL(NBC),YUR(NBC))
           ELSE                        !too many b.c.'s
            WRITE (OUNIT,30)           !display warning message
            WRITE (OUNIT,40)
            CALL PAUSE('to continue...',1)
           END IF 
C
C       get detailed geometrical information for diagonal b.c.
        ELSE IF (MREALS(IBNDRY) .EQ. -2) THEN       
          IF (NBC .LT. BCMAX) THEN
            NBC=NBC+1                  !update number of  b.c.'s  
            BCGEOM(NBC)=DIAGNL         !set geometry
            XLL(NBC)= GETINT(2,2,NX-1,' Enter first X value ')
            YLL(NBC)= GETINT(2,2,NY-1,' Enter first Y value ')
            XUR(NBC)= GETINT(1,1,4,'Enter quadrant: 1, 2, 3, or 4')
C           make sure that diagonal stays within the lattice
            XRDIST=NX-1-XLL(NBC)
            XLDIST=XLL(NBC)-2
            YUDIST=NY-1-YLL(NBC)
            YLDIST=YLL(NBC)-2
            IF (XUR(NBC) .EQ. 1) THEN
                UPLIM=MIN(XRDIST,YUDIST)
            ELSE IF (XUR(NBC) .EQ. 2) THEN
                UPLIM=MIN(XLDIST,YUDIST)
            ELSE IF (XUR(NBC) .EQ. 3) THEN
                UPLIM=MIN(XLDIST,YLDIST)
            ELSE IF (XUR(NBC) .EQ. 4) THEN
                UPLIM=MIN(XRDIST,YLDIST)
            END IF
            YUR(NBC)= GETINT(UPLIM,0,UPLIM,
     +         'Enter length of 45 degree line (max is default)')
            BCVAL(NBC)=GETFLT(FMIN,FMIN,FMAX,' Enter field value ')
            CALL SETBC(P,BCGEOM(NBC),TYPE,BCVAL(NBC),XLL(NBC),XUR(NBC),
     +                 YLL(NBC),YUR(NBC))
           ELSE                       !too many b.c.'s
            WRITE (OUNIT,30)          !display warning message
            WRITE (OUNIT,40)
            CALL PAUSE('to continue...',1)
C
           END IF 
C
        ELSE IF (MREALS(IBNDRY) .EQ. -3) THEN     !review bc
            IF (NBC .EQ. 0) THEN
                WRITE(OUNIT,2)                    !nothing to review
                CALL PAUSE('to continue...',1)
            ELSE
                CALL REVIEW(P)       
            END IF
C
        ELSE IF (MREALS(IBNDRY) .EQ. -4) THEN     !display bc
            IF (NBC .EQ. 0) THEN
                WRITE(OUNIT,2)                    !nothing to display
                CALL PAUSE('to continue...',1)
            ELSE
               IF (GTERM) THEN                    !graphics   
                  CALL GRFBC
               ELSE
                  CALL DSPLBC(P)                  !no graphics
               END IF
            END IF
C
        ELSE IF (MREALS(IBNDRY) .EQ. -5) THEN     !go back to Main menu
            RETURN             
        END IF
C
      GOTO 1
C
30    FORMAT ( ' You have entered the maximum number of boundary'
     +         ' conditions allowed')
40    FORMAT ( ' You can add more only if you delete others first'
     +         ' using the REVIEW option')
2     FORMAT (' No interior boundary conditions have been entered')
C       
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SETBC(P,GEOM,TYPE,PVALUE,X1,X2,Y1,Y2) 
C given an interior boundary condition whose geometry is given by GEOM,
C X1, X2, Y1, Y2; set P array elements to PVALUE and set BNCCND to TYPE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E6'
C Input variables:
      REAL P(NXMAX,NYMAX)                !field values (I/O)
      INTEGER GEOM                       !geometry type of bound cond
      INTEGER TYPE                       !type of bound cond
      REAL PVALUE                        !value of P on bound
      INTEGER X1,X2,Y1,Y2                !specify location of bound cond
C Local variables:
      INTEGER IX,IY,I                    !lattice indices
      INTEGER XSLOPE,YSLOPE              !slopes for diagonal bc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (GEOM .EQ. RECTGL) THEN
         DO 10 IX=X1,X2
            P(IX,Y1)=PVALUE        !lower side
            P(IX,Y2)=PVALUE        !upper side
            BNDCND(IX,Y1)=TYPE
            BNDCND(IX,Y2)=TYPE
10       CONTINUE
         DO 20 IY=Y1,Y2
            P(X1,IY)=PVALUE        !left side
            P(X2,IY)=PVALUE        !right side
            BNDCND(X1,IY)=TYPE
            BNDCND(X2,IY)=TYPE
20       CONTINUE
C
C     for diagonal, Y2 is length and X2 is angle
      ELSE IF (GEOM .EQ. DIAGNL) THEN
          IF ((X2 .EQ. 1) .OR. (X2 .EQ. 2)) YSLOPE=1
          IF ((X2 .EQ. 3) .OR. (X2 .EQ. 4)) YSLOPE=-1
          IF ((X2 .EQ. 1) .OR. (X2 .EQ. 4)) XSLOPE=1
          IF ((X2 .EQ. 2) .OR. (X2 .EQ. 3)) XSLOPE=-1
          DO 30 I=0,Y2
             IY=Y1+I*YSLOPE
             IX=X1+I*XSLOPE
             P(IX,IY)=PVALUE
C display boundary conditions (P) as small letters
C all non-bound conditions are displayed as '-', regardless of value
C so user may clearly see boundary conditions; always written to OUNIT
             BNDCND(IX,IY)=TYPE
30        CONTINUE
      END IF
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DSPLBC(P)
C display boundary conditions (P) as small letters
C all non-bound conditions are displayed as '-', regardless of value
C so user may clearly see boundary conditions; always written to OUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E6'
      INCLUDE 'IO.ALL'
C Input variables:
      REAL P(NXMAX,NYMAX)                !field values
C Local variables:
      INTEGER IX,IY                      !lattice indices
      INTEGER PTEMP                      !field at current lattice site
      CHARACTER*1 FIELD(NXMAX)           !field as character data
      CHARACTER*80 BLNK                  !blanks for centering in X
      CHARACTER*1 ASKII,NEGASK(0:25)     !charac data for display
      DATA BLNK /' '/
      DATA ASKII/'-'/
      DATA NEGASK/'a','b','c','d','e','f','g','h','i','j','k','l','m',
     +           'n','o','p','q','r','s','t','u','v','w','x','y','z'/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL CLEAR
      DO 20 IY=1,YCNTR                      !center output
         WRITE (OUNIT,*) ' '
20    CONTINUE
C
      DO 100 IY=NY,1,-1
         DO 50 IX=1,NX                                                  
            IF (BNDCND(IX,IY) .EQ. 0) THEN  !set non-b.c. to '-'
               FIELD(IX)=ASKII
            ELSE
               PTEMP=NINT(P(IX,IY))         !change b.c. to ascii values
               IF (PTEMP .GT. FMAX) PTEMP=FMAX
               IF (PTEMP .LT. FMIN) PTEMP=FMIN
               FIELD(IX)=NEGASK(PTEMP)
            END IF
50       CONTINUE    
         IF (XSKIP) THEN
           WRITE (OUNIT,10) BLNK(1:XCNTR),(FIELD(IX),IX=1,NX)
         ELSE 
           WRITE (OUNIT,15) BLNK(1:XCNTR),(FIELD(IX),IX=1,NX)
         END IF
         IF (YSKIP) WRITE (OUNIT,*) ' '
10       FORMAT (1X,A,100(A1,1X))
15       FORMAT (1X,A,100(A1))
100   CONTINUE
      CALL PAUSE(' to go back to menu...',1)
C
      RETURN
      END  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFBC
C display boundary conditions as dashed lines so user may inspect
C boundary conditions; always written to OUNIT using graphics package
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'   
      INCLUDE 'PARAM.E6'
      INCLUDE 'GRFDAT.ALL'
C Local variables:
      REAL BCX(5),BCY(5)            !corners of boundary conditions
      INTEGER I,IX,IY               !level index, lattice indices
      INTEGER IBC                   !index on boundary conditions
      REAL XSLOPE, YSLOPE           !slope of diagonal BC
      INTEGER SCREEN                !send to terminal
      INTEGER PAPER                 !make a hardcopy
      INTEGER FILE                  !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NPLOT=1                              !how many plots
      IPLOT=1
C 
      ILINE=2
      YMAX=NY                              !axes data
      YMIN=1.
      XMIN=1.
      XMAX=NX
      Y0VAL=XMIN
      X0VAL=YMIN
      NXTICK=5
      NYTICK=5                             
C
      LABEL(1)='NX'                        !descriptions
      LABEL(2)='NY'
      TITLE='Boundary conditions'
C                
      CALL GTDEV(SCREEN)                   !device nomination
      CALL GMODE                           !change to graphics mode
      CALL LNLNAX                          !draw axes
C 
C     display interior boundary conditions as dashed lines
      DO 200 IBC=1,NBC                     !loop over b.c.
         IF (BCGEOM(IBC) .EQ. DIAGNL) THEN  
           NPOINT=2
           IF ((XUR(IBC) .EQ. 1) .OR. (XUR(IBC) .EQ. 2)) YSLOPE=1
           IF ((XUR(IBC) .EQ. 3) .OR. (XUR(IBC) .EQ. 4)) YSLOPE=-1
           IF ((XUR(IBC) .EQ. 1) .OR. (XUR(IBC) .EQ. 4)) XSLOPE=1
           IF ((XUR(IBC) .EQ. 2) .OR. (XUR(IBC) .EQ. 3)) XSLOPE=-1
           BCX(1)=XLL(IBC)
           BCX(2)=XLL(IBC)+XSLOPE*YUR(IBC) 
           BCY(1)=YLL(IBC)
           BCY(2)=YLL(IBC)+YSLOPE*YUR(IBC)
         ELSE IF (BCGEOM(IBC) .EQ. RECTGL) THEN
           NPOINT=5
           BCX(1)=XLL(IBC)
           BCX(2)=BCX(1)
           BCX(3)=XUR(IBC)
           BCX(4)=BCX(3)
           BCY(1)=YLL(IBC)
           BCY(2)=YUR(IBC)
           BCY(3)=BCY(2)  
           BCY(4)=BCY(1)
           BCX(5)=BCX(1)
           BCY(5)=BCY(1)
         END IF
         CALL XYPLOT(BCX,BCY)          !plot boundaries 
200   CONTINUE
C
C     end graphing session
      CALL GPAGE(SCREEN)   !end graphics package
      CALL TMODE           !switch to text mode
C 
      RETURN
      END                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE REVIEW(P)
C routine to allow user to review interior boundary conditions;
C they can 1) leave as is , 2) delete, 3) alter field value
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E6'
      INCLUDE 'IO.ALL'
C Input/output variables:
      REAL P(NXMAX,NYMAX)        !field values        
C Local variables:
      INTEGER IBC,JBC            !index on bound cond
      INTEGER GETINT,CHOICE      !get user choice
      LOGICAL DELETE(BCMAX)      !delete this bound cond? 
      INTEGER NDEL               !number deleted
      INTEGER TYPE               !type of boundary condition
      REAL GETFLT                !get float from user
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NDEL=0 
C
      DO 100 IBC=1,NBC                  
         DELETE(IBC)=.FALSE.          !initialize
C
C        write out description of this bound cond
         WRITE (OUNIT,5) IBC
         IF (BCGEOM(IBC) .EQ. DIAGNL) THEN  
            WRITE (OUNIT,10) XLL(IBC),YLL(IBC),XUR(IBC),YUR(IBC)
         ELSE IF (BCGEOM(IBC) .EQ. RECTGL) THEN
            WRITE (OUNIT,15) XLL(IBC),YLL(IBC),XUR(IBC),YUR(IBC)
         END IF 
         WRITE (OUNIT,20) BCVAL(IBC)
C
         CHOICE=GETINT(1,1,3,
     +   'Do you want to 1)keep it 2)delete it 3) change field value?')
C
         IF (CHOICE .EQ. 1) THEN       !do nothing
            CONTINUE                               
         ELSE IF (CHOICE .EQ. 2) THEN  !delete it
            NDEL=NDEL+1
            DELETE(IBC)=.TRUE.
            TYPE= NON                 !releases b.c.
            CALL SETBC(P,BCGEOM(IBC),TYPE,BCVAL(IBC),XLL(IBC),XUR(IBC),
     +                 YLL(IBC),YUR(IBC))
         ELSE IF (CHOICE .EQ. 3) THEN  !get new field value
            BCVAL(IBC)=GETFLT(FMIN,FMIN,FMAX,' Enter field value ')
         END IF
100   CONTINUE
C
C   get rid of spaces left in bound cond set when one or more is deleted
      IF (NDEL .NE. 0) THEN            
       IBC=0                      !JBC labels old set of bound cond
       DO 200 JBC=1,NBC           !IBC labels new set of bound cond
          IF (DELETE(JBC)) THEN
             CONTINUE
          ELSE
             IBC=IBC+1            !next empty place
             XLL(IBC)=XLL(JBC)    !move up parameter values to 
             YLL(IBC)=YLL(JBC)    !next empty place
             XUR(IBC)=XUR(JBC)
             YUR(IBC)=YUR(JBC)
             BCGEOM(IBC)=BCGEOM(JBC)
             BCVAL(IBC)=BCVAL(JBC)
          END IF
200    CONTINUE
       NBC=NBC-NDEL               !update number of b.c.'s
      END IF
C            
      DO 300 IBC=1,NBC
C        reset boundary conditions in case some of deleted 
C        bound cond wrote over others
C        also adjusts those for which the field value was changed
         TYPE=DRCHLT
         CALL SETBC(P,BCGEOM(IBC),TYPE,BCVAL(IBC),XLL(IBC),XUR(IBC),
     +                 YLL(IBC),YUR(IBC))
300   CONTINUE
      RETURN
C
5     FORMAT (' Boundary Condition Number ',I2)
10    FORMAT (' Diagonal with first coordinates=(',I2,',',I2,
     +        '),  quadrant =',I2,', and length='I2)
15    FORMAT (' Rectangle with lower left at (',I2,',',I2,
     +        ') and upper right at (',I2,',',I2,')')
20    FORMAT (' The field is fixed at ',1PE12.5)
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE READIN(P)
C read in field values and boundary conditions;
C reads in only files written by subroutine WRTOUT, and files written
C with same lattice size
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E6'
      INCLUDE 'MENU.ALL'
      INCLUDE 'MAP.E6'
C Output variables:
      REAL P(NXMAX,NYMAX)           !field values
C Local variables:
      INTEGER IX,IY                 !lattice indices
      CHARACTER*80 JUNK             !first line if of no interest
      LOGICAL SUCESS                !did we find a file to open?
      INTEGER CHOICE                !user answer
      INTEGER MX,MY,MOMEGA          !parameters from PFILE
C Function:
      INTEGER YESNO                 !get yesno input from user
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
5       FORMAT (A)
C
C       lattice sizes must match; if they don't, allow for other options
        READ (DUNIT,*) MX,MY,MOMEGA          
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
        DO 100 IX=1,NX                            
          READ (DUNIT,110) (P(IX,IY),IY=1,NY)
100     CONTINUE
110     FORMAT (5(2X,1PE14.7))
        CALL FLCLOS(PFILE,DUNIT)            !close file
      END IF
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE REASK
C prompt again for initial P configuration; 
C called only if reading in from a file is requested and failed
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'
      INCLUDE 'MAP.E6'
      INCLUDE 'PARAM.E6'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     redisplay initial Field Value Menu, but disallow choice 1
      MPRMPT(45)='1) (not allowed)'
      MLOLIM(48)=2
      MINTS(48)=3
      MREALS(48)=3
      CALL ASK(44,52)                  
C
      PTYPE=ABS(MREALS(IPTYPE))          !set parameter choices
      PINIT=MREALS(IPINIT)
C
      MPRMPT(45)='1) Read in initial values from a file'   !reset menu
      MLOLIM(48)=1
C            
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRMOUT(MUNIT,NLINES)
C write out parameter summary to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E6'
C Input variables:
      INTEGER MUNIT                   !fortran unit number
      INTEGER NLINES                  !number of lines sent to terminal
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) CALL CLEAR
C
      WRITE (MUNIT,5)
      IF (MUNIT .EQ. GUNIT) THEN
         WRITE (MUNIT,*) NX,NY,OMEGA
      ELSE
         WRITE (MUNIT,10) NX,NY,OMEGA
         WRITE (MUNIT,12) UPPER,LOWER,LEFT,RIGHT
         IF (SUBLAT) WRITE (MUNIT,15) NXLL,NYLL,NXUR,NYUR
         WRITE (MUNIT,*) ' '
         WRITE (MUNIT,30)
         WRITE (MUNIT,40)
      END IF
      NLINES=4
      IF (SUBLAT) NLINES=NLINES+1
C      
5     FORMAT (' Output from example 6: Solutions to Laplace''s Equation'
     +        ' in two dimensions')
10    FORMAT (' NX =',I3,5X,' NY =',I3,5X,'  Omega = ',F6.3)
12    FORMAT (' field values on upper, lower, left, and right are ',
     +        4(2X,F5.2))
15    FORMAT (' Sublattice defined by  (',I3,',',I3,') and (',
     +        I3,',',I3,')')
30    FORMAT (8X,'Iter',11X,'Energy',13X,'Delta E',12X,'max Delta P')
40    FORMAT (8X,'----',11X,'------',13X,'-------',12X,'-----------')
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,NITER,ENERGY,DELE,DELMAX,NLINES)
C output text results to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:      
      INTEGER NITER               !number of iterations
      REAL DELMAX,OLDMAX          !keep track of changes in P
      REAL ENERGY,DELE            !current and delta energy
      INTEGER MUNIT               !fortran unit number
      INTEGER NLINES              !number of lines sent to terminal(I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     if screen is full, clear screen and retype headings          
      IF ((MOD(NLINES,TRMLIN-6) .EQ. 0) 
     +                          .AND. (MUNIT .EQ. OUNIT)) THEN
         CALL PAUSE('to continue...',1)
         CALL CLEAR
         WRITE (MUNIT,30)
         WRITE (MUNIT,40)
      END IF
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+1
C             
      WRITE (MUNIT,50) NITER,ENERGY,DELE,DELMAX
50    FORMAT (8X,I4,3(8X,1PE12.5))
30    FORMAT (8X,'Iter',11X,'Energy',13X,'Delta E',12X,'max Delta P')
40    FORMAT (8X,'----',11X,'------',13X,'-------',12X,'-----------')
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DISPLY(P,MUNIT)
C display field (P) as ascii characters: bound cond are small letters,
C non bound cond are in capitals; output sent to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E6'
      INCLUDE 'IO.ALL'
C Input variables:
      REAL P(NXMAX,NYMAX)                !field values
      INTEGER MUNIT                      !unit we're writing to
C Local variables:
      INTEGER IX,IY                      !lattice indices
      INTEGER PTEMP                      !field at current lattice site
      CHARACTER*1 FIELD(NXMAX)           !field as character data
      CHARACTER*80 BLNK                  !blanks for centering in X
      CHARACTER*1 ASKII(0:25),NEGASK(0:25)!charac data for display
      DATA BLNK /' '/
      DATA ASKII/'A','B','C','D','E','F','G','H','I','J','K','L','M',
     +           'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA NEGASK/'a','b','c','d','e','f','g','h','i','j','k','l','m',
     +           'n','o','p','q','r','s','t','u','v','w','x','y','z'/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) THEN
        DO 20 IY=1,YCNTR                   !center output
           WRITE (OUNIT,*) ' '
20      CONTINUE
      ELSE
        WRITE (MUNIT,*) ' Field values are:'          
      END IF
C
      DO 100 IY=NY,1,-1
         DO 50 IX=1,NX
            PTEMP=NINT(P(IX,IY))           !convert field value to ascii
            IF (PTEMP .GT. FMAX) PTEMP=FMAX   !keep things in bounds
            IF (PTEMP .LT. FMIN) PTEMP=FMIN
            IF (BNDCND(IX,IY) .EQ. 0) FIELD(IX)=ASKII(PTEMP)
            IF (BNDCND(IX,IY) .NE. 0) FIELD(IX)=NEGASK(PTEMP)
50       CONTINUE    
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
      SUBROUTINE GRFOUT(DEVICE,P,PGRF,MX,MY,E)
C display contours of the field P
C     the field values must be in an array that is exactly NX by NY;
C     this can be accomplished with implicit dimensioning which 
C     requires that PGRF and its dimensions be passed to this routine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'   
      INCLUDE 'PARAM.E6'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE                !which device
      REAL P(NXMAX,NYMAX)           !field values
      INTEGER MX,MY                 !NX and NY in disguise
      REAL PGRF(MX,MY)              !field values for graphing
      REAL E                        !energy
C Local variables:
      REAL BCX(5),BCY(5)            !corners of boundary conditions
      INTEGER I,IX,IY               !level index, lattice indices
      REAL PMAX,PMIN                !largest and smallest field value
      INTEGER IBC                   !index on boundary conditions
      REAL XSLOPE, YSLOPE           !slope of diagonal BC
      CHARACTER*9 CMIN,CMAX,CE      !data as characters
      INTEGER LMIN,LMAX,ELEN        !string lengths
      INTEGER SCREEN                !send to terminal
      INTEGER PAPER                 !make a hardcopy
      INTEGER FILE                  !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PMAX=0                        !set scale 
      PMIN=P(1,1)                                 
      DO 10 IX=1,MX
         DO 20 IY=1,MY
           IF (P(IX,IY) .GT. PMAX) PMAX=P(IX,IY)
           IF (P(IX,IY) .LT. PMIN) PMIN=P(IX,IY)
           PGRF(IX,IY)=P(IX,IY)     !load field into PGRF
20       CONTINUE
10    CONTINUE
C
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     calculate parameters for graphing
      NPLOT=1                       !how many plots
      IPLOT=1
C 
      YMAX=MY                              !axes data
      YMIN=1.
      XMIN=1.
      XMAX=MX
      Y0VAL=XMIN
      X0VAL=YMIN
      NXTICK=5
      NYTICK=5                             
C
      LABEL(1)='NX'                        !descriptions
      LABEL(2)='NY'
      CALL CONVRT(PMIN,CMIN,LMIN)
      CALL CONVRT(PMAX,CMAX,LMAX)
      CALL CONVRT(E,CE,ELEN)
      INFO='Pmin ='//CMIN(1:LMIN)//'  Pmax='//
     +         CMAX(1:LMAX)//'  E='//CE(1:ELEN)
      TITLE='Solutions to Laplace''s Equation'
C 
      CALL GTDEV(DEVICE)                   !device nomination
      IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
      CALL LNLNAX                          !draw axes
C 
C     display interior boundary conditions as dashed lines
      DO 200 IBC=1,NBC                     !loop over b.c.
         IF (BCGEOM(IBC) .EQ. DIAGNL) THEN  
           NPOINT=2
           IF ((XUR(IBC) .EQ. 1) .OR. (XUR(IBC) .EQ. 2)) YSLOPE=1
           IF ((XUR(IBC) .EQ. 3) .OR. (XUR(IBC) .EQ. 4)) YSLOPE=-1
           IF ((XUR(IBC) .EQ. 1) .OR. (XUR(IBC) .EQ. 4)) XSLOPE=1
           IF ((XUR(IBC) .EQ. 2) .OR. (XUR(IBC) .EQ. 3)) XSLOPE=-1
           BCX(1)=XLL(IBC)
           BCX(2)=XLL(IBC)+XSLOPE*YUR(IBC) 
           BCY(1)=YLL(IBC)
           BCY(2)=YLL(IBC)+YSLOPE*YUR(IBC)
         ELSE IF (BCGEOM(IBC) .EQ. RECTGL) THEN
           NPOINT=5
           BCX(1)=XLL(IBC)
           BCX(2)=BCX(1)
           BCX(3)=XUR(IBC)
           BCX(4)=BCX(3)
           BCY(1)=YLL(IBC)
           BCY(2)=YUR(IBC)
           BCY(3)=BCY(2)  
           BCY(4)=BCY(1)
           BCX(5)=BCX(1)
           BCY(5)=BCY(1)
         END IF
         ILINE=2
         CALL XYPLOT(BCX,BCY)          !plot boundaries 
200   CONTINUE
C
      ILINE=1
      CALL CONTOR(PGRF,MX,MY,PMIN,PMAX,NLEV)
C 
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)   !end graphics package
      IF (DEVICE .EQ. SCREEN) CALL TMODE         !switch to text mode
C 
100   FORMAT (/,' Patience, please; output going to a file.') 
C                                                         
      RETURN
      END                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE WRTOUT(P)
C write out field (P) and boundary conditions to GUNIT for reading back
C in as initial conditions or for graphing with an external package
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'MENU.ALL'
      INCLUDE 'PARAM.E6'
C Input variables:
      REAL P(NXMAX,NYMAX)           !field values
C Local variables:
      INTEGER IX,IY                 !lattice indices
      INTEGER NLINES                !lines written in PRMOUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL FLOPEN(GNAME,GUNIT)                  !open file
      MSTRNG(MINTS(IGNAME))=GNAME
      CALL PRMOUT(GUNIT,NLINES)                 !write out header
C
      DO 100 IX=1,NX                            !write out field
        WRITE (GUNIT,10) (P(IX,IY),IY=1,NY)
100   CONTINUE
10    FORMAT (5(2X,1PE14.7))
C
      DO 200 IX=1,NX                            !write out bound cond
        WRITE (GUNIT,20) (BNDCND(IX,IY),IY=1,NY)
200   CONTINUE
20    FORMAT (40(1X,I1))
      CALL FLCLOS(GNAME,GUNIT)                  !close it up
C
      RETURN
      END
