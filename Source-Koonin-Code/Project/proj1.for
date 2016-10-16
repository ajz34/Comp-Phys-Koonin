CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM PROJ1
C     Project 1: Scattering by a central 6-12 potential
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop; execute once for each set of param
         CALL PARAM       !get input parameters
         CALL ARCHON      !calculate scattering angles
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON  
C calculate scattering angles as a function of the impact parameter
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables
      INCLUDE 'PARAM.P1'
      INCLUDE 'IO.ALL'
C Local variables:
      REAL RMIN                 !radius of closest approach
      REAL ANGLE(0:MAX)         !final scattering angle
      REAL IMPACT(0:MAX)        !impact parameter
      REAL INT1, INT2           !values of the two integrals
      INTEGER IB                !index of impact param
      INTEGER NLINES            !number of lines printed so far
      INTEGER SCREEN            !send to terminal
      INTEGER PAPER             !make a hardcopy
      INTEGER FILE              !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     output summary of parameters
      IF (TTERM) CALL PRMOUT(OUNIT,NLINES)
      IF (TFILE) CALL PRMOUT(TUNIT,NLINES)
      IF (GFILE) CALL PRMOUT(GUNIT,NLINES)
C
      DO 25 IB=0,(NB-1)                      !loop over impact param
C
        IMPACT(IB)=BMIN+(IB*DB)              !current value of B
        CALL THETA(IMPACT(IB),INT1,INT2,RMIN)!calculate THETA for each B
        ANGLE(IB)=(2*IMPACT(IB)*(INT1-INT2)*180)/PI !convert to degrees
C
C       text output                           
        IF (TTERM) CALL TXTOUT(OUNIT,NLINES,IMPACT(IB),RMIN,ANGLE(IB))
        IF (TFILE) CALL TXTOUT(TUNIT,NLINES,IMPACT(IB),RMIN,ANGLE(IB))
C
   25 CONTINUE
C
      IF (TTERM) CALL PAUSE('to continue...',1)     
      IF (TTERM) CALL CLEAR
C
      IF (GTERM) CALL GRFOUT(SCREEN,ANGLE,IMPACT)   !graphics output
      IF (GFILE) CALL GRFOUT(FILE,ANGLE,IMPACT)
      IF (GHRDCP) CALL GRFOUT(PAPER,ANGLE,IMPACT)
C      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE THETA(B,INT1,INT2,RMIN)      
C  searches to find RMIN (closest approach)   
C  integrates FNI1 from B    to RMAX to find INT1
C  integrates FN12 from RMIN to RMAX to find INT2 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P1'
C Input/Output variables:
      REAL B                     !impact parameter (input)
      REAL INT1, INT2            !values for the two integrals (output)
      REAL RMIN                  !radius of closes approach (output)
C Local variables:
      REAL R                     !current radius
      REAL DR                    !current search step size
      REAL U                     !U=SQRT(R-B) or SQRT(R-RMIN)
      REAL UMAX                  !UMAX=SQRT(RMAX-B) or SQRT(RMAX-RMIN)
      REAL H                     !radial step for quadrature
      REAL SUM                   !sum for quadrature 
      INTEGER IU                 !index on U
      REAL FNV                   !scattering potential (function)
      REAL FNI1,FNI2             !first and second integrands (function)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     define functions
      FNV(R)      = 4*(R**(-12)-R**(-6))
      FNI1(B,R)   = 1/R**2/SQRT(1-(B/R)**2)
      FNI2(B,R,E) = 1/R**2/SQRT(1-(B/R)**2-FNV(R)/E)
C
      RMIN=RMAX                      !inward search for turning point
      DR = 0.2
1     IF (DR.GT.TOLR) THEN                                      
         RMIN=RMIN-DR                !simple search
         IF ((1-((B/RMIN)**2)-(FNV(RMIN)/E)).LT.0.0) THEN
           RMIN=RMIN +DR             !RMIN is where FNI2 is infinite
           DR=DR/2                   !i.e., the denominator=0
         END IF
      GO TO 1
      END IF
C
      SUM=0                          !first integral by rectangle rule
      UMAX=SQRT(RMAX-B)              !change variable to U=SQRT(R-B)
      H=UMAX/NPTS                    !to remove singularity at R=B
      DO 10 IU=1,NPTS
            U=H*(IU-0.5)
            R=U**2+B                 !change back to R to eval integrand
            SUM=SUM+U*FNI1(B,R)
   10 CONTINUE
      INT1=2*H*SUM
C
      SUM=0                          !second integral by rectangle rule
      UMAX=SQRT(RMAX-RMIN)           !change of variable to 
      H=UMAX/NPTS                    !to remove singularity at R=RMIN
      DO 22 IU=1,NPTS
            U=H*(IU-0.5)
            R=U**2+RMIN              !change back to R to eval integrand
            SUM=SUM+U*FNI2(B,R,E)
   22 CONTINUE
      INT2=2*H*SUM
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
      INCLUDE 'PARAM.P1'
C Local parameters:
      CHARACTER*80 DESCRP        !program description
      DIMENSION DESCRP(22)
      INTEGER NHEAD,NTEXT,NGRF   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL SETUP                !get environment parameters
C 
      !display header screen     
      DESCRP(1)= 'PROJECT 1'
      DESCRP(2)= 'Scattering by a 6-12 potential'
      NHEAD=2
C 
C     text output description
      DESCRP(3)= 'impact parameter, radius of closest approach, '
      DESCRP(4)= 'and the angle of deflection'
      NTEXT=2
C 
C     graphics output description
      DESCRP(5)= 'deflection function'
     +   //' (scattering angle vs. impact parameter)'
      NGRF=1
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRF)
C
C     calculate constants                
      PI=4*ATAN(1.0)
C 
      CALL MENU                          !setup constant part of menu
C                
      MTYPE(13)=FLOAT
      MPRMPT(13)= 'Enter incident energy E (scaled units)'
      MTAG(13)= 'Incident energy (scaled units)'
      MLOLIM(13)=.01
      MHILIM(13)=1000.
      MREALS(13)=1.0
C                
      MTYPE(14)=FLOAT
      MPRMPT(14)= 'Enter the minimum impact parameter (scaled units)'
      MTAG(14)= 'Minimum impact parameter Bmin (scaled units)'
      MLOLIM(14)=.0
      MHILIM(14)=2.5
      MREALS(14)=0.1
C 
      MTYPE(15)=FLOAT
      MPRMPT(15)= 'Enter the maximum impact parameter (scaled units)'
      MTAG(15)= 'Maximum impact parameter Bmax (scaled units)'
      MLOLIM(15)=.0
      MHILIM(15)=7.
      MREALS(15)=2.4
C 
      MTYPE(16)=NUM
      MPRMPT(16)= 'Enter number of values for impact parameter'
      MTAG(16)= 'Number of values for impact parameter'
      MLOLIM(16)=1.
      MHILIM(16)=MAX
      MINTS(16)=20
C 
      MTYPE(17)=SKIP
      MREALS(17)=35.
C                    
      MTYPE(38)=FLOAT
      MPRMPT(38)= 'Enter tolerance for turning point search'
      MTAG(38)= 'Turning point search tolerance'
      MLOLIM(38)=.000005
      MHILIM(38)=.1
      MREALS(38)=.0001
C                   
      MTYPE(39)=FLOAT
      MPRMPT(39)= 'Enter radius at which V(r)<<E (scaled units)'
      MTAG(39)= 'Maximum radius RMAX (scaled units)'
      MLOLIM(39)=0.
      MHILIM(39)=10.0
      MREALS(39)=2.5           
C 
      MTYPE(40)=NUM
      MPRMPT(40)= 'Enter number of quadrature points'
      MTAG(40)= 'Number of quadrature points'
      MLOLIM(40)=20.
      MHILIM(40)=300.
      MINTS(40)=40
C 
      MTYPE(41)=SKIP
      MREALS(41)=60.
C 
      MSTRNG(MINTS(75))= 'proj1.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C 
      MSTRNG(MINTS(86))= 'proj1.grf'
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
      INCLUDE 'PARAM.P1'
      INCLUDE 'MAP.P1'
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
      E=MREALS(IE)
      BMIN=MREALS(IBMIN)
      BMAX=MREALS(IBMAX)
      NB=MINTS(INB)
      RMAX=MREALS(IRMAX)
      TOLR=MREALS(ITOLR)
      NPTS=MINTS(INPTS)
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
C     files may have been renamed
      MSTRNG(MINTS(ITNAME))=TNAME
      MSTRNG(MINTS(IGNAME))=GNAME
C
      CALL PCHECK            !check 0<BMIN<BMAX<=RMAX
      DB=(BMAX-BMIN)/(NB-1)  !calculate step in B                
      CALL CLEAR
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PCHECK
C to ensure that 0<BMIN<BMAX<=RMAX
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:                      
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P1'
      INCLUDE 'MAP.P1'
C Function:
      REAL GETFLT                      !gets floating point number 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
40    IF (BMAX .LE. BMIN) THEN
         WRITE (OUNIT,10) BMAX,BMIN
         WRITE (OUNIT,20)
         MREALS(IBMIN)=GETFLT(MREALS(IBMIN),MLOLIM(IBMIN),MHILIM(IBMIN),
     +                 'Enter BMIN')                           
         MREALS(IBMAX)=GETFLT(MREALS(IBMAX),BMIN,MHILIM(IBMAX),
     +                 'Enter BMAX')
         BMAX=MREALS(IBMAX)
         BMIN=MREALS(IBMIN)
      GO TO 40
      END IF
C
50    IF (RMAX .LE. BMAX) THEN
         WRITE (OUNIT,30) RMAX,BMAX
         WRITE (OUNIT,20)
         !prompt for new values
         MREALS(IBMAX)=GETFLT(MREALS(IBMAX),BMIN,MHILIM(IBMAX),
     +                 'Enter BMAX')
         BMAX=MREALS(IBMAX)
         MREALS(IRMAX)=GETFLT(MREALS(IRMAX),MLOLIM(IRMAX),MHILIM(IRMAX),
     +                 'Enter RMAX')
         RMAX=MREALS(IRMAX)
      GO TO 50
      END IF    
C      
10    FORMAT (/,' BMAX (=',F6.3,')',' is less than BMIN (=',
     +        F6.3,'),')
20    FORMAT ('  but it should be the other way round; '
     +         ' reenter parameter values')
30    FORMAT (/,' RMAX (=',F6.3,')',' is less than BMAX (=',F6.3,')')
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRMOUT(MUNIT,NLINES)
C outputs parameter summary to the specified unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P1'
C Passed variables
      INTEGER MUNIT             !unit number for output (input)
      INTEGER NLINES            !number of lines written so far (output)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) CALL CLEAR
C
      WRITE (MUNIT,19)
      WRITE (MUNIT,21)
      WRITE (MUNIT,27)  E
      WRITE (MUNIT,23)  RMAX,NPTS
      WRITE (MUNIT,25)  TOLR
      WRITE (MUNIT,29)  BMIN,BMAX,NB
      WRITE (MUNIT,30)
      WRITE (MUNIT,19)
C     
C     different data headers for graphics and text output
      IF (MUNIT .NE. GUNIT) THEN
        WRITE (MUNIT,31)
        WRITE (MUNIT,32)
        WRITE (MUNIT,33)
      ELSE
        WRITE (MUNIT,35)
      END IF
C
      NLINES = 11
C
19    FORMAT  (' ')
21    FORMAT  (' Output from project 1:',
     +         '  Scattering by a central 6-12 potential ')
23    FORMAT  (' Rmax=',F6.3,5X,'number of quadrature points =',I4)
25    FORMAT  (' Turning point tolerance = ',E12.5)
27    FORMAT  (' Incident energy = 'F8.3)
29    FORMAT  (' Bmin = ',F6.3,5X,'Bmax = ',F6.3,
     +         5X,'number of B values = ',I3)
30    FORMAT  (' Lengths are in scaled units, angles in degrees')
31    FORMAT (15X,'Impact',14X,'Closest',16X,'Scattering')
32    FORMAT (13X,'Parameter',13X,'Approach'17X,'Angle')
33    FORMAT (13X,'---------',13X,'--------'17X,'-----')
35    FORMAT  (2X,'Impact parameter',2X,'Scattering Angle')
C            
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,NLINES,B,RMIN,ANGLE)
C writes results for one impact parameter to the requested unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Passed variables:
      INTEGER MUNIT                !output unit specifier (input)
      INTEGER NLINES               !number of lines printed so far (I/0)
      REAL B                       !impact parameter (input)
      REAL RMIN                    !radius of closest approach (input)
      REAL ANGLE                   !scattering angle(input)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     if screen is full, clear screen and retype headings          
      IF ((MOD(NLINES,TRMLIN-4) .EQ. 0) 
     +                          .AND. (MUNIT .EQ. OUNIT)) THEN
         CALL PAUSE('to continue...',1)
         CALL CLEAR
         WRITE (MUNIT,31)
         WRITE (MUNIT,32)
         WRITE (MUNIT,33)
         NLINES=3
      END IF
C
      WRITE (MUNIT,35) B,RMIN,ANGLE
C
C     keep track of printed lines only for terminal output
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+1
C
31    FORMAT (15X,'Impact',14X,'Closest',16X,'Scattering')
32    FORMAT (13X,'Parameter',13X,'Approach'17X,'Angle')
33    FORMAT (13X,'---------',13X,'--------'17X,'-----')
35    FORMAT  (15X,F6.3,15X,F6.3,15X,F9.3)
C
      RETURN     
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFOUT(DEVICE,ANGLE,IMPACT)
C graphs scattering angle vs. impact parameter
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables                                                   
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P1'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE            !which device is being used?
      REAL ANGLE(0:MAX)         !final scattering angle
      REAL IMPACT(0:MAX)        !impact parameter
C Local variables
      INTEGER IB                !index on impact parameter
      CHARACTER*9 CE,CR         !energy, RMAX as a character string
      INTEGER LEN,RLEN          !length of string
      INTEGER SCREEN            !send to terminal
      INTEGER PAPER             !make a hardcopy
      INTEGER FILE              !send to a file
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
          YMIN=0.
          DO 20 IB=0,NB-1
             IF (ANGLE(IB) .GT. YMAX) YMAX=ANGLE(IB)
             IF (ANGLE(IB) .LT. YMIN) YMIN=ANGLE(IB)
20        CONTINUE
          XMIN=BMIN
          XMAX=BMAX
          Y0VAL=BMIN
          X0VAL=0.
C 
          NPOINT=NB
C 
          ILINE=1                      !line and symbol styles
          ISYM=1
          IFREQ=1
          NXTICK=5
          NYTICK=5
C
          CALL CONVRT(E,CE,LEN)             !titles and labels
          CALL CONVRT(RMAX,CR,RLEN)
          INFO='RMAX = '//CR(1:RLEN)
          TITLE = 'Scattering from 6-12 potential, Energy='//CE
          LABEL(1)= 'Impact parameter (scaled units)'
          LABEL(2)= 'Scattering angle (degrees)'
C 
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
          CALL LNLNAX                          !draw axes
      END IF
C                                                      
C     output results
      IF (DEVICE .EQ. FILE) THEN
          WRITE (GUNIT,70) (IMPACT(IB),ANGLE(IB),IB=1,NB)
      ELSE
          CALL XYPLOT (IMPACT,ANGLE)
      END IF
C                                                 
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)  !close graphing package
      IF (DEVICE .EQ. SCREEN) CALL TMODE        !switch to text mode
C
70    FORMAT (2(5X,E11.3))
100   FORMAT (/,'Patience, please; output going to a file.')
C
      RETURN                                                    
      END
                       
