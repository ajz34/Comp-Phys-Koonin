CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C file GRAPHIT.HI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LNLNAX
C draws linear-linear axes 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables
       INCLUDE 'GRFDAT.ALL'
       REAL XSIZE,YSIZE              !size of axes in inches
C Local variables                                                       
       REAL XSTP,YSTP                !length between ticks in user units
       REAL VSIZE,HSIZE              !XSIZE and YSIZE
       COMMON/SIZE/VSIZE,HSIZE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      write title
       IF (IPLOT .EQ. 1) CALL BANNER(TITLE)
C      set subplot area
       CALL SUBPLT(IPLOT,NPLOT,XSIZE,YSIZE)
C      label axes
       CALL LABELS(LABEL)
C      write informational message
       HSIZE=XSIZE
       VSIZE=YSIZE
       CALL LEGEND
C      draw linear-linear axes
       XSTP=(XMAX-XMIN)/NXTICK
       YSTP=(YMAX-YMIN)/NYTICK
       CALL GRAF(XMIN,XSTP,XMAX,YMIN,YSTP,YMAX)
C            
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LGLNAX
C draws log-linear axes 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables
       INCLUDE 'GRFDAT.ALL'
       REAL XSIZE,YSIZE           !size of axes in inches
C Local variables
       REAL YCYCLE                !inches per cycle on y axis 
       REAL XSTEP                 !length in user units/length in inches
       REAL VSIZE,HSIZE              !XSIZE and YSIZE
       COMMON/SIZE/VSIZE,HSIZE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      write title
       IF (IPLOT .EQ. 1) CALL BANNER(TITLE)
C      set subplot area
       CALL SUBPLT(IPLOT,NPLOT,XSIZE,YSIZE)
C      label axes
       CALL LABELS(LABEL)
C      write informational message
       HSIZE=XSIZE
       VSIZE=YSIZE
       CALL LEGEND
C      draw log-linear axes
       YCYCLE=YSIZE/(ALOG10(YMAX)-ALOG10(YMIN))
       XSTEP=(XMAX-XMIN)/XSIZE
       CALL YLOG(XMIN,XSTEP,YMIN,YCYCLE)
C
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LGLGAX
C draws log-linear axes 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables
       INCLUDE 'GRFDAT.ALL'
C Local variables
       REAL XSIZE,YSIZE              !size of axes in inches
       REAL XCYCLE,YCYCLE            !inches per cycle
       REAL VSIZE,HSIZE              !XSIZE and YSIZE
       COMMON/SIZE/VSIZE,HSIZE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      write title                            
       IF (IPLOT .EQ. 1) CALL BANNER(TITLE)
C      set subplot area
       CALL SUBPLT(IPLOT,NPLOT,XSIZE,YSIZE)
C      label axes
       CALL LABELS(LABEL)
C      write informational message
       HSIZE=XSIZE
       VSIZE=YSIZE
       CALL LEGEND
C      draw log-log axes
       XCYCLE=XSIZE/(ALOG10(XMAX)-ALOG10(XMIN))
       YCYCLE=YSIZE/(ALOG10(YMAX)-ALOG10(YMIN))
       CALL LOGLOG(XMIN,XCYCLE,YMIN,YCYCLE)
C
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE LEGEND
C write information at the top left of the plot
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'GRFDAT.ALL'                                    
C Passed variables:
       REAL XSIZE,YSIZE             
C Local variables:
       INTEGER LEN,LENTRU           !length of char string
       REAL VSIZE,HSIZE             !XSIZE and YSIZE
       COMMON/SIZE/VSIZE,HSIZE
       DATA INFO/' '/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      prints message at top of plotting area
       LEN=LENTRU(INFO)
       IF (LEN .GT. 0) THEN
          CALL MESSAG(INFO(1:LEN),LEN,HSIZE*.05,VSIZE*1.03)
       END IF
       RETURN 
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE XYPLOT(X,Y)
C plots xy data
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'GRFDAT.ALL'
       REAL X                        !independent variable data array
       REAL Y                        !dependent variable data array
       DIMENSION X(NPOINT),Y(NPOINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      set line type
       IF (ILINE .EQ. 2) CALL DOT
       IF (ILINE .EQ. 3) CALL DASH
       IF (ILINE .EQ. 4) CALL CHNDSH
       IF (ILINE .EQ. 5) IFREQ=-IFREQ    !no line at all
C
C      set symbol type
       IF (ISYM .EQ. 1) CALL MARKER(16)  !circle
       IF (ISYM .EQ. 2) CALL MARKER(2)   !triangle
       IF (ISYM .EQ. 3) CALL MARKER(0)   !square
       IF (ISYM .EQ. 4) CALL MARKER(4)   !cross
       IF (ISYM .EQ. 5) THEN             !point
           CALL SCLPIC(.3)
           CALL MARKER(15)
       END IF
C
C      plot
       CALL CURVE(X,Y,NPOINT,IFREQ)
C
C      reset line type and frequency
       IF (ILINE .EQ. 2) CALL RESET('DOT')
       IF (ILINE .EQ. 3) CALL RESET('DASH')
       IF (ILINE .EQ. 4) CALL RESET('CHNDSH')
       IF (ILINE .EQ. 5) IFREQ=-IFREQ
C
       RETURN
       END    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE CONTOR(Z,MX,MY,ZMIN,ZMAX,NCONT)
C drawn NCONT lines for data contained in Z
C lines alternate between solid and dashed; 
C if possible, solid lines are labeled
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'GRFDAT.ALL'
C Passed variables:
       INTEGER MX,MY              !dimensions of Z
       REAL Z(MX,MY)              !data
       REAL ZMIN,ZMAX             !limits on data
       INTEGER NCONT              !number of contours to draw
C Local variables:
       REAL ZINCR                 !increment between contours
       REAL WORK(4000)            !work space for DISSPLA
       COMMON WORK
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       CALL FRAME                 !draw a box around area
       CALL BCOMON(4000)          !pass information about workspace
       ZINCR=(ZMAX-ZMIN)/NCONT    !increments between contours
       CALL CONMAK(Z,MX,MY,ZINCR) !make the contour lines
C
C      provide two line types
       CALL CONLIN(0,'SOLID','LABELS',1,10)  
       CALL CONLIN(1,'DOT','NOLABELS',1,10)  
       CALL CONTUR(2,'LABELS','DRAW')          !draw contour lines
C
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE GTDEV(DEVICE)
C sets device for graphics output to DEVICE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Passed variables:
      INTEGER DEVICE                    !device flag
      INTEGER SCREEN                    !send to terminal
      INTEGER PAPER                     !make a hardcopy
      INTEGER FILE                      !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      reset output device
       CALL IOMGR(0,-102)
C      4014 tektronix screen at 9600 baud
       IF (DEVICE .EQ. SCREEN) CALL TEKALL(4014,960,0,0,0) 
C      postscript Printer with default values
       IF (DEVICE .EQ. PAPER) CALL PSCRPT(0,0,0) 
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE GPAGE(DEVICE)
C end graphics page
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       INTEGER DEVICE                !which device is it?
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       CALL ENDPL(0)
       CALL CLEAR
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The following subroutines are not called directly by the user
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE BANNER(TITLE)
C prints title to top of graphics page
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       CHARACTER*60 TITLE             !title to be printed
C Functions:
       INTEGER LENTRU                 !returns string length
C Local variables:
       INTEGER LENTTL                 !length of title                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      treat title as a graph by itself so that it is centered
       CALL AREA2D(7.5,9.)
       CALL HEIGHT(.14)
       LENTTL=LENTRU(TITLE)      
       CALL HEADIN(TITLE(1:LENTTL),LENTTL,1.,1)
       CALL ENDGR(0)
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE LABELS(LABEL)
C labels both x and y axes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       CHARACTER*60 LABEL(2)         !x and y labels
C Functions:
       INTEGER LENTRU                !returns string length
C Local variables:
       INTEGER LENX,LENY             !length of labels
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       LENX=LENTRU(LABEL(1))
       LENY=LENTRU(LABEL(2))
       CALL XNAME(LABEL(1)(1:LENX),LENX)
       CALL YNAME(LABEL(2)(1:LENY),LENY)
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE SUBPLT(IPLOT,NPLOT,XSIZE,YSIZE)
C defines subplotting area 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The choice of subplot size and placement is dependent on the
C device and graphics package; those given here are for an 8.5 X 11
C inch page using DISSPLA        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
       INTEGER NPLOT             !total number of plots on page
       INTEGER IPLOT             !number of current plot         
       REAL XSIZE,YSIZE          !size of plotting area in inches
C Local variables:               
       REAL XORIG,YORIG          !location of lower left corner (inches)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      finish last plot (sets DISSPLA back to level 1)
       IF (IPLOT .NE. 1) CALL ENDGR(0)
C
C      define subplotting area (in inches); whole page is 8.5 x 11
       IF (NPLOT .EQ. 1) THEN
           XSIZE=7.5
           YSIZE=9.
           CALL HEIGHT(.14)
       ELSE IF (NPLOT .EQ. 2) THEN
           YSIZE=4.25
           XSIZE=6.5
           XORIG=1.
           IF (IPLOT .EQ. 1) YORIG=5.75
           IF (IPLOT .EQ. 2) YORIG=.75
           CALL HEIGHT(.125)
       ELSE IF (NPLOT .EQ. 3) THEN
           YSIZE=2.80
           XSIZE=6.5
           XORIG=1.
           IF (IPLOT .EQ. 1) YORIG=7.20
           IF (IPLOT .EQ. 2) YORIG=3.9
           IF (IPLOT .EQ. 3) YORIG=.6
           CALL HEIGHT(.125)
       ELSE IF (NPLOT .EQ. 4)THEN
           XSIZE=3.25
           YSIZE=4.25
           CALL HEIGHT(.11)
           IF (IPLOT .EQ. 1) THEN
               XORIG=.5
               YORIG=5.75
           ELSE IF (IPLOT .EQ. 2) THEN
               XORIG=4.75
               YORIG=5.75
           ELSE IF (IPLOT .EQ. 3) THEN
               XORIG=.5
               YORIG=.5
           ELSE IF (IPLOT .EQ. 4) THEN
               XORIG=4.75
               YORIG=.5
           END IF
        END IF
C 
C      use default origin if there is only one plot
       IF (NPLOT .NE. 1) CALL PHYSOR(XORIG,YORIG)
C
       CALL AREA2D(XSIZE,YSIZE)
C
       RETURN
       END
