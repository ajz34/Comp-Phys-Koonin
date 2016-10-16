CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  file UTIL.for
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C displays header and description of output to screen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:                          
      INCLUDE 'IO.ALL'
C Input varaibles:
      CHARACTER*(*) DESCRP(20)    !description of program and output
      INTEGER NHEAD,NTEXT,NGRAPH  !number of lines for each description
C Local variables:
      INTEGER N                   !current line number
      INTEGER LENGTH              !true length of character strings
      INTEGER NBLNKS              !num of blanks needed to center string
      CHARACTER*80 BLANKS         !array of blanks for centering
C Function:
      INTEGER LENTRU              !true length of character string 
      DATA BLANKS/' '/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL CLEAR                  !vertically center output
      DO 190 N=1,(TRMLIN-18-NHEAD-NGRAPH-NTEXT)/2
         WRITE (OUNIT,20)
190   CONTINUE
C
C     write out constant part of header
      WRITE (OUNIT,40)
      WRITE (OUNIT,50)
      WRITE (OUNIT,60)
      WRITE (OUNIT,80)
      WRITE (OUNIT,20)
      WRITE (OUNIT,20)
C
C     write out chapter dependent section of the header
      DO 140 N=1,NHEAD+NTEXT+NGRAPH
         IF (N .EQ. NHEAD+1) THEN !text output header
             WRITE (OUNIT,110) 
         END IF
         IF (N .EQ. NHEAD+NTEXT+1) THEN
             WRITE (OUNIT,115)    !graphics output header
         END IF
         LENGTH=LENTRU(DESCRP(N)) !horizontally center output
         NBLNKS=(80-LENGTH)/2
         WRITE (OUNIT,120) BLANKS(1:NBLNKS),DESCRP(N)(1:LENGTH)
140   CONTINUE
C
      CALL PAUSE('to begin the program...',1)
      CALL CLEAR
C
20    FORMAT (' ')
40    FORMAT (/,30X,'COMPUTATIONAL PHYSICS')
50    FORMAT (/,32X,'(FORTRAN VERSION)')
60    FORMAT (/,20X,'by Steven E. Koonin and Dawn C. Meredith')
80    FORMAT (/,14X,
     +   'Copyright 1989, Benjamin/Cummings Publishing Company')
110   FORMAT (/,30X, 'Text output displays')
115   FORMAT (/,28X, 'Graphics output displays')
120   FORMAT (A,A)
C
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MENU
C sets up the part of the menu that is the same for all programs
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     main menu
      MTYPE(1)=CLRTRM
C             
      MTYPE(2)=MTITLE
      MPRMPT(2)='MAIN MENU'
      MLOLIM(2)=2
      MHILIM(2)=1
C 
      MTYPE(3)=MTITLE
      MPRMPT(3)='1) Change physical parameters'
      MLOLIM(3)=0
      MHILIM(3)=0
C 
      MTYPE(4)=MTITLE
      MPRMPT(4)='2) Change numerical parameters'
      MLOLIM(4)=0
      MHILIM(4)=0
C               
      MTYPE(5)=MTITLE
      MPRMPT(5)='3) Change output parameters'
      MLOLIM(5)=0
      MHILIM(5)=0
C               
      MTYPE(6)=MTITLE
      MPRMPT(6)='4) Display physical and numerical parameters'
      MLOLIM(6)=0.
      MHILIM(6)=0.
C             
      MTYPE(7)=MTITLE
      MPRMPT(7)='5) Display output parameters'
      MLOLIM(7)=0.
      MHILIM(7)=0.
C                             
      MTYPE(8)=MTITLE                         
      MPRMPT(8)='6) Run program'
      MLOLIM(8)=0
      MHILIM(8)=0
C             
      MTYPE(9)=MTITLE
      MPRMPT(9)='7) Stop program'
      MLOLIM(9)=0
      MHILIM(9)=1
C 
      MTYPE(10)=MCHOIC
      MPRMPT(10)= 'Make a menu choice'
      MTAG(10)='11 36 61 91 94 99 99'
      MLOLIM(10)=1
      MHILIM(10)=7
      MINTS(10)=6
      MREALS(10)=-6
C 
C     physical parameters
      MTYPE(11)=CLRTRM
C              
      MTYPE(12)=TITLE
      MPRMPT(12)= 'PHYSICAL PARAMETERS'
      MLOLIM(12)=2.
      MHILIM(12)=1.
C                
      MTYPE(35)=SKIP
      MREALS(35)=1.
C 
C     numerical parameters
      MTYPE(36)=CLRTRM
C              
      MTYPE(37)=TITLE
      MPRMPT(37)= 'NUMERICAL PARAMETERS'
      MLOLIM(37)=2.
      MHILIM(37)=1.
C                    
      MTYPE(60)=SKIP
      MREALS(60)=1.
C 
C     output menu
      MTYPE(61)=CLRTRM
C              
      MTYPE(62)=MTITLE
      MPRMPT(62)= 'OUTPUT MENU'
      MLOLIM(62)=0.
      MHILIM(62)=1.
C                
      MTYPE(63)=MTITLE
      MPRMPT(63)='1) Change text output parameters'
      MLOLIM(63)=0.
      MHILIM(63)=0.
C 
      MTYPE(64)=MTITLE
      MPRMPT(64)='2) Change graphics output parameters'
      MLOLIM(64)=0.
      MHILIM(64)=0.
C                
      MTYPE(65)=MTITLE
      MPRMPT(65)='3) Return to main menu'
      MLOLIM(65)=0.
      MHILIM(65)=1.   
C               
      MTYPE(66)=MCHOIC         
      MPRMPT(66)= 'Make menu choice and press Return'
      MTAG(66)='71 81 01'
      MLOLIM(66)=1.
      MHILIM(66)=3.
      MINTS(66)=3.
C 
C     text output parameters
      MTYPE(71)=CLRTRM
C              
      MTYPE(72)=TITLE
      MPRMPT(72)= 'TEXT OUTPUT PARAMETERS'
      MLOLIM(72)=2.
      MHILIM(72)=1.
C                   
      MTYPE(73)=BOOLEN
      MPRMPT(73)= 'Do you want text output displayed on screen?'
      MTAG(73)= 'Text output to screen'
      MINTS(73)=TXTTRM
C              
      MTYPE(74)=NOSKIP
      MPRMPT(74)= 'Do you want text output sent to a file?'
      MTAG(74)= 'Text output to file'
      MREALS(74)=76.
      MINTS(74)=TXTFIL
C            
      MTYPE(75)=CHSTR
      MPRMPT(75)= 'Enter name of file for text output'
      MTAG(75)= 'File name for text output'
      MLOLIM(75)=1.
      MHILIM(75)=12.
      MINTS(75)=1
      MSTRNG(MINTS(75))= 'cmphys.txt'
C                      
      MTYPE(80)=SKIP
      MREALS(80)=61.
C 
C     graphics output parameters
      MTYPE(81)=CLRTRM
C              
      MTYPE(82)=TITLE
      MPRMPT(82)= 'GRAPHICS OUTPUT PARAMETERS'
      MLOLIM(82)=2.
      MHILIM(82)=1.
C 
      MTYPE(83)=BOOLEN
      MPRMPT(83)= 'Do you want graphics sent to the terminal?'
      MTAG(83)= 'Graphics output to terminal'
      MINTS(83)=GRFTRM
C 
      MTYPE(84)=BOOLEN
      MPRMPT(84)= 'Do you want graphics sent to the hardcopy device?'
      MTAG(84)= 'Graphics output to hardcopy device'
      MINTS(84)=GRFHRD
C                                                    
      MTYPE(85)=NOSKIP
      MPRMPT(85)= 'Do you want data for graphing sent to a file?'
      MTAG(85)= 'Data for graphing sent to file'
      MREALS(85)=87.
      MINTS(85)=GRFFIL
C               
      MTYPE(86)=CHSTR
      MPRMPT(86)= 'Enter name of file for graphics data'
      MTAG(86)= 'File for graphics data'
      MLOLIM(86)=2.
      MHILIM(86)=12.
      MINTS(86)=2.
      MSTRNG(MINTS(86))= 'cmphys.grf'
C 
      MTYPE(90)=SKIP
      MREALS(90)=61.
C                
C     printing numerical and physical parameters
      MTYPE(91)=PPRINT
      MLOLIM(91)=11.
      MHILIM(91)=60.
C 
      MTYPE(92)=SKIP
      MREALS(92)=1.
C            
C     printing output parameters
      MTYPE(94)=PPRINT
      MLOLIM(94)=71.
      MHILIM(94)=90.
C 
      MTYPE(95)=SKIP
      MREALS(95)=1.
C
      RETURN
      END         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ASK(START,END)
C executes menu items from START to END;
C see Appendix A for a description of the menu
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'        
      INCLUDE 'IO.ALL'          
C Input variables:        
      INTEGER START,END         !starting/ending menu items to execute
C Local variables:
      INTEGER I                 !current menu item
      INTEGER ILOW,IHIGH        !integer limits for NUM type
      INTEGER NUMSKP            !number of blank lines to print
      INTEGER ICHOIC            !current menu choice
C Functions
      CHARACTER*40 CHARAC       !character input
      REAL GETFLT               !real input
      INTEGER GETINT            !integer input
      INTEGER PARSE             !determines menu branching
      INTEGER YESNO             !boolean input
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      I=START
1000  CONTINUE
         IF (MTYPE(I) .EQ. FLOAT) THEN
            MREALS(I) = 
     +         GETFLT(MREALS(I),MLOLIM(I),MHILIM(I),MPRMPT(I))
C                    
         ELSE IF ( MTYPE(I) .EQ. NUM) THEN
            ILOW = MLOLIM(I)
            IHIGH = MHILIM(I)
            MINTS(I) = GETINT(MINTS(I), ILOW, IHIGH, MPRMPT(I))
C 
         ELSE IF (MTYPE(I) .EQ. BOOLEN) THEN
            MINTS(I) = YESNO(MINTS(I), MPRMPT(I))
C 
         ELSE IF (MTYPE(I) .EQ. CHSTR) THEN                           
            MSTRNG(MINTS(I)) = 
     +        CHARAC(MSTRNG(MINTS(I)),INT(MHILIM(I)),MPRMPT(I))
C 
         ELSE IF ( MTYPE(I) .EQ. MCHOIC) THEN
            ILOW = MLOLIM(I)
            IHIGH = MHILIM(I)
            ICHOIC = GETINT(MINTS(I), ILOW, IHIGH, MPRMPT(I))
C           if MREALS is > 0, save ICHOIC and change default
            IF (MREALS(I) .GT. 0) THEN
                MREALS(I)=REAL(ICHOIC)
                MINTS(I)=ICHOIC
C           if MREALS is < 0, save ICHOIC but leave default the same
            ELSE IF (MREALS(I) .LT. 0) THEN
                MREALS(I)=-REAL(ICHOIC)
            END IF
            I = PARSE (MTAG(I), ICHOIC) - 1
C 
         ELSE IF (MTYPE(I) .EQ. TITLE .OR. MTYPE(I) .EQ. MTITLE) THEN
            NUMSKP = MLOLIM(I)
            CALL PRBLKS(NUMSKP)
            WRITE (OUNIT, 10) MPRMPT(I)
            NUMSKP = MHILIM(I)
            CALL PRBLKS(NUMSKP)
C             
         ELSE IF (MTYPE(I) .EQ. YESKIP) THEN
            MINTS(I) = YESNO(MINTS(I), MPRMPT(I))
            IF (MINTS(I) .NE. 0) THEN
               I = MREALS(I) - 1
            END IF
C 
         ELSE IF (MTYPE(I) .EQ. NOSKIP) THEN
            MINTS(I) = YESNO(MINTS(I), MPRMPT(I))
            IF (MINTS(I) .EQ. 0) THEN
               I = MREALS(I) - 1
            END IF
C 
         ELSE IF (MTYPE(I) .EQ. SKIP) THEN
            I = MREALS(I) - 1
C 
         ELSE IF (MTYPE(I) .EQ. WAIT) THEN
            WRITE(OUNIT, 10) MPRMPT(I)
            CALL PAUSE('to continue',1)
C 
         ELSE IF (MTYPE(I) .EQ. CLRTRM) THEN
              CALL CLEAR
C 
         ELSE IF (MTYPE(I) .EQ. QUIT) THEN
              I=END
C 
         ELSE IF (MTYPE(I) .EQ. PPRINT) THEN
            ILOW = MLOLIM(I)
            IHIGH = MHILIM(I)      
            CALL CLEAR
            CALL PRTAGS(ILOW,IHIGH)
            CALL PAUSE('to see the Main Menu...',1)
            CALL CLEAR
C 
         END IF
C        display info about defaults
         IF (I .EQ. 1) THEN
             WRITE (OUNIT,*) ' '
             WRITE (OUNIT,100)
             WRITE (OUNIT,101)
         END IF
C
      I = I+1            
      IF (I .LE. END) GO TO 1000
C 
10    FORMAT( 1X, A )
11    FORMAT( 1X, A, 1PE11.3 )
12    FORMAT( 1X, A, I6 )
100   FORMAT (' To accept the default value [in brackets] for any item')
101   FORMAT (' just press Return at the prompt')
C
      RETURN
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRTAGS(START,END)
C prints menu prompts and default values for items START to END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'        
      INCLUDE 'IO.ALL'          
C Input variables:
      INTEGER START,END         !limiting indices of printed menu items 
C Local variables:
      INTEGER I                 !menu items index
      INTEGER NUMSKP            !number of lines to skip
      INTEGER INDEX             !subindex for menu items
      INTEGER PLEN              !length of prompt
      INTEGER ICHOIC            !menu/parameter choice
C Functions:
      INTEGER LENTRU            !true length of character string
      INTEGER PARSE             !menu choice
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      I=START
1000  CONTINUE
         IF (MTYPE(I) .EQ. FLOAT) THEN
            WRITE (OUNIT, 11) MTAG(I), MREALS(I)
C 
         ELSE IF (MTYPE(I) .EQ. NUM) THEN
            WRITE (OUNIT, 12) MTAG(I), MINTS(I)
C 
         ELSE IF (MTYPE(I) .EQ. BOOLEN) THEN
            CALL PRYORN(MTAG(I), MINTS(I))
C 
         ELSE IF (MTYPE(I) .EQ. CHSTR) THEN
            WRITE( OUNIT, 13) MTAG(I), MSTRNG(MINTS(I))
C 
         ELSE IF (MTYPE(I) .EQ. TITLE) THEN
            NUMSKP = MLOLIM(I)
            CALL PRBLKS(NUMSKP)
            WRITE (OUNIT, 10) MPRMPT(I)
            NUMSKP = MHILIM(I)
            CALL PRBLKS(NUMSKP)
C 
         ELSE IF (MTYPE(I) .EQ. YESKIP) THEN
            CALL PRYORN(MTAG(I), MINTS(I))
            IF (MINTS(I) .NE. 0 .AND. MREALS(I) .GT. I) THEN
               I = MREALS(I) - 1
            END IF
C 
         ELSE IF (MTYPE(I) .EQ. NOSKIP) THEN
            CALL PRYORN(MTAG(I), MINTS(I))
            IF (MINTS(I) .EQ. 0 .AND. MREALS(I) .GT. I) THEN
               I = MREALS(I) - 1
            END IF
C 
         ELSE IF (MTYPE(I) .EQ. SKIP) THEN
            IF (MREALS(I) .GT. I) I=MREALS(I) - 1  !don't skip backwards
C 
         ELSE IF (MTYPE(I) .EQ. MCHOIC)  THEN
            IF (MREALS(I) .GT. 0) THEN
C            for menu choices that are parameter choices, print out 
C            choice, but first you must find it
             DO 20 INDEX=I-MHILIM(I),I-1
                IF ((I+MREALS(I)-MHILIM(I)-1) .EQ. INDEX) THEN
                  PLEN=LENTRU(MPRMPT(INDEX))
                  WRITE (OUNIT,10) MPRMPT(INDEX)(4:PLEN)
                END IF
20           CONTINUE
            END IF
            IF (MREALS(I) .NE. 0) THEN
               !branch to chosen parameter
               ICHOIC=ABS(INT(MREALS(I)))
               I=PARSE(MTAG(I),ICHOIC)-1
            END IF
C
         END IF
      I = I+1                                           
      IF (I .LE. END) GO TO 1000
C 
10    FORMAT( 1X, A )
11    FORMAT( 1X, A, 1PE11.3 )
12    FORMAT( 1X, A, I6 )  
13    FORMAT( 1X, A, 5X, A)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRYORN(PMPT,YORN)             
C print a 'yes' or 'no'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'             
C Input variables:
      INTEGER  YORN                !1 == 'YES', '0' == 'NO'
      CHARACTER*(*) PMPT           !string to print before y/n
C Functions:
      INTEGER LENTRU               !actual length of prompt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (YORN .EQ. 0) THEN
         WRITE(OUNIT, 10) PMPT(1:LENTRU(PMPT))
      ELSE
         WRITE(OUNIT, 11) PMPT(1:LENTRU(PMPT))
      END IF
10    FORMAT( 1X, A, ': no')
11    FORMAT( 1X, A, ': yes')
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRBLKS(NUMLIN) 
C prints NUMLIN blank lines on terminal
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'               
C Passed variables:
      INTEGER NUMLIN                 !number of blank lines to print
C Local variables:
      INTEGER I                      !dummy index
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 1000 I=1,NUMLIN
         WRITE( OUNIT,*) ' '
1000  CONTINUE
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PAUSE(PHRASE,NSKIP)
C gives user time to read screen by waiting for dummy input;
C allows for printing of PHRASE to screen;
C skips NSKIP lines before printing PHRASE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Passed variables:
      CHARACTER*(*) PHRASE             !phrase to be printed
      INTEGER NSKIP                    !number of lines to skip
C Local variables:
      CHARACTER*1 DUMMY                !dummy variable
      INTEGER ISKIP                    !NSKIP index
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 10 ISKIP=1,NSKIP              !skip lines
         WRITE (OUNIT,5)
10    CONTINUE
5     FORMAT (' ')
      WRITE (OUNIT,15) PHRASE          !write phrase
      READ (IUNIT,20) DUMMY            !wait for dummy input
15    FORMAT (' Press return ',A)
20    FORMAT (A1)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FLOPEN(FNAME,FUNIT)
C opens a new file, unless one by the same name already exists
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:
      CHARACTER*(*) FNAME  !file name
      INTEGER FUNIT        !unit number
C Local variables:
      LOGICAL OPN          !is the file open?
      LOGICAL EXST         !does it exist?
      CHARACTER*40 CHARAC  !function that return character input
      INTEGER LENTRU
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
10    INQUIRE(FILE=FNAME,EXIST=EXST,OPENED=OPN)
C                  
        IF (OPN) RETURN
C
        IF (EXST) THEN
            WRITE (OUNIT,20) FNAME(1:LENTRU(FNAME))
20          FORMAT (' Output file ',A,' already exists')
            FNAME=CHARAC(FNAME,12, 'Enter another filename')
        ELSE
            OPEN(UNIT=FUNIT,FILE=FNAME,STATUS='NEW')
            RETURN
        END IF
      GOTO 10
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FLOPN2(FNAME,FUNIT,SUCESS)
C opens an existing file for input data 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:
      CHARACTER*(*) FNAME  !file name
      INTEGER FUNIT        !unit number
      LOGICAL SUCESS       !did we find an existing file to open?
C Local variables:
      LOGICAL OPN          !is the file open?
      LOGICAL EXST         !does it exist?
      CHARACTER*40 CHARAC  !function that return character input
      INTEGER CHOICE       !choice for continuing
C Functions:
      INTEGER YESNO        !get yes or no input
      INTEGER LENTRU
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
10    INQUIRE(FILE=FNAME,EXIST=EXST,OPENED=OPN)
C                  
        IF ((.NOT. EXST) .OR. (OPN)) THEN
           WRITE (OUNIT,20) FNAME(1:LENTRU(FNAME))
20         FORMAT(' Input file ',A,' does not exist or is already open')
           CHOICE=YESNO(1,' Would you like to try another file name?')
C
           IF (CHOICE .EQ. 0) THEN
               SUCESS=.FALSE.
               RETURN      !leave without opening file for reading
           ELSE        
             FNAME=CHARAC(FNAME,12, 'Enter another filename')
           END IF
        ELSE
            OPEN(UNIT=FUNIT,FILE=FNAME,STATUS='OLD')
            SUCESS=.TRUE.
            RETURN
        END IF
      GOTO 10
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FLCLOS(FNAME,FUNIT)
C checks on file status of file, and closes if open
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:
      CHARACTER*(*) FNAME  !file name
      INTEGER FUNIT        !unit number
C Local variables:
      LOGICAL OPN          !is the file open
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INQUIRE(FILE=FNAME,OPENED=OPN)
      IF (OPN) CLOSE(UNIT=FUNIT)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE FINISH 
C closes files and stops execution
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL FLCLOS(TNAME,TUNIT)
      CALL FLCLOS(GNAME,GUNIT)
      STOP
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FLTDEF(XPRMPT,X)
C prints prompt for floating number
C and displays default X in a format dictated by size of X
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:
      CHARACTER*(*) XPRMPT             !prompt string
      REAL X                           !default value
C Function:
      INTEGER LENTRU                   !true length of string
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     positive numbers (leave no room for a sign)
      IF (X .GT. 0) THEN
       IF ((ABS(X) .LT. 999.49) .AND. (ABS(X) .GE. 99.949)) THEN
          WRITE (OUNIT,5) XPRMPT(1:LENTRU(XPRMPT)),X
       ELSE IF ((ABS(X) .LT. 99.949) .AND. (ABS(X) .GE. 9.9949)) THEN
          WRITE (OUNIT,10) XPRMPT(1:LENTRU(XPRMPT)),X
       ELSE IF ((ABS(X) .LT. 9.9949) .AND. (ABS(X) .GE. .99949)) THEN
          WRITE (OUNIT,15) XPRMPT(1:LENTRU(XPRMPT)),X
       ELSE IF ((ABS(X) .LT. .99949) .AND. (ABS(X) .GE. .099949)) THEN
          WRITE (OUNIT,20) XPRMPT(1:LENTRU(XPRMPT)),X
       ELSE 
          WRITE (OUNIT,25) XPRMPT(1:LENTRU(XPRMPT)),X
       END IF
C 
C     negative numbers (leave room for the sign)
      ELSE
       IF ((ABS(X) .LT. 999.49) .AND. (ABS(X) .GE. 99.949)) THEN
          WRITE (OUNIT,105) XPRMPT(1:LENTRU(XPRMPT)),X
       ELSE IF ((ABS(X) .LT. 99.949) .AND. (ABS(X) .GE. 9.9949)) THEN
          WRITE (OUNIT,110) XPRMPT(1:LENTRU(XPRMPT)),X
       ELSE IF ((ABS(X) .LT. 9.9949) .AND. (ABS(X) .GE. .99949)) THEN
          WRITE (OUNIT,115) XPRMPT(1:LENTRU(XPRMPT)),X
       ELSE IF ((ABS(X) .LT. .99949) .AND. (ABS(X) .GE. .099949)) THEN
          WRITE (OUNIT,120) XPRMPT(1:LENTRU(XPRMPT)),X
       ELSE 
          WRITE (OUNIT,125) XPRMPT(1:LENTRU(XPRMPT)),X
       END IF
                          
      END IF
C              
5     FORMAT (1X,A,1X,'[',F4.0,']')
10    FORMAT (1X,A,1X,'[',F4.1,']')
15    FORMAT (1X,A,1X,'[',F4.2,']')
20    FORMAT (1X,A,1X,'[',F4.3,']')
25    FORMAT (1X,A,1X,'[',1PE8.2,']')
105   FORMAT (1X,A,1X,'[',F5.0,']')
110   FORMAT (1X,A,1X,'[',F5.1,']')
115   FORMAT (1X,A,1X,'[',F5.2,']')
120   FORMAT (1X,A,1X,'[',F5.3,']')
125   FORMAT (1X,A,1X,'[',1PE9.2,']')
C 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTDEF(KPRMPT,K)
C prints prompt for integer input from screen
C and default value in appropriate format
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:
      CHARACTER *(*) KPRMPT       !prompt string
      INTEGER K                   !default values
C Function:
      INTEGER LENTRU              !true length of string
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     positive numbers (leave no room for a sign)
      IF (K .GE. 0 ) THEN
         IF ((IABS(K) .LE. 9999) .AND. (IABS(K) .GE. 1000)) THEN
            WRITE (OUNIT,10) KPRMPT(1:LENTRU(KPRMPT)),K
         ELSE IF ((IABS(K) .LE. 999) .AND. (IABS(K) .GE. 100)) THEN
            WRITE (OUNIT,20) KPRMPT(1:LENTRU(KPRMPT)),K
         ELSE IF ((IABS(K) .LE. 99) .AND. (IABS(K) .GE. 10)) THEN
            WRITE (OUNIT,30) KPRMPT(1:LENTRU(KPRMPT)),K
         ELSE IF ((IABS(K) .LE. 9) .AND. (IABS(K) .GE. 0)) THEN
            WRITE (OUNIT,40) KPRMPT(1:LENTRU(KPRMPT)),K
         ELSE
            WRITE (OUNIT,50) KPRMPT(1:LENTRU(KPRMPT)),K
         END IF
C
C     negative numbers (leave room for the sign)
      ELSE
         IF ((IABS(K) .LE. 9999) .AND. (IABS(K) .GE. 1000)) THEN
            WRITE (OUNIT,110) KPRMPT(1:LENTRU(KPRMPT)),K
         ELSE IF ((IABS(K) .LE. 999) .AND. (IABS(K) .GE. 100)) THEN
            WRITE (OUNIT,120) KPRMPT(1:LENTRU(KPRMPT)),K
         ELSE IF ((IABS(K) .LE. 99) .AND. (IABS(K) .GE. 10)) THEN
            WRITE (OUNIT,130) KPRMPT(1:LENTRU(KPRMPT)),K
         ELSE IF ((IABS(K) .LE. 9) .AND. (IABS(K) .GE. 1)) THEN
            WRITE (OUNIT,140) KPRMPT(1:LENTRU(KPRMPT)),K
         ELSE
            WRITE (OUNIT,150) KPRMPT(1:LENTRU(KPRMPT)),K
         END IF           
      END IF
C  
10    FORMAT (1X,A,1X,'[',I4,']')
20    FORMAT (1X,A,1X,'[',I3,']')
30    FORMAT (1X,A,1X,'[',I2,']')
40    FORMAT (1X,A,1X,'[',I1,']')
50    FORMAT (1X,A,1X,'[',I10,']')
110   FORMAT (1X,A,1X,'[',I5,']')
120   FORMAT (1X,A,1X,'[',I4,']')
130   FORMAT (1X,A,1X,'[',I3,']')
140   FORMAT (1X,A,1X,'[',I2,']')
150   FORMAT (1X,A,1X,'[',I10,']')
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE CONVRT(X,STRING,LEN)
C converts a real number x to a character variable string of length LEN
C for printing; the format is chosen according to the value of X,
C taking roundoff into account
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
       CHARACTER*9 STRING                   !routine output
       REAL X                               !routine input
       INTEGER LEN                          !string length
C Function 
       INTEGER LENTRU                       !gets string length
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     positive numbers (leave no room for a sign)
      IF (X .GT. 0) THEN
       IF ((ABS(X) .LT. 999.4) .AND. (ABS(X) .GE. 99.94)) THEN
          WRITE(STRING,5) X 
       ELSE IF ((ABS(X) .LT. 99.94) .AND. (ABS(X) .GE. 9.994)) THEN
          WRITE(STRING,10)  X
       ELSE IF ((ABS(X) .LT. 9.994) .AND. (ABS(X) .GE. .9994)) THEN
          WRITE (STRING,15) X
       ELSE IF ((ABS(X) .LT. .9994) .AND. (ABS(X) .GE. .09994)) THEN
          WRITE (STRING,20) X
       ELSE 
          WRITE (STRING,25) X
       END IF
C
C     negative numbers (leave room for the sign)
      ELSE
       IF ((ABS(X) .LT. 999.4) .AND. (ABS(X) .GE. 99.94)) THEN
          WRITE(STRING,105) X 
       ELSE IF ((ABS(X) .LT. 99.94) .AND. (ABS(X) .GE. 9.994)) THEN
          WRITE(STRING,110)  X
       ELSE IF ((ABS(X) .LT. 9.994) .AND. (ABS(X) .GE. .9994)) THEN
          WRITE (STRING,115) X
       ELSE IF ((ABS(X) .LT. .9994) .AND. (ABS(X) .GE. .09994)) THEN
          WRITE (STRING,120) X
       ELSE 
          WRITE (STRING,125) X
       END IF
      END IF
C      
      LEN=LENTRU(STRING)
C
5     FORMAT (F4.0)
10    FORMAT (F4.1)
15    FORMAT (F4.2)
20    FORMAT (F4.3)
25    FORMAT (1PE8.2)
105   FORMAT (F5.0)
110   FORMAT (F5.1)
115   FORMAT (F5.2)
120   FORMAT (F5.3)
125   FORMAT (1PE9.2)
C
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE ICNVRT(I,STRING,LEN)
C converts an integer I to a character variable STRING for 
C printing; the format is chosen according to the value of I
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
       CHARACTER*9 STRING                   !routine output
       INTEGER I                            !routine input
       INTEGER LEN                          !length of string
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     positive numbers (leave no room for a sign)
      IF (I .GE. 0) THEN
       IF ((ABS(I) .LE. 9) .AND. (ABS(I) .GE. 0)) THEN
          WRITE(STRING,5) I 
          LEN=1
       ELSE IF ((ABS(I) .LE. 99) .AND. (ABS(I) .GE. 10)) THEN
          WRITE(STRING,10)  I
          LEN=2
       ELSE IF ((ABS(I) .LE. 999) .AND. (ABS(I) .GE. 100)) THEN
          WRITE (STRING,15) I
          LEN=3
       ELSE IF ((ABS(I) .LE. 9999) .AND. (ABS(I) .GE. 1000)) THEN
          WRITE (STRING,20) I
          LEN=4
       ELSE 
          WRITE (STRING,25) REAL(I)
          LEN=8
       END IF
C
C     negative numbers (leave room for the sign)
      ELSE
       IF ((ABS(I) .LE. 9) .AND. (ABS(I) .GE. 1)) THEN
          WRITE(STRING,105) I 
          LEN=2
       ELSE IF ((ABS(I) .LE. 99) .AND. (ABS(I) .GE. 10)) THEN
          WRITE(STRING,110)  I
          LEN=3
       ELSE IF ((ABS(I) .LE. 999) .AND. (ABS(I) .GE. 100)) THEN
          WRITE (STRING,115) I
          LEN=4
       ELSE IF ((ABS(I) .LE. 9999) .AND. (ABS(I) .GE. 1000)) THEN
          WRITE (STRING,120) I
          LEN=5
       ELSE 
          WRITE (STRING,125) REAL(I)
          LEN=9
       END IF
      END IF
C
5     FORMAT (I1)
10    FORMAT (I2)
15    FORMAT (I3)
20    FORMAT (I4)
25    FORMAT (1PE8.2)
105   FORMAT (I2)
110   FORMAT (I3)
115   FORMAT (I4)
120   FORMAT (I5)
125   FORMAT (1PE9.2)
C
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER FUNCTION PARSE(STRING,CHOICE)
C determines branching in menu list
C
C breaks STRING (of the form 'nn nn nn nn nn nn ....') into pieces, and
C returns the integer value represented by the CHOICE group of digits
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Input variables:
      CHARACTER*(*) STRING         !string to look at
      INTEGER CHOICE               !specific number to look at
C Local variables:
      INTEGER IPOS                 !current character position in string
      INTEGER IGROUP               !current group of digits in string
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IPOS=1
      DO 20 IGROUP = 1,CHOICE-1
40      IF ( STRING(IPOS:IPOS) .NE. ' ' ) THEN
            IPOS = IPOS+1
        GOTO 40  
        END IF
        IPOS=IPOS+1
20    CONTINUE
      READ( STRING(IPOS:IPOS+2),10) PARSE
10    FORMAT( I2 )
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER FUNCTION LENTRU(CHARAC)
C finds the true length of a character string by searching
C backward for first nonblank character
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Input variables:
      CHARACTER *(*) CHARAC    !string whose length we are finding
C Local variables:
      INTEGER ISPACE           !ascii value of a blank
      INTEGER I                !index of elements in CHARAC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ISPACE=ICHAR(' ')
      I=LEN(CHARAC)
10    IF (ICHAR(CHARAC(I:I)) .EQ. ISPACE) THEN
         I=I-1
      ELSE
         LENTRU=I
         RETURN
      END IF  
      IF (I .GT. 0) GOTO 10
      LENTRU=0
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION GETFLT(X,XMIN,XMAX,XPRMPT)
C get a floating point number GETFLT; make sure it is between XMIN
C and XMAX and prompt with XPRMPT
C
C If your compiler accepts (FMT=*) to an internal unit, comment out 
C lines 3 and 5, and uncomment lines 2 and 4
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:
      CHARACTER*(*) XPRMPT         !prompt
      REAL X                       !default value
      REAL XMIN,XMAX               !limits on input
C Local variables:
      CHARACTER*40 STRING          !internal unit
C Function
      INTEGER LENTRU               !returns true length of string
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     prompt for float, display default value
10    CALL FLTDEF(XPRMPT,X)        
      READ (IUNIT,35,ERR=10) STRING
C  
C     accept default value X if STRING is empty
      IF (LENTRU(STRING) .EQ. 0) THEN   
          GETFLT=X
      ELSE
C2          READ (UNIT=STRING,FMT=*,ERR=10) GETFLT
3         READ (UNIT=STRING,FMT=1,ERR=10) GETFLT
1         FORMAT (E9.2)
      END IF
C
C     make sure GETFLT is between XMIN and XMAX
40    IF ((GETFLT .LT. XMIN) .OR. (GETFLT .GT. XMAX)) THEN
50       WRITE (OUNIT,60) XMIN,XMAX
         READ (IUNIT,35,ERR=50) STRING
         IF (LENTRU(STRING) .EQ. 0) THEN
             GETFLT=X
         ELSE                        
C4             READ (UNIT=STRING,FMT=*,ERR=50) GETFLT
5             READ (UNIT=STRING,FMT=1,ERR=50) GETFLT
         END IF
      GOTO 40
      END IF
C
35    FORMAT (A40)
60    FORMAT (' Try again: input outside of range = [',1PE11.3,
     +        1PE11.3,']')
100   FORMAT (1PE9.2)

      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER FUNCTION GETINT(K,KMIN,KMAX,KPRMPT)
C get an integer value GETINT;
C check that it lies between KMIN and KMAX and prompt with KPRMPT
C
C This function allows input of integers in a natural way (i.e., without
C preceding blanks or decimal points) even though we cannot use list
C directed READ (i.e., FMT=*) from internal units
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:
      CHARACTER*(*) KPRMPT         !string prompt
      INTEGER K                    !default value
      INTEGER KMIN,KMAX            !upper and lower limits
C Local variables:
      CHARACTER*40 STRING          !internal unit
      REAL TEMP                    !temp var to allow for easier input
C Functions:
      INTEGER LENTRU               !returns true length of string
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     prompt for input; display default
10    CALL INTDEF(KPRMPT,K)        
      READ (IUNIT,35,ERR=10) STRING
C
C     accept default value K if STRING is empty
      IF (LENTRU(STRING) .EQ. 0) THEN
          GETINT=K
      ELSE
C         change the integer into a real number
          STRING=STRING(1:LENTRU(STRING))//'.'
C         read the real number from string
          READ (UNIT=STRING,FMT=1,ERR=10) TEMP
1         FORMAT(F7.0)
C         change it to an integer
          GETINT=INT(TEMP)
      END IF
C
C     check that GETINT lies between KMIN and KMAX
40    IF ((GETINT .LT. KMIN) .OR. (GETINT .GT. KMAX)) THEN          
50       WRITE (OUNIT,60) KMIN,KMAX
         READ (IUNIT,35,ERR=50) STRING
         IF (LENTRU(STRING) .EQ. 0) THEN
             GETINT=K
         ELSE                        
             STRING=STRING(1:LENTRU(STRING))//'.'
             READ (UNIT=STRING,FMT=1,ERR=50) TEMP
             GETINT=INT(TEMP)
         END IF
      GOTO 40
      END IF            
C
35    FORMAT (A40)
60    FORMAT (' Try again: input is outside of range = [',I6,I6,']')
100   FORMAT (I10)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CHARACTER*40 FUNCTION CHARAC(C,CLNGTH,CPRMPT)
C gets character string CHARAC no longer than CLNGTH from the screen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:
      CHARACTER*(*) C                 !default value
      CHARACTER*(*) CPRMPT            !prompt
      INTEGER CLNGTH                  !max length
C Local variables:
      CHARACTER*40 STRING             !internal unit
C Functions:
      INTEGER LENTRU                  !returns true length of string
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     data can't be longer than 40 characters due to fixed format 
      IF (CLNGTH .GT. 40) CLNGTH=40
C
C     prompt for string; display default value C 
10    WRITE (OUNIT,20) CPRMPT(1:LENTRU(CPRMPT)),C(1:LENTRU(C))
      READ (IUNIT,35,ERR=10) STRING
C
C     accept default value C if STRING is empty
      IF (LENTRU(STRING) .EQ. 0) THEN                         
          CHARAC=C
      ELSE
          READ (STRING,35,ERR=10) CHARAC
      END IF                                       
C                 
C     find the true length of the input; verify that it is not too long
40    IF (LENTRU(CHARAC) .GT. CLNGTH) THEN
50       WRITE (OUNIT,60) CLNGTH
         READ (IUNIT,35,ERR=50) STRING
         IF (LENTRU(STRING) .EQ. 0) THEN
             CHARAC=C
         ELSE                        
             READ (STRING,35,ERR=50) CHARAC
         END IF
      GOTO 40
      END IF
C
20    FORMAT (1X,A,1X,'[',A,']')
35    FORMAT (A40)
60    FORMAT (' Try again: string is too long, maximum length = ',I2)
C
      RETURN
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER FUNCTION YESNO(BINARY,PROMPT)
C obtains YESNO from the screen; value is 0 for no, 1 for yes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input parameters:
      CHARACTER*(*) PROMPT        !prompt
      INTEGER BINARY              !default value
C Local variables:  
      CHARACTER*3 STRING          !internal unit
C Functions:
      INTEGER LENTRU              !returns true length of string
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
1000  CONTINUE
C     write prompt and display default values
      IF (BINARY .EQ. 1) WRITE(OUNIT,10) PROMPT(1:LENTRU(PROMPT))
      IF (BINARY .EQ. 0) WRITE(OUNIT,11) PROMPT(1:LENTRU(PROMPT))
C 
      READ (IUNIT,20,ERR=1000) STRING
C 
C     accept default value; check that input is 'y' or 'n'
      IF (LENTRU(STRING) .EQ. 0) THEN
         YESNO = BINARY
      ELSE IF (STRING(1:1) .EQ. 'y' .OR. STRING(1:1) .EQ. 'Y') THEN
         YESNO = 1
      ELSE IF (STRING(1:1) .EQ. 'n' .OR. STRING(1:1) .EQ. 'N') THEN
         YESNO = 0
      ELSE
           WRITE (OUNIT,200)
           GOTO 1000
      END IF
C 
10    FORMAT(1X,A,1X,'[yes]')
11    FORMAT(1X,A,1X,'[no]')
20    FORMAT(A)
200   FORMAT (' Try again, answer must be yes or no')
C
      RETURN
      END   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      LOGICAL FUNCTION LOGCVT(IJK)
C converts 1 to true and anything else to false
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER IJK                    !input
      IF (IJK .EQ. 1) THEN
         LOGCVT=.TRUE.
      ELSE
         LOGCVT=.FALSE.
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION RANNOS(DSEED)
C returns a uniformly distributed random number between 0 and 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION         DSEED
      DOUBLE PRECISION         D2P31M,D2P31
      DATA             D2P31M/2147483647.D0/
      DATA             D2P31 /2147483711.D0/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DSEED = MOD(16807.D0*DSEED,D2P31M)
      RANNOS = DSEED / D2P31
      RETURN
      END
                                                               
