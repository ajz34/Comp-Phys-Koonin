C STRIP.FOR
C Reads in a FORTRAN code and deletes all comments 
C starting with the non-fortran-77 standard comment delimiter: "!";
C also deletes all blank lines
C
       CHARACTER*80 LINE
       CHARACTER*80 INFILE,OUTFIL                             
       CHARACTER*80 ANSWER,FIRST
       LOGICAL EXST,OPN,FOUND
       INTEGER EXCLM,SPACE
       INTEGER TERM,OUNIT,IUNIT
       INTEGER I,ASCII,NONBLK
C
       EXCLM=ICHAR('!')
       SPACE=ICHAR(' ')
C               
C units for keyboard, INFILE, and OUTFIL (change these if necessary)
       TERM=5
       IUNIT=10
       OUNIT=20
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
10     CONTINUE       
C       
C get name of input file and open file
C
       PRINT *, 'ENTER NAME OF FORTRAN PROGRAM'
       READ (TERM,60)  INFILE
C
20     INQUIRE(FILE=INFILE,EXIST=EXST,OPENED=OPN)
          IF (EXST .EQV. .FALSE.) THEN
            PRINT *, 'FILE DOES NOT EXIST'
            PRINT *, 'ENTER ANOTHER NAME'
            READ (TERM,60)  INFILE
          ELSE IF (OPN .EQV. .FALSE.) THEN
            OPEN(UNIT=IUNIT,FILE=INFILE,STATUS='OLD')
            OPN=.TRUE.
          END IF
       IF (OPN .EQV. .FALSE.) GOTO 20
C
C get name of output file and open 
C
       PRINT *, 'ENTER NAME OF NEW FILE'
       READ (TERM,60)    OUTFIL
C
C this line is for a VAX only
C       OPEN(UNIT=OUNIT,FILE=OUTFIL,STATUS='NEW',CARRIAGECONTROL='LIST')
C this line is for any other machine
       OPEN(UNIT=OUNIT,FILE=OUTFIL,STATUS='NEW')
C
C read each line of input file and search for "!"
35     READ(IUNIT,60,END=50) LINE   
         FOUND=.FALSE.
         I=0
         NONBLK=0
45       IF ((FOUND .EQV. .FALSE.) .AND. (I .LT. 80)) THEN
            I=I+1
            ASCII=ICHAR(LINE(I:I))
            IF (ASCII .EQ. EXCLM) THEN
               FOUND=.TRUE.
            ELSE IF (ASCII .NE. SPACE) THEN
               NONBLK=I
            END IF
         GOTO 45
         END IF 
C
C        print up to last nonblank character, exclude "!"
         IF (NONBLK .GT. 0) WRITE (OUNIT,60) LINE(1:NONBLK)
C
       GOTO 35
50     CONTINUE
C    
       CLOSE(UNIT=OUNIT)
       CLOSE(UNIT=IUNIT)
C
C allow for another file
       PRINT *, 'DO YOU WISH TO STANDARDIZE ANOTHER FILE? [Y]'
       READ (TERM,60)  ANSWER
C
       FIRST=ANSWER(1:1)
       IF ((FIRST .EQ. 'N') .OR. (FIRST .EQ. 'n')) STOP       
C
       GOTO 10
C
60     FORMAT (A)
       END
