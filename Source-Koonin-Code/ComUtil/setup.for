C file SETUP.FOR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE SETUP 
C allows users to supply i/o parameters for their computing environment
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'IO.ALL'   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      fortran unit numbers for i/o
C      unit for text output to a file
       TUNIT=10                          
C      unit for graphics output to file
       GUNIT=20
C      unit for input from keyboard
       IUNIT=5
C      unit for output to screen
       OUNIT=6
C      unit for input of data
       DUNIT=11
C
C      how many lines and columns of text fit on your screen?
       TRMLIN=24
       TRMWID=80
C
C      default output parameters 
C      There are five forms of output provided, here you are choosing
C      which forms of output you will want MOST of the time (any 
C      combination is possible), you always have the option to change
C      your mind at run time.
C      0=no   1=yes
C      do you want text sent to the screen?
       TXTTRM=1
C      do you want text sent to a file?
       TXTFIL=1
C      do you want graphics sent to the screen?
       GRFTRM=0
C      do you want graphics sent to a hardcopy device?
       GRFHRD=0
C      do you want graphics data sent to a file?
       GRFFIL=0
C
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE CLEAR
C clears text screen by sending an escape sequence;
C check your terminal manual for the correct sequence 
C THIS IS NOT AN ESSENTIAL ROUTINE - YOU CAN LEAVE IT BLANK
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
        INCLUDE 'IO.ALL'   
C Local variables:
        CHARACTER*1 ESC1(4),ESC2(6)     !escape characters 
        INTEGER I,I1,I2                 !index of escape sequence arrays
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C VT200 terminal; text mode
        ESC1(1)=CHAR(27)                !<ESC>[2J
        ESC1(2)=CHAR(91)
        ESC1(3)=CHAR(50)
        ESC1(4)=CHAR(74)
        I1=4

        ESC2(1)=CHAR(27)                !<ESC>11;1f
        ESC2(2)=CHAR(91)
        ESC2(3)=CHAR(49)
        ESC2(4)=CHAR(59)
        ESC2(5)=CHAR(49)
        ESC2(6)=CHAR(102)
        I2=6
C
C TEK4010
C        ESC1(1)=CHAR(27)               !<ESC><FF>
C        ESC1(2)=CHAR(12)
C        I1=2
C        I2=0
C
C PST (Prime)
C        ESC1(1)=CHAR(27)                !<ESC>?
C        ESC1(2)=CHAR(63)
C        I1=2
C        I2=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        WRITE(OUNIT,10) (ESC1(I),I=1,I1)         
C        WRITE(OUNIT,10) (ESC2(I),I=1,I2)
10      FORMAT (1X,6A1)  
        RETURN
        END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GMODE
C switches terminal from text to graphics mode
C by writing hardware dependent escape sequences to the terminal
C This routine contains the escape sequence for a Graphon terminal
C to switch between vt200 and tek4014 modes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'   
C Local variables:
      CHARACTER*1 ESC(2)                    
      ESC(1)=CHAR(27)                       !ascii codes for <ESC> 1
      ESC(2)=CHAR(49)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      WRITE(OUNIT,10) ESC(1),ESC(2)         
10    FORMAT (1X,2A1)  
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TMODE
C switches terminal from graphics to text mode
C by writing hardware dependent escape sequences to the terminal
C This routine contains the escape sequence for a Graphon terminal
C to switch between tek4014 and vt200 modes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'   
C Local variables:       
        CHARACTER*1 ESC(2)                
        ESC(1)=CHAR(27)                   !ascii codes for <ESC> 2  
        ESC(2)=CHAR(50)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      WRITE(OUNIT,10) ESC(1),ESC(2)
10    FORMAT (1X,2A1)  
      RETURN
      END
