CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM EXMPL8
C   Example 8: Monte Carlo simulation of the two-dimensional Ising Model
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company, Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop/ execute once for each set of param
        CALL PARAM        !get input from screen
        CALL ARCHON       !calculate the thermodynamic quantities
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON  
C calculates the thermodynamic quantities:  energy, magnetization
C susceptibility, and specific heat at constant field and intrctn strngth
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E8'
      INCLUDE 'IO.ALL'
C Local variables:
      INTEGER SPIN(MAXX,MAXY)    !spin configuration
C     all of the thermodynamic quant have 2 indices 
C     (not all array elements are used, e.g. CHI(sweep,value))
C     first index is the level: sweep, group, or total
C     second index is the value: quantity, quant**2, or sigma**2
      REAL MAG(3,3)              !magnetization
      REAL ENERGY(3,3)           !energy
      REAL CB(3,3)               !specific heat
      REAL CHI(3,3)              !susceptibility
      INTEGER ITHERM             !thermalization index
      INTEGER ITER               !iteration index
      REAL ACCPT                 !acceptance ratio
      INTEGER IX,IY              !horiz and vert indices
      INTEGER NLINES             !number of lines printed to terminal
      INTEGER MORE,IGRP          !how many more groups, group index
      INTEGER ISWEEP             !sweep index
      INTEGER SWEEP,GROUP,TOTAL  !which level of calculation
      INTEGER VALUE,SQUARE,SIGSQ !which quantity
C Functions:
      REAL GETFLT                !get floating point number from screen
      INTEGER YESNO              !get yes/no answer from screen
      LOGICAL LOGCVT             !change from 1,0 to true and false
      INTEGER GETINT             !get integer data from screen
      REAL RANNOS                !generates a random number
      DATA SWEEP,GROUP,TOTAL/1,2,3/
      DATA VALUE,SQUARE,SIGSQ/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     output summary of parameters
      IF (TTERM) CALL PRMOUT(OUNIT,NLINES)
      IF (TFILE) CALL PRMOUT(TUNIT,NLINES)
      IF (GFILE) CALL PRMOUT(GUNIT,NLINES)
C                  
      DO 5 IX=1,NX              !random initial spin configuration
         DO 6 IY=1,NY
              IF (RANNOS(DSEED) .GT. .5) THEN
                 SPIN(IX,IY)=1
              ELSE
                 SPIN (IX,IY)=-1
              END IF
6        CONTINUE
5     CONTINUE
C
      DO 10 ITHERM=1,NTHERM     !thermalize init config
         CALL METROP(SPIN,ACCPT)
         IF (TTERM) WRITE (OUNIT,7) ITHERM,ACCPT
         IF (TFILE) WRITE (TUNIT,7) ITHERM,ACCPT
7        FORMAT (5X,' thermalization sweep ',I4,
     +           ', acceptance ratio =',F5.3) 
10    CONTINUE            
      IF (TTERM) CALL PAUSE('to begin summing...',1)
C
      CALL ZERO(TOTAL,MAG,ENERGY,CHI,CB)  !zero total averages
      MORE=NGROUP
C
15    CONTINUE                                !allow for more groups
        DO 20 IGRP=NGROUP-MORE+1,NGROUP       !loop over groups
           CALL ZERO(GROUP,MAG,ENERGY,CHI,CB) !zero group averages
C
           IF ((TTERM) .AND. (.NOT. GTERM) .AND. (.NOT. TERSE)) 
     +              CALL TITLES(OUNIT,NLINES)
           IF ((TFILE) .AND. (.NOT. TERSE)) CALL TITLES(TUNIT,NLINES)
C
           DO 30 ITER=1,NFREQ*NSIZE
              CALL METROP(SPIN,ACCPT)     !make a sweep of the lattice
C                               
              IF (MOD(ITER,NFREQ) .EQ. 0) THEN  !include in averages 
                  ISWEEP=ITER/NFREQ               !which sweep is it
                  CALL ZERO(SWEEP,MAG,ENERGY,CHI,CB) !zero sweep totals
                  CALL SUM(SPIN,MAG,ENERGY)!sweep totals, add to group
C
                  IF (GTERM) THEN          !display data
                     CALL DISPLY(OUNIT,SPIN,MAG,
     +                       ENERGY,ACCPT,IGRP,ISWEEP)
                  ELSE IF ((TTERM) .AND. (.NOT. TERSE)) THEN            
                     CALL  
     +               SWPOUT(OUNIT,MAG,ENERGY,ACCPT,IGRP,ISWEEP,NLINES)
                  END IF
                  IF ((TFILE) .AND. (.NOT. TERSE)) CALL
     +               SWPOUT(TUNIT,MAG,ENERGY,ACCPT,IGRP,ISWEEP,NLINES)
                  IF ((GFILE) .OR. (GHRDCP)) CALL          
     +               DISPLY(GUNIT,SPIN,MAG,ENERGY,ACCPT,IGRP,ISWEEP)
C
              END IF
30         CONTINUE  
           CALL AVERAG(MAG,ENERGY,CHI,CB,IGRP) !calc total averages
20      CONTINUE
C
        MORE=GETINT(10,0,1000,'How many more groups?')
      IF (MORE .GT. 0) THEN
           NGROUP=NGROUP+MORE
           NLINES=0                                                     
           IF ((TTERM) .AND. (.NOT. TERSE))CALL CLEAR
           GOTO 15
      END IF
C     
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE METROP(SPIN,ACCPT)
C make one sweep of the lattice using the Metropolis algorithm 
C to generate a new configuration
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E8'
C Input/Output variables:
      INTEGER SPIN(MAXX,MAXY)     !spin configuration (I/O)
      REAL ACCPT                  !acceptance ratio (output)
C Local variables:
      INTEGER IX,IY               !horiz and vert indices
      INTEGER IXM1,IXP1,IYM1,IYP1 !indices of nearest neighbors
      INTEGER NNSUM               !sum of nearest neighbors
C Function:
      REAL RANNOS                 !generates a random number
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ACCPT=0.                          !zero acceptance ratio
      DO 10 IX=1,NX              
        IXP1=IX+1                       !nearest neighbors
        IF (IX .EQ. NX) IXP1=1          !with periodic b.c.
        IXM1=IX-1
        IF (IX .EQ. 1) IXM1=NX
C     
        DO 20 IY=1,NY
          IYP1=IY+1                     !nearest neighbors
          IF (IY .EQ. NY) IYP1=1        !with periodic b.c.
          IYM1=IY-1
          IF (IY .EQ. 1) IYM1=NY
C                                       !term to weight new config
          NNSUM=SPIN(IX,IYP1)+SPIN(IX,IYM1)+SPIN(IXP1,IY)+SPIN(IXM1,IY)
C
          IF (RANNOS(DSEED) .LT. RATIO(NNSUM,SPIN(IX,IY))) THEN
             SPIN(IX,IY)=-SPIN(IX,IY)   !flip the spin
             ACCPT=ACCPT+1              !update accept count
          END IF
C
20      CONTINUE
10    CONTINUE
      ACCPT=ACCPT/NSPIN                 !make it a ratio
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SUM(SPIN,MAG,ENERGY)
C calculate magnetization and energy for this sweep
C add these values to the group averages
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E8'
C Input/output variables:
      INTEGER SPIN(MAXX,MAXY)    !spin configuration (input)
C     all of the thermodynamic quant have 2 indices
C     (not all array elements are used, e.g. CHI(sweep,value))
C     first index is the level: sweep, group, or total
C     second index is the quantity: value, square, or sigma**2
      REAL MAG(3,3)              !magnetization (I/O)
      REAL ENERGY(3,3)           !energy (I/O)
C Local variables:
      INTEGER PAIRS              !interaction sum
      INTEGER SWEEP,GROUP,TOTAL  !which level of calculation
      INTEGER VALUE,SQUARE,SIGSQ !which quantity
      INTEGER IX,IY              !horiz and vert indices
      INTEGER IXM1,IYM1          !neighbor indices
      DATA SWEEP,GROUP,TOTAL/1,2,3/
      DATA VALUE,SQUARE,SIGSQ/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PAIRS=0                      !zero pair sum
      DO 10 IY=1,NY
         IYM1=IY-1                 !neighbor just below
         IF (IY .EQ. 1) IYM1=NY    !periodic b.c.
         DO 20 IX=1,NX
           IXM1=IX-1               !neighbor to the left
           IF (IX .EQ. 1) IXM1=NX  !periodic b.c.
C
C          this method of summing pairs does not count twice
           PAIRS=PAIRS+SPIN(IX,IY)*(SPIN(IX,IYM1)+SPIN(IXM1,IY))
C          magnetization is the sum of the spins  (Eq. 8.21a)
           MAG(SWEEP,VALUE)=MAG(SWEEP,VALUE)+SPIN(IX,IY)
C
20       CONTINUE
10    CONTINUE
C
      MAG(SWEEP,SQUARE)=MAG(SWEEP,VALUE)**2
      ENERGY(SWEEP,VALUE)=-J*PAIRS-B*MAG(SWEEP,VALUE)   !Eq 8.18
      ENERGY(SWEEP,SQUARE)=ENERGY(SWEEP,VALUE)**2
C
C     add sweep contributions to group sums
      MAG(GROUP,VALUE)=MAG(GROUP,VALUE)+MAG(SWEEP,VALUE)
      MAG(GROUP,SQUARE)=MAG(GROUP,SQUARE)+MAG(SWEEP,SQUARE)
      ENERGY(GROUP,VALUE)=ENERGY(GROUP,VALUE)+ENERGY(SWEEP,VALUE)
      ENERGY(GROUP,SQUARE)=ENERGY(GROUP,SQUARE)+ENERGY(SWEEP,SQUARE)
C
      RETURN
      END      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE AVERAG(MAG,ENERGY,CHI,CB,IGROUP) 
C find group averages from group sums and add these to total averages;
C calculate uncertainties and display results
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E8'
      INCLUDE 'IO.ALL'
C Input/Output variables:
C     all of the thermodynamic quant have 2 indices
C     (not all array elements are used, e.g. CHI(sweep,value))
C     first index is the level: sweep, group, or total
C     second index is the value: quantity, quant**2, or sigma**2
      REAL MAG(3,3)              !magnetization
      REAL ENERGY(3,3)           !energy
      REAL CB(3,3)               !specific heat
      REAL CHI(3,3)              !susceptibility
      INTEGER IGROUP             !group index (input)
C Local variables:
      REAL M,MSIG1,MSIG2         !magnetization and uncertainties
      REAL E,ESIG1,ESIG2         !energy and uncertainties
      REAL SUS,SUSSIG            !susceptibility and uncertainty
      REAL C,CSIG                !specific heat and uncertainty
      INTEGER IQUANT             !index the quantity
      INTEGER SWEEP,GROUP,TOTAL  !which level of calculation
      INTEGER VALUE,SQUARE,SIGSQ !which quantity
      DATA SWEEP,GROUP,TOTAL/1,2,3/
      DATA VALUE,SQUARE,SIGSQ/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     calculate group averages and uncertainties from group sums
      DO 10 IQUANT=VALUE,SQUARE
         MAG(GROUP,IQUANT)=MAG(GROUP,IQUANT)/NSIZE
         ENERGY(GROUP,IQUANT)=ENERGY(GROUP,IQUANT)/NSIZE
10    CONTINUE
      CHI(GROUP,VALUE)=MAG(GROUP,SQUARE)-MAG(GROUP,VALUE)**2
      MAG(GROUP,SIGSQ)=CHI(GROUP,VALUE)/NSIZE
      IF (MAG(GROUP,SIGSQ) .LT. 0.) MAG(GROUP,SIGSQ)=0.
      CB(GROUP,VALUE)=ENERGY(GROUP,SQUARE)-ENERGY(GROUP,VALUE)**2
      ENERGY(GROUP,SIGSQ)=CB(GROUP,VALUE)/NSIZE
      IF (ENERGY(GROUP,SIGSQ) .LT. 0.) ENERGY(GROUP,SIGSQ)=0.
      CHI(GROUP,SQUARE)=CHI(GROUP,VALUE)**2
      CB(GROUP,SQUARE)=CB(GROUP,VALUE)**2
C
C     add group averages to total sums
      DO 20 IQUANT=VALUE,SIGSQ
         MAG(TOTAL,IQUANT)=MAG(TOTAL,IQUANT)+MAG(GROUP,IQUANT)
         ENERGY(TOTAL,IQUANT)=ENERGY(TOTAL,IQUANT)+ENERGY(GROUP,IQUANT)
         CHI(TOTAL,IQUANT)=CHI(TOTAL,IQUANT)+CHI(GROUP,IQUANT)
         CB(TOTAL,IQUANT)=CB(TOTAL,IQUANT)+CB(GROUP,IQUANT)
20    CONTINUE
C
C     find total averages using total sums accumulated so far
      M=MAG(TOTAL,VALUE)/IGROUP
      MSIG1=(MAG(TOTAL,SQUARE)/IGROUP-M**2)/IGROUP/NSIZE
      IF (MSIG1 .LT. 0) MSIG1=0.
      MSIG1=SQRT(MSIG1)
      MSIG2=SQRT(MAG(TOTAL,SIGSQ))/IGROUP
C
      E=ENERGY(TOTAL,VALUE)/IGROUP
      ESIG1=(ENERGY(TOTAL,SQUARE)/IGROUP-E**2)/IGROUP/NSIZE
      IF (ESIG1 .LT. 0) ESIG1=0.
      ESIG1=SQRT(ESIG1)
      ESIG2=SQRT(ENERGY(TOTAL,SIGSQ))/IGROUP
C
      SUS=CHI(TOTAL,VALUE)/IGROUP
      SUSSIG=(CHI(TOTAL,SQUARE)/IGROUP-SUS**2)/IGROUP
      IF (SUSSIG .LT. 0.) SUSSIG=0.
      SUSSIG=SQRT(SUSSIG)
C
      C=CB(TOTAL,VALUE)/IGROUP
      CSIG=(CB(TOTAL,SQUARE)/IGROUP-C**2)/IGROUP
      IF (CSIG .LT. 0.) CSIG=0.
      CSIG=SQRT(CSIG)
C
C     write out summary
      IF (TTERM) CALL TXTOUT(MAG,ENERGY,CB,CHI,E,ESIG1,ESIG2,
     +    M,MSIG1,MSIG2,SUS,SUSSIG,C,CSIG,IGROUP,OUNIT)
      IF (TFILE) CALL TXTOUT(MAG,ENERGY,CB,CHI,E,ESIG1,ESIG2,
     +    M,MSIG1,MSIG2,SUS,SUSSIG,C,CSIG,IGROUP,TUNIT)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ZERO(ILEVEL,MAG,ENERGY,CHI,CB)
C zero sums for ILEVEL thermodynamic values
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Input/Output variables:
      INTEGER ILEVEL             !which level to zero  (input)
C     all of the thermodynamic quant have 2 indices
C     (not all array elements are used, e.g. CHI(sweep,value))
C     first index is the level: sweep, group, or total
C     second index is the value: quantity, quant**2, or sigma**2
      REAL MAG(3,3)              !magnetization (output)
      REAL ENERGY(3,3)           !energy (output)
      REAL CB(3,3)               !specific heat (output)
      REAL CHI(3,3)              !susceptibility (output)
C Local variable:                                                              
      INTEGER IQUANT             !which quantity
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 10 IQUANT=1,3 
           MAG(ILEVEL,IQUANT)=0.
           ENERGY(ILEVEL,IQUANT)=0.
           CHI(ILEVEL,IQUANT)=0.
           CB(ILEVEL,IQUANT)=0.
10    CONTINUE
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
      INCLUDE 'PARAM.E8'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)         
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get environment parameters
      CALL SETUP                
C 
C     display header screen     
      DESCRP(1)= 'EXAMPLE 8'
      DESCRP(2)= 'Monte Carlo simulation of the 2-D Ising Model'
      DESCRP(3)= 'using the Metropolis algorithm'
      NHEAD=3
C 
C     text output description
      DESCRP(4)= 'acceptance rate, energy, magnetization, specific heat'
      DESCRP(5)= 'and susceptibility (all are values per spin)'
      NTEXT=2
C 
C     graphics output description
      DESCRP(6)= 'spin configuration (blank=-1; X=+1)'
      NGRAPH=1
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C 
C     setup menu arrays, beginning with constant part
      CALL MENU                      
C                
      MTYPE(13)=FLOAT                                        
      MPRMPT(13)='Enter value for magnetic field (units of kT)'
      MTAG(13)='Magnetic field (units of kT)'
      MLOLIM(13)=-20.
      MHILIM(13)=20.
      MREALS(13)=0.
C
      MTYPE(14)=FLOAT
      MPRMPT(14)='Enter value for interaction strength (units of kT)'
      MTAG(14)='interaction strength (units of kT)'
      MLOLIM(14)=-20.
      MHILIM(14)=20.
      MREALS(14)=.3
C
      MTYPE(15)=SKIP
      MREALS(15)=35.
C 
      MTYPE(38)=NUM
      MPRMPT(38)= 'Enter number of X lattice points'
      MTAG(38)= 'Number of X lattice points'
      MLOLIM(38)=2.
      MHILIM(38)=MAXX
      MINTS(38)=20
C 
      MTYPE(39)=NUM
      MPRMPT(39)= 'Enter number of Y lattice points'
      MTAG(39)= 'Number of Y lattice points'
      MLOLIM(39)=2.
      MHILIM(39)=MAXY
      MINTS(39)=20
C 
      MTYPE(40)=NUM
      MPRMPT(40)= 'Integer random number seed for init fluctuations'
      MTAG(40)= 'Random number seed'
      MLOLIM(40)=1000.
      MHILIM(40)=99999.
      MINTS(40)=54767
C 
      MTYPE(41)=NUM
      MPRMPT(41)= 'Number of thermalization sweeps'
      MTAG(41)= 'Thermalization sweeps'
      MLOLIM(41)=0
      MHILIM(41)=1000
      MINTS(41)=20
C 
      MTYPE(42)=NUM
      MPRMPT(42)=  'Enter sampling frequency (to avoid correlations)'
      MTAG(42)= 'Sampling frequency'
      MLOLIM(42)=1
      MHILIM(42)=100
      MINTS(42)=5
C 
      MTYPE(43)=NUM
      MPRMPT(43)= 'Number of samples in a group'
      MTAG(43)= 'Group size'
      MLOLIM(43)=1
      MHILIM(43)=1000
      MINTS(43)=10
C                   
      MTYPE(44)=NUM
      MPRMPT(44)= 'Enter number of groups'
      MTAG(44)= 'Number of groups'
      MLOLIM(44)=1
      MHILIM(44)=1000
      MINTS(44)=10
C 
      MTYPE(45)=SKIP
      MREALS(45)=60.
C 
      MSTRNG(MINTS(75))= 'exmpl8.txt'
C                      
      MTYPE(76)=BOOLEN
      MPRMPT(76)='Do you want the short version of the output?'
      MTAG(76)='Short version of output'
      MINTS(76)=0
C                      
      MTYPE(77)=SKIP
      MREALS(77)=80.
C                                                 
      MSTRNG(MINTS(86))= 'exmpl8.grf'
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
C performs checks on parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:                           
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'                                            
      INCLUDE 'PARAM.E8'
C Local variables:
      INTEGER IF              !possible values for sum of neighb. spins  
C map between menu indices and parameters
      INTEGER IB,IJ,INX,INY,IDSEED,INTHRM,INFREQ,INSIZE,INGRP,ITERSE
      PARAMETER (IB     = 13)
      PARAMETER (IJ     = 14)
      PARAMETER (INX    = 38)
      PARAMETER (INY    = 39)
      PARAMETER (IDSEED = 40)
      PARAMETER (INTHRM = 41)
      PARAMETER (INFREQ = 42)
      PARAMETER (INSIZE = 43)
      PARAMETER (INGRP  = 44)
      PARAMETER (ITERSE = 76)
C Functions:
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
C     set new parameter values   
C     physical and numerical 
      B=MREALS(IB)
      J=MREALS(IJ)
      NX=MINTS(INX)
      NY=MINTS(INY)
      DSEED=DBLE(MINTS(IDSEED))
      NTHERM=MINTS(INTHRM)
      NFREQ=MINTS(INFREQ)
      NSIZE=MINTS(INSIZE)
      NGROUP=MINTS(INGRP)
C
C     text output
      TTERM=LOGCVT(MINTS(ITTERM))
      TFILE=LOGCVT(MINTS(ITFILE))
      TNAME=MSTRNG(MINTS(ITNAME))
      TERSE=LOGCVT(MINTS(ITERSE))
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
      CALL CLEAR
C
C     calculate derivative parameters
      NSPIN=NX*NY
      DO 10 IF=-4,4,2    !ratio of prob.; not all matrix elem are used
            RATIO(IF,-1)=EXP(2*(J*IF+B))
            RATIO(IF,1) =1./RATIO(IF,-1)
10    CONTINUE
C
C     calculate parameters for best looking text display
      IF (2*NX .LE. TRMWID) THEN
           XSKIP=.TRUE.              !skip spaces in x
           XCNTR=(TRMWID-2*NX)/2     !how to center display
      ELSE
           XSKIP=.FALSE.
           XCNTR=(TRMWID-NX)/2
      END IF
      IF (XCNTR .LT. 1) XCNTR=1

      IF (2*NY .LE. TRMLIN-5) THEN
          YSKIP=.TRUE.               !skip lines in y
          YCNTR=(TRMLIN-2*NY)/2-3    !how to center display
      ELSE
          YSKIP=.FALSE.
          YCNTR=(TRMLIN-NY)/2-3
      END IF
      IF (YCNTR .LT. 0) YCNTR=0
C           
      RETURN               
      END    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRMOUT(MUNIT,NLINES)
C write out parameter summary to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E8'
C Input variables:
      INTEGER MUNIT                   !fortran unit number
      INTEGER NLINES                  !number of lines sent to terminal
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) CALL CLEAR
C
      WRITE (MUNIT,5)
      WRITE (MUNIT,6)
      WRITE (MUNIT,7) B
      WRITE (MUNIT,8) J
      WRITE (MUNIT,10) NX,NY
      WRITE (MUNIT,15) NTHERM
      WRITE (MUNIT,20) NFREQ,NSIZE
      WRITE (MUNIT,*) ' '
      NLINES=8
C      
5     FORMAT (' Output from example 8:')
6     FORMAT (' Monte Carlo Simulation of the'
     +        ' 2-D Ising Model using the Metropolis Algorithm')
7     FORMAT (' Magnetic field (units of kT) =',1PE12.5)
8     FORMAT (' Interaction strength (units of kT) =',1PE12.5)
10    FORMAT (' NX =',I3,5X,' NY =',I3)
15    FORMAT (' number of thermalization sweeps =',I4)
20    FORMAT (' sweep frequency = ',I4,' group size =',I4)
C                               
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DISPLY(MUNIT,SPIN,MAG,ENERGY,ACCPT,IGROUP,ISWEEP)
C display spin configuration (spin=-1 is a blank, spin=1 is an X)
C and write out data for this sweep
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E8'
      INCLUDE 'IO.ALL'
C Input variables:
      INTEGER SPIN(MAXX,MAXY)    !spin configuration
C     all of the thermodynamic quant have 2 indices
C     (not all array elements are used, e.g. CHI(sweep,value))
C     first index is the level: sweep, group, or total
C     second index is the value: quantity, quant**2, or sigma**2
      REAL MAG(3,3)              !magnetization
      REAL ENERGY(3,3)           !energy
      REAL ACCPT                 !acceptance ratio
      INTEGER MUNIT              !unit we're writing to
      INTEGER ISWEEP,IGROUP      !sweep and group index
C Local variables:
      INTEGER SWEEP,GROUP,TOTAL  !which level of calculation
      INTEGER VALUE,SQUARE,SIGSQ !which quantity
      INTEGER IX,IY              !lattice indices
      CHARACTER*1 CSPIN(MAXX)    !spin as character data
      CHARACTER*80 BLNK          !blanks for centering in X
      DATA BLNK /' '/
      DATA SWEEP,GROUP,TOTAL/1,2,3/
      DATA VALUE,SQUARE,SIGSQ/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) THEN
        CALL CLEAR
        DO 20 IY=1,YCNTR                   !center output
           WRITE (OUNIT,*) ' '
20      CONTINUE
      END IF
      WRITE (MUNIT,11) IGROUP,ISWEEP,NSIZE,ACCPT,
     +          ENERGY(SWEEP,VALUE)/NSPIN,MAG(SWEEP,VALUE)/NSPIN
11    FORMAT (' group ',I3,', sweep ',I3,' out of ',I3,5X,
     +        ' accpt =',F5.3,'  Energy = ',F7.3,'  Mag =',F6.3)
C
      DO 100 IY=NY,1,-1                    !change +-1 to X and blank
         DO 50 IX=1,NX
            IF (SPIN(IX,IY) .EQ. 1) THEN
                CSPIN(IX)='X'
            ELSE                                                        
                CSPIN(IX)=' '
            END IF
50       CONTINUE    
C        write out a line at a time (no centering done for TUNIT)
         IF (MUNIT .EQ. TUNIT) THEN    
                  WRITE (TUNIT,16) (CSPIN(IX),IX=1,NX)
         ELSE
           IF (XSKIP) THEN
             WRITE (MUNIT,10) BLNK(1:XCNTR),(CSPIN(IX),IX=1,NX)
           ELSE 
             WRITE (MUNIT,15) BLNK(1:XCNTR),(CSPIN(IX),IX=1,NX)
           END IF
           IF (YSKIP) WRITE (MUNIT,*) ' '
         END IF
10       FORMAT (1X,A,100(A1,1X))                    
15       FORMAT (1X,A,100(A1))
16       FORMAT (1X,100(A1))
100   CONTINUE
      IF (MUNIT .EQ. OUNIT) CALL PAUSE('to continue...',0)
C
      RETURN
      END  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SWPOUT(MUNIT,MAG,ENERGY,ACCPT,IGROUP,ISWEEP,NLINES)
C and write out data for this sweep
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E8'
      INCLUDE 'IO.ALL'
C Input variables:
C     all of the thermodynamic quant have 2 indices
C     (not all array elements are used, e.g. CHI(sweep,value))
C     first index is the level: sweep, group, or total
C     second index is the value: quantity, quant**2, or sigma**2
      REAL MAG(3,3)              !magnetization
      REAL ENERGY(3,3)           !energy
      REAL ACCPT                 !acceptance ratio
      INTEGER MUNIT              !unit we're writing to
      INTEGER ISWEEP,IGROUP      !sweep and group index
      INTEGER NLINES             !lines written to terminal (I/O)
C Local variables:
      INTEGER SWEEP,GROUP,TOTAL  !which level of calculation
      INTEGER VALUE,SQUARE,SIGSQ !which quantity
      DATA SWEEP,GROUP,TOTAL/1,2,3/
      DATA VALUE,SQUARE,SIGSQ/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE (MUNIT,11) IGROUP,ISWEEP,NSIZE,ACCPT,
     +       ENERGY(SWEEP,VALUE)/NSPIN,MAG(SWEEP,VALUE)/NSPIN
11    FORMAT (7X,I3,7X,I3,'/',I3,7X,F5.3,7X,F9.5,7X,F9.3)
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+1
      RETURN
      END  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MAG,ENERGY,CB,CHI,E,ESIG1,ESIG2,
     +    M,MSIG1,MSIG2,SUS,SUSSIG,C,CSIG,IGROUP,MUNIT)
C write out averages and uncertaintes to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E8'
      INCLUDE 'IO.ALL'
C Input variables:
C     all of the thermodynamic quant have 2 indices
C     (not all array elements are used, e.g. CHI(sweep,value))
C     first index is the level: sweep, group, or total
C     second index is the value: quantity, quant**2, or sigma**2
      REAL MAG(3,3)              !magnetization
      REAL ENERGY(3,3)           !energy
      REAL CB(3,3)               !specific heat
      REAL CHI(3,3)              !susceptibility
      INTEGER IGROUP             !group index
      REAL M,MSIG1,MSIG2         !magnetization and uncertainties
      REAL E,ESIG1,ESIG2         !energy and uncertainties
      REAL SUS,SUSSIG            !susceptibility and uncertainty
      REAL C,CSIG                !specific heat and uncertainty
      INTEGER MUNIT              !which unit number
C Local variables:
      INTEGER SWEEP,GROUP,TOTAL  !which level of calculation
      INTEGER VALUE,SQUARE,SIGSQ !which quantity
      DATA SWEEP,GROUP,TOTAL/1,2,3/
      DATA VALUE,SQUARE,SIGSQ/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE (MUNIT,30) IGROUP,NGROUP
      WRITE (MUNIT,32)
      WRITE (MUNIT,33)
      WRITE (MUNIT,35) 
     +  ENERGY(GROUP,VALUE)/NSPIN,SQRT(ENERGY(GROUP,SIGSQ))/NSPIN,
     +  MAG(GROUP,VALUE)/NSPIN,SQRT(MAG(GROUP,SIGSQ))/NSPIN,
     +  CHI(GROUP,VALUE)/NSPIN,CB(GROUP,VALUE)/NSPIN
      WRITE (MUNIT,40) E/NSPIN,ESIG1/NSPIN,ESIG2/NSPIN,
     +  M/NSPIN,MSIG1/NSPIN,MSIG2/NSPIN,
     +  SUS/NSPIN,SUSSIG/NSPIN,C/NSPIN,CSIG/NSPIN
      WRITE (MUNIT,*) ' '
      IF ((MUNIT .EQ. OUNIT) .AND. (.NOT. TERSE)) 
     +        CALL PAUSE(' to continue...',1)
C
30    FORMAT ('  Group ',I3,' (out of ',I4,') averages')
32    FORMAT (14X,'Energy',13X,'Magnetization',5X,'Susceptibility',
     +        2X,'Specific Heat')
33    FORMAT (14X,'------',13X,'-------------',5X,'--------------',
     +        2X,'-------------')
35    FORMAT ('  group ',2(1X,F7.3,'+-',F5.3,6X),
     +        2(2X,F6.3,7X))
40    FORMAT ('  total ',2(1X,F7.3,'+-',F5.3,'/',F5.3),
     +        2(1X,F6.3,'+-',F6.3))
C
      RETURN
      END                                   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TITLES(MUNIT,NLINES)
C write out text data titles
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input/Output variables:
      INTEGER MUNIT                   !output unit (input)
      INTEGER NLINES                  !number of lines written (I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) CALL CLEAR
C
      WRITE (MUNIT,10)
      WRITE (MUNIT,11)
10    FORMAT (6X,'Group',3X,'Sweep/Out of',5X,'Accpt',9X,'Energy',
     +       6X,'Magnetization')      
11    FORMAT (6X,'-----',3X,'------------',5X,'-----',9X,'------',
     +       6X,'-------------')      
C
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+2
C
      RETURN                                                            
      END
                                                                    
