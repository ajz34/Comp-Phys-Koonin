CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM EXMPL3
C     Example 3: Bound states in a one-dimensional potential
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company, Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop; execute once for each set of param
         CALL PARAM       !get input parameters
         CALL ARCHON      !solve time-independent Schroedinger equation
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON
C solves the time-independent Schroedinger equation for 1-dimensional
C potential; allows for more than one eigenvalue search
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E3'
C Local variables:
      REAL ENERGY(MAXLEV)            !array of eigenvalues
      INTEGER LEFTTP(MAXLEV),RGHTTP(MAXLEV)   !array of turning points
      INTEGER NFOUND                 !number of levels found
      REAL PSI(MAXLAT)               !wave function
      INTEGER ITER                   !continue finding levels?
      INTEGER NODES                  !number of nodes
      REAL E                         !trial eigenenergy 
      INTEGER LTP,RTP                !classical turning points (lattice)
      REAL DE,DESAVE                 !step in energy search
      INTEGER NSTP                   !number of steps taken in search
      INTEGER SCREEN                 !send to terminal
      INTEGER PAPER                  !make a hardcopy
      INTEGER FILE                   !send to a file
C Functions:
      REAL GETFLT                    !obtain real input
      INTEGER YESNO                  !obtain yes or no input  
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     output summary of parameters                    
      IF (TTERM) CALL PRMOUT(OUNIT) 
      IF (TFILE) CALL PRMOUT(TUNIT) 
      IF (GFILE) CALL PRMOUT(GUNIT)
C
      ITER=1                          !initialize searching variables
      NFOUND=0
      E=EBOTTM                        !reasonable starting values
      DE=DEBOT
C                                                    
5     IF (ITER .EQ . 1) THEN          !loop over energy levels
C                                      
          E=GETFLT(E,EMIN,EMAX,'Enter energy')     !get starting values
          DE=GETFLT(DE,DEMIN,DEMAX,'Enter step in energy search')       
          DESAVE=DE                   !save for next level
C                                                            
          IF (TTERM) CALL CTITL(OUNIT)!output titles
          IF (TFILE) CALL CTITL(TUNIT)
C
          CALL SEARCH(E,DE,LTP,RTP,NODES,NSTP,PSI) !search for energy
C
C         prompt for new E and DE if too many steps taken in search
          IF (NSTP .GE. MAXSTP) THEN
              WRITE (OUNIT,*) ' '
              WRITE (OUNIT,*) ' Eigenstate not yet found'
              ITER=YESNO(1,
     +          'Do you want to try again with new E or DE?')
          ELSE
C         or output eigenvalue information
            IF (TTERM) CALL EOUT(E,LTP,RTP,NODES,OUNIT)
            IF (TFILE) CALL EOUT(E,LTP,RTP,NODES,TUNIT)
            IF (GFILE) CALL EOUT(E,LTP,RTP,NODES,GUNIT)
C           save energy and turning points
            IF (NFOUND .EQ. MAXLEV) THEN
               NFOUND=0
               WRITE (OUNIT,*) ' Storage capacity for levels exceeded'
               WRITE (OUNIT,*) ' Old information will be overwritten'
            END IF                        
            NFOUND=NFOUND+1
            ENERGY(NFOUND)=E
            LEFTTP(NFOUND)=LTP
            RGHTTP(NFOUND)=RTP
C
            IF (TTERM) CALL PAUSE('to continue...',1)
            IF (TTERM) CALL CLEAR
C
            IF (GTERM) THEN             !graphics output
               CALL GRFOUT(SCREEN,ENERGY,LEFTTP,RGHTTP,NFOUND,PSI)
               CALL CLEAR
            END IF
            IF (GFILE) CALL GRFOUT(FILE,ENERGY,LEFTTP,RGHTTP,NFOUND,PSI)
            IF (GHRDCP) 
     +          CALL GRFOUT(PAPER,ENERGY,LEFTTP,RGHTTP,NFOUND,PSI)
C
C           check if another level is to be found
            WRITE (OUNIT,*) ' '
            WRITE (OUNIT,*) ' '
            ITER=YESNO(1,'Do you want to find another level?')
C
          END IF
C
          E=E+DESAVE          !reasonable starting values for next level
          DE=DESAVE
C
      GOTO 5
      END IF        
C              
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SEARCH(E,DE,LTP,RTP,NODES,NSTP,PSI)
C search for one eigenvalue beginning with energy E and step DE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E3'
C Passed variables:
      INTEGER NODES          !number of nodes (output)
      REAL E                 !trial eigenenergy, true eigenenergy (I/O)
      INTEGER LTP,RTP        !classical turning points (lattice)(output)
      REAL DE                !step in energy search (input)
      INTEGER NSTP           !number of steps taken in search (output)
      REAL PSI(MAXLAT)       !wave function (output)
C Local variables:
      LOGICAL SECANT         !are we doing secant search yet?
      REAL F,FOLD,FSTART     !values of the mismatch
      REAL EOLD              !last value of the energy
      INTEGER NCROSS         !classically forbidden regions
      LOGICAL FIRST          !signal first step in search
      INTEGER NLINES         !number of lines sent to screen
C Functions:
      REAL GETFLT            !obtain real input
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FIRST=.TRUE.           !first call to NUMERV for this level
      SECANT=.FALSE.         !start off with simple search
      F=TOLF                 !dummy value to get into loop
      NLINES=2
      NSTP=0                  
C
C    search for eigenenergy until F is small enough or too many
C    steps are taken
10    IF ((ABS(F) .GE. TOLF) .AND. (NSTP .LT. MAXSTP)) THEN
          NSTP=NSTP+1
C         integ Schroedinger equation; find discontinuity in derivative
          CALL NUMERV (F,NODES,E,LTP,RTP,NCROSS,FIRST,PSI)
          IF (FIRST .EQV. .TRUE.) THEN
              FSTART=F
              FIRST=.FALSE.
          END IF
          IF (F*FSTART  .LE. 0.) SECANT=.TRUE.    
          IF ((SECANT .EQV. .TRUE.) .AND. (F .NE. FOLD)) THEN
              DE=-F*(E-EOLD)/(F-FOLD)
          END IF
C
          IF (TTERM) CALL TXTOUT(OUNIT,E,DE,F,NODES,NCROSS,NLINES)
          IF (TFILE) CALL TXTOUT(TUNIT,E,DE,F,NODES,NCROSS,NLINES)
C
          EOLD=E             !update values for next step
          FOLD=F                                                     
          E=E+DE
C                                
C         keep energy greater than VMIN=-1
          IF (E .LT. -1.) THEN
              WRITE (OUNIT,*) ' Energy is less than -1.'
              E=GETFLT(-.99,EMIN,EMAX,' Enter new energy')
              DE=GETFLT(DE,DEMIN,DEMAX,'Enter step in energy search')
          END IF
C
C         end search if DE is smaller than machine precision         
          IF ((ABS(DE) .LT. DEMIN) .AND. (ABS(F) .GT. TOLF)) THEN
             WRITE (OUNIT,20)
20           FORMAT('0',21X,'Delta E is too small, search must stop')
             F=TOLF/2             !dummy value to end search
          END IF
C
      GOTO 10  
      END IF 
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NUMERV (F,NODES,E,LTP,RTP,NCROSS,FIRST,PSI)
C subroutine to integrate time independent Schroedinger equation
C with energy=E
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E3'
C Passed variables:
      REAL F                     !size of mismatch (output)
      INTEGER NODES              !number of nodes (output)
      REAL E                     !trial eigenvalue (input)
      INTEGER LTP,RTP            !classical turning points (output)
      INTEGER NCROSS             !classically forbidden regions (output)
      LOGICAL FIRST              !signal first step in search (output)
      REAL PSI(MAXLAT)           !wave function (output)
C Local variables:
      INTEGER IMATCH             !matching lattice point
      REAL C                     !useful constant
      INTEGER IX,KX              !X indices
      REAL KI,KIM1,KIP1          !terms in Numerov algorithm
      REAL NORM                  !norm of wave function
      REAL PMMTCH,PPMTCH         !PSI at match point
      LOGICAL FLIP               !do we need to flip sign of PSI
      Real PSIMAX                !max value of wave function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (FIRST .EQV. .TRUE.) IMATCH=0   !signals search for IMATCH
      C=(DX*DX/12)*GAMMA*GAMMA           !evaluate constant
C                                                                        
C     find IMATCH once per energy level by looking for entry
C     into classically forbidden regions
      IX=1
      KIM1=C*(E-V(IX))
1     IF (IMATCH .EQ. 0) THEN
         IX=IX+1
         KI=C*(E-V(IX))    !neg value of K indicates class frbdn region
         IF ((KI*KIM1 .LT. 0) .AND. (KIM1 .GT. 0)) IMATCH=IX
C        if other procedure fails to find IMATCH, set IMATCH=NPTS-10
         IF (IX .EQ. NPTS-10) IMATCH=IX
         KIM1=KI
      GOTO 1
      END IF
C
      PSI(1)=0                     !left hand side bound. cond.
      PSI(2)=9.999999E-10          
      KIM1=C*(E-V(1))              !initial K*K values
      KI=C*(E-V(2))
C
C     Numerov algorithm; S=0; integrate until enter class. frbdn. region
      DO 10 IX=2,IMATCH
         KIP1=C*(E-V(IX+1))
         PSI(IX+1)=(PSI(IX)*(2.-10.*KI)-PSI(IX-1)*(1+KIM1))/(1+KIP1)
         KIM1=KI                  !roll values of K*K
         KI=KIP1
C
C        if PSI grows too large rescale all previous points
         IF (ABS(PSI(IX+1)) .GT. (1.0E+10)) THEN
             DO 20 KX=1,IX+1
                PSI(KX)=PSI(KX)*9.999999E-06
20           CONTINUE
         END IF
10    CONTINUE
      PMMTCH=PSI(IMATCH)           !save value for normalization
C
      PSI(NPTS)=0                  !rhs boundary conditions
      PSI(NPTS-1)=9.999999E-10
      KIP1=C*(E-V(NPTS))           !initial K*K values
      KI=C*(E-V(NPTS-1))
C
C     Numerov algorithm, S=0;integrate from rhs to IMATCH
      DO 30 IX=NPTS-1,IMATCH+1,-1   
           KIM1=C*(E-V(IX-1))
           PSI(IX-1)=(PSI(IX)*(2.-10.*KI)-PSI(IX+1)*(1+KIP1))/(1+KIM1)
           KIP1=KI                 !roll values of K*K
           KI=KIM1
C 
C          if PSI grows too large rescale all previous points
           IF (ABS(PSI(IX-1)) .GT. (1.0E+10)) THEN 
              DO 40 KX=NPTS-1,IX-1,-1
                   PSI(KX)=PSI(KX)*9.999999E-06
40            CONTINUE
           END IF
30    CONTINUE                                    
C
      KIM1=C*(E-V(IMATCH-1))       !finds values needed for log deriv
      PPMTCH=(PSI(IMATCH)*(2-10*KI)-PSI(IMATCH+1)*(1+KIP1))/(1+KIM1)
C
      NORM=PMMTCH/PSI(IMATCH)      !norm PSI right to PSI left
      DO 5 IX=IMATCH,NPTS                                     
           PSI(IX)=PSI(IX)*NORM
5     CONTINUE
      PPMTCH=PPMTCH*NORM
C                                                       
C     find nodes, turning points, entry into classically forb regions
C     and determine if overall sign must be flipped
      CALL DETAIL(IMATCH,NODES,E,LTP,RTP,NCROSS,PMMTCH,FLIP,FIRST,PSI)
C
C     find maximum PSI value and flip sign if necessary
      PSIMAX=ABS(PSI(1))           
      DO 7 IX=1,NPTS
           IF (ABS(PSI(IX)) .GT. PSIMAX) PSIMAX=ABS(PSI(IX))
           IF (FLIP .EQV. .TRUE.) PSI(IX)=-PSI(IX)
7     CONTINUE
      IF (FLIP) PPMTCH=-PPMTCH
C
      F=(PSI(IMATCH-1)-PPMTCH)/PSIMAX   !evaluate matching condition
C
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE 
     +    DETAIL(IMATCH,NODES,E,LTP,RTP,NCROSS,PMMTCH,FLIP,FIRST,PSI)
C subroutine to calculates nodes, turning points, entry
C into classically forbidden regions, and whether or not the
C overall sign of PSI must be flipped to make F continuous function of E
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E3'
C Passed variables:
      INTEGER IMATCH            !matching lattice point (input)
      INTEGER NODES             !number of nodes (output)
      REAL E                    !trial eigenvalue (input)
      INTEGER LTP,RTP           !classical turning points (output)
      INTEGER NCROSS            !classically forbidden regions (output)
      REAL PMMTCH               !PSI at match point (input)
      LOGICAL FIRST             !signal first step in search(input)
      LOGICAL FLIP              !do we need to flip sign of PSI (output)
      REAL PSI(MAXLAT)          !wave function (input)
C Local variables:
      LOGICAL LSAME             !does PSI left retain its sign?
      INTEGER IX                !X index
      REAL K,KLAST              !wavenumbers 
      INTEGER NLOLD,NLFT        !left nodes
      INTEGER NROLD,NRT         !right nodes
      REAL PSIR,PSIL            !sign of PSI which must be constant
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     find number of nodes on each side
      NLOLD=NLFT                !save old values first
      NROLD=NRT
      NLFT=0                    !zero sums                 
      NRT=0
      DO 8 IX=2,IMATCH-1        !left side
          IF (PSI(IX)*PSI(IX-1) .LT. 0) NLFT=NLFT+1
8     CONTINUE
      IF (PSI(IMATCH-1)*PMMTCH .LT. 0) NLFT=NLFT+1
      DO 10 IX=IMATCH+1,NPTS        !right side
          IF (PSI(IX)*PSI(IX-1) .LT. 0) NRT=NRT+1             
10    CONTINUE
C
C     check for change in node number
C      if one side of the wave function has gained a node,
C      that side must maintain the same sign between steps
C      in order that F be a continuous function of Energy.
C      if no new node has appeared, keep same sign same
C      as last time
      IF ((NLOLD .NE. NLFT) .AND. (NRT .EQ. NROLD)) THEN
         IF (LSAME .EQV. .FALSE.) LSAME=.TRUE.
      ELSE IF ((NROLD .NE. NRT) .AND. (NLOLD .EQ. NLFT)) THEN
         IF (LSAME .EQV. .TRUE.) LSAME=.FALSE.
      END IF
C
      IF (FIRST .EQV. .TRUE.) LSAME=.TRUE.    !initialize variables
      IF (FIRST .EQV. .TRUE.) PSIL=1.
C                                                           
C     now, finally, determine if the sign needs to be flipped
      FLIP=.FALSE.
      IF ((LSAME) .AND. (PSIL*PSI(2) .LT. 0.)) FLIP=.TRUE.
      IF (( .NOT. LSAME ) .AND. (PSIR*PSI(NPTS-1) .LT. 0.))
     +    FLIP=.TRUE.  
C
C     save values for next step
      PSIL=PSI(2)
      PSIR=PSI(NPTS-1)
      IF (FLIP) PSIL=-PSIL
      IF (FLIP) PSIR=-PSIR
C
      NODES=NLFT+NRT             !total number of nodes
C 
C     find leftmost turning point
      LTP=0
      IX=0
110   IF (LTP .EQ. 0) THEN
         IX=IX+1
         IF ((E-V(IX)) .GE. 0) LTP=IX
      GOTO 110
      END IF
C
C     find rightmost turning point
      RTP=0
      IX=NPTS+1
120   IF (RTP .EQ. 0) THEN
         IX=IX-1
         IF ((E-V(IX)) .GE. 0) RTP=IX
      GOTO 120
      END IF
C
C     do we integrate into classically forbidden regions?
      KLAST=(E-V(LTP))
      NCROSS=0
      DO 90 IX=LTP+1,RTP-1
         K=(E-V(IX))
         IF ((K*KLAST .LE. 0) .AND. (K .LT. 0)) NCROSS=NCROSS+1
         KLAST=K
90    CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE POTNTL
C sets up array for current potential;
C all potentials have a minimum value of -1 and max of 1;
C if you change the potential, also change XMIN and XMAX in INIT,
C and estimate DEBOT and EBOTTM (one fifth of the expected spacing
C at the bottom of the well, and two DEBOT's below the expected ground
C state energy).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E3'
C Local variables:
      INTEGER ILEFT,IWID,IRIGHT     !lattice values for width, left
      REAL VSCALE                   !scaling factor for renorm of parab
      INTEGER IX                    !labels lattice points
      REAL VMIN                     !minimum value of potential 
      REAL CURVE                    !curvature for smoothing bump
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     define limits on X and DX
      IF (POT .EQ. SQUARE) THEN
         VXMIN=XMIN1
         VXMAX=XMAX1
      ELSE IF (POT .EQ. PARAB) THEN                  
         VXMIN=XMIN2
         VXMAX=XMAX2
      ELSE IF (POT .EQ. LENRD) THEN
         VXMIN=XMIN3
         VXMAX=XMAX3
      END IF           
C
      DX=(VXMAX-VXMIN)/(NPTS-1)
C
C     setup X (space coordinate) array
      DO 5 IX=1,NPTS
         X(IX)=VXMIN+(IX-1)*DX
5     CONTINUE             
C                                
C     setup V (potential) array
      IF (POT .EQ. SQUARE) THEN             !square well
         DO 10 IX=1,NPTS
            V(IX)=-1.
10       CONTINUE
         VMAX=1.
         DEBOT=.2*(3.1415926/(4*GAMMA))**2!energy spacing at well bottom
         IF (BUMP .EQ. 1) THEN            !bump in well bottom
            ILEFT=NINT((LEFT-XMIN1)/DX)+1 !locate bump edges
            IRIGHT=NPTS+NINT((RIGHT-XMAX1)/DX)
            RIGHT=X(IRIGHT)
            LEFT=X(ILEFT)
            CURVE=-.025/((RIGHT-LEFT)/2.)**2
C           top of bump is curved to avoid discontinuities
            DO 20 IX=ILEFT,IRIGHT
               V(IX)=-1.+HEIGHT*(1.+CURVE*(X(IX)-(RIGHT+LEFT)/2.)**2)
20          CONTINUE
         END IF
C
      ELSE IF (POT .EQ. PARAB) THEN    !parabolic well             
         DO 30 IX=1,NPTS
            V(IX)=-(1.-.5*X(IX)*X(IX))
30       CONTINUE
         VMAX=1. 
         DEBOT=.2/(SQRT(2.)*GAMMA)     !energy spacing at well bottom
         IF (FLAT .EQ. 1) THEN
C           the parabolic potential is flattened symmetrically about x=0
            IWID=INT(PWID/DX)
            IF (MOD(IWID,2) .NE. MOD(NPTS,2)) IWID=IWID-1  
            ILEFT=(NPTS-IWID)/2+1
            VMIN=V(ILEFT)
            DO 40 IX=ILEFT,ILEFT+IWID-1
               V(IX)=VMIN
40          CONTINUE
C
            VSCALE=2./(V(1)-VMIN)    !normalize so that min value is -1
            DO 50 IX=1,NPTS
               V(IX)=VSCALE*(V(IX)-VMIN)-1.
50          CONTINUE
         END IF
C 
      ELSE IF (POT .EQ. LENRD) THEN    !Lennard-Jones potential
         DO 60 IX=1,NPTS  
            V(IX)=4*(X(IX)**(-12.)-X(IX)**(-6))
60       CONTINUE
         VMAX=V(1)
         DEBOT=1.0692/GAMMA !energy spacing at well bottom (see Ex. 1.8)
C
      END IF           
C 
C     educated guesses for ground state energy and energy spacing;
C     with these values, it will take three steps for F to change sign
      EBOTTM=-1.+2.5*DEBOT
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
      INCLUDE 'PARAM.E3'
C Local parameters:
      CHARACTER*80 DESCRP        !program description
      DIMENSION DESCRP(22)
      INTEGER NHEAD,NTEXT,NGRF   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL SETUP                !get environment parameters
C                                               
C     display header screen     
      DESCRP(1)= 'EXAMPLE 3'
      DESCRP(2)= 'Bound states in a one-dimensional potential'
      NHEAD=2
C
C     text output description
      DESCRP(3)= 'during search: energy, delta energy, discontinuity,'
      DESCRP(4)='nodes, and entrance into classically forbidden regions'
      DESCRP(5)= 'after search: energy, classical turning points,'
     +  //' number of nodes'
      NTEXT=3
C
C     graphics output description
      DESCRP(6)= 'potential, eigenvalues, and eigenfunctions' 
      NGRF=1
C                          
C     setup constant values for the potentials 
      !all potentials have VMIN=-1, VMAX=1
      XMIN1=-2.
      XMAX1=2.
      XMIN2=-2.
      XMAX2=2.
      XMIN3=(2.*(SQRT(2.)-1.))**(.1666666)      
      XMAX3=1.9 
C
C     reasonable limits on energy and energy increment
      EMIN=-.99999
      EMAX=100.
      DEMIN=5.E-06
      DEMAX=1.
C
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRF)
C
      CALL MENU                !setup constant part of header
C
      MTYPE(13)=FLOAT
      MPRMPT(13)= 'Enter Gamma=sqrt[2m(a**2)V/hbar**2] (dimensionless)'
      MTAG(13)= 'Gamma=sqrt[2m(a**2)V/hbar**2] (dimensionless)'
      MLOLIM(13)=0.001
      MHILIM(13)=150.
      MREALS(13)=30.0
C                
      MTYPE(14)=TITLE
      MPRMPT(14)= 'POTENTIAL FUNCTION MENU'
      MLOLIM(14)=2.
      MHILIM(14)=1.
C
      MTYPE(15)=MTITLE                   
      MPRMPT(15)='1) Square-well potential'
      MLOLIM(15)=0.
      MHILIM(15)=0.
C             
      MTYPE(16)=MTITLE
      MPRMPT(16)='2) Parabolic-well potential'
      MLOLIM(16)=0.
      MHILIM(16)=0.
C
      MTYPE(17)=MTITLE
      MPRMPT(17)='3) Lennard-Jones potential'
      MLOLIM(17)=0.
      MHILIM(17)=1.
C
      MTYPE(18)=MCHOIC
      MPRMPT(18)= 'Make menu choice and press Return'
      MTAG(18)='19 24 26'
      MLOLIM(18)=1.
      MHILIM(18)=3.
      MINTS(18)=3
      MREALS(18)=3.
C 
      MTYPE(19)=NOSKIP
      MPRMPT(19)= 'Do you want a bump in the square well?'
      MTAG(19)= 'Square Well Bump:'
      MREALS(19)=35.
      MINTS(19)=0
C              
      MTYPE(20)=FLOAT
      MPRMPT(20)= 'Enter the left hand edge of bump (scaled units)'
      MTAG(20)= 'Left hand edge of bump (scaled units)'
      MLOLIM(20)=XMIN1
      MHILIM(20)=XMAX1                                 
      MREALS(20)=-.1
C               
      MTYPE(21)=FLOAT
      MPRMPT(21)= 'Enter the right hand edge of bump (scaled units)'
      MTAG(21)= 'Right hand edge of bump (scaled units)'
      MLOLIM(21)=XMIN1
      MHILIM(21)=XMAX1 
      MREALS(21)=.1
C              
      MTYPE(22)=FLOAT
      MPRMPT(22)= 'Enter the height of the bump (scaled units)'
      MTAG(22)= 'Height of the bump (scaled units)'
      MLOLIM(22)=0.
      MHILIM(22)=20
      MREALS(22)=1.
C                
      MTYPE(23)=SKIP
      MREALS(23)=35.
C
      MTYPE(24)=NOSKIP
      MPRMPT(24)= 'Do you want a flattened bottom in parabolic-well'
      MTAG(24)= 'Parabolic Flattened area:'
      MREALS(24)=35.
      MINTS(24)=0
C                        
      MTYPE(25)=FLOAT
      MPRMPT(25)= 'Enter length of the flattened area (scaled units)'
      MTAG(25)= 'Length of the flattened area (scaled units)'
      MLOLIM(25)=0
      MHILIM(25)=XMAX2-XMIN2
      MREALS(25)=.1
C
      MTYPE(26)=SKIP
      MREALS(26)=35.
C                
      MTYPE(38)=FLOAT
      MPRMPT(38)= 'Enter the matching tolerance'
      MTAG(38)= 'Matching tolerance'
      MLOLIM(38)=5.E-06
      MHILIM(38)=1.000
      MREALS(38)=.00005     
C               
      MTYPE(39)=NUM
      MPRMPT(39)= 'Enter number of lattice points'
      MTAG(39)= 'Number of lattice points'
      MLOLIM(39)=50
      MHILIM(39)=MAXLAT
      MINTS(39)=160
C                  
      MTYPE(40)=NUM
      MPRMPT(40)= 'Enter maximum number of steps in search'
      MTAG(40)= 'Maximum number of search steps'
      MLOLIM(40)=10
      MHILIM(40)=300
      MINTS(40)=50
C
      MTYPE(41)=SKIP
      MREALS(41)=60.
C
      MSTRNG(MINTS(75))= 'exmpl3.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C            
      MSTRNG(MINTS(86))= 'exmpl3.grf'
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
      INCLUDE 'PARAM.E3'
C Local variables:
C map between menu items and parameters
      INTEGER IGAMMA,ITOLR,INPTS,IMXSTP 
      INTEGER IPOT,IBUMP,ILEFT,IHEIGH,IRIGHT,IFLAT,IPWID
      PARAMETER (IGAMMA = 13 )
      PARAMETER (IPOT   = 18 )
      PARAMETER (IBUMP  = 19 )
      PARAMETER (ILEFT  = 20 )
      PARAMETER (IRIGHT = 21 )
      PARAMETER (IHEIGH = 22 )
      PARAMETER (IFLAT  = 24 )
      PARAMETER (IPWID  = 25 )
      PARAMETER (ITOLR  = 38 )
      PARAMETER (INPTS  = 39 )
      PARAMETER (IMXSTP = 40 )
C Function:
      LOGICAL LOGCVT          !converts 1 and 0 to true and false 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get input from terminal
      CALL CLEAR                  
      CALL ASK(1,ISTOP)
      CALL CLEAR
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
C     physical and numerical:
      GAMMA=MREALS(IGAMMA)                                       
      TOLF=MREALS(ITOLR)
      NPTS=MINTS(INPTS)
      MAXSTP=MINTS(IMXSTP)
C      
C     potential parameters
      POT=MREALS(IPOT)
      BUMP=MINTS(IBUMP)
      LEFT=MREALS(ILEFT)
      RIGHT=MREALS(IRIGHT)
      HEIGHT=MREALS(IHEIGH)
      FLAT=MINTS(IFLAT)
      PWID=MREALS(IPWID)
C
      CALL POTNTL  !setup potential array
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
      !files may have been renamed
      MSTRNG(MINTS(ITNAME))=TNAME
      MSTRNG(MINTS(IGNAME))=GNAME
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE PRMOUT(MUNIT) 
C outputs parameter summary to the specified unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'IO.ALL'
       INCLUDE 'PARAM.E3'
C Passed variables:               
       INTEGER MUNIT            !unit number for output
       REAL VAR(4)              !initial values of coord and momenta
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (MUNIT .EQ. OUNIT) CALL CLEAR
C 
       WRITE (MUNIT,2)
       WRITE (MUNIT,4)
       WRITE (MUNIT,6) GAMMA
       WRITE (MUNIT,8) TOLF
       WRITE (MUNIT,12) NPTS
C
       IF (POT .EQ. SQUARE) THEN
          WRITE (MUNIT,14)
          IF (BUMP .EQ. 1) WRITE (MUNIT,16) LEFT,RIGHT,HEIGHT
       ELSE IF (POT .EQ. PARAB) THEN
          WRITE (MUNIT,18)
          IF (FLAT .EQ. 1) WRITE (MUNIT,20) -PWID,PWID
       ELSE IF (POT .EQ. LENRD) THEN
          WRITE (MUNIT,22)
       END IF
C
       WRITE (MUNIT,24)VXMIN,VXMAX,DX
       WRITE (MUNIT,26)
       WRITE (MUNIT,2)
C
2      FORMAT (' ')
4      FORMAT (' Output from Example 3: Time Independent Schroedinger',
     +         ' Equation')
6      FORMAT (' Gamma =',F9.3)
8      FORMAT (' Matching tolerance = ', 1PE10.3)
12     FORMAT (' Number of lattice points =', I5)
14     FORMAT (' Square Well potential')
16     FORMAT (' with bump from ',F6.3,' to ',F6.3,' of height =',F8.3)
18     FORMAT (' Parabolic potential')
20     FORMAT (' with flat bottom from ',F6.3,' to ',F6.3)
22     FORMAT (' Lennard-Jones potential')
24     FORMAT (' Xmin=',F6.3,' Xmax=',F6.3,' Dx=',1PE10.3)
26     FORMAT (' Energies and lengths are in scaled units')
C
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE CTITL(MUNIT)
C writes title for text output to requested device
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
      INTEGER MUNIT                 !unit to which we are writing
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE (MUNIT,30)
      WRITE (MUNIT,35)
C
30    FORMAT 
     +     (10X,'Energy',14X,'De',16X,'F',10x,'Nodes',5x,'Frbddn')
35    FORMAT 
     +     (10X,'------',14X,'--',16X,'-',10X,'-----',5X,'------')
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,E,DE,F,NODES,NCROSS,NLINES)
C writes results for one level to the requested unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E3'
C Input variables:
      INTEGER MUNIT                !unit to which we are writing
      REAL F                       !values of the mismatch
      INTEGER NODES                !number of nodes
      REAL E                       !trial eigenenergy 
      REAL DE                      !step in energy search
      INTEGER NCROSS               !classically forbidden regions
      CHARACTER*3 FORBDN           !have we integrated into forb regions
      INTEGER NLINES               !number of lines sent to screen (I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     if screen is full, clear screen and retype headings          
      IF ((MOD(NLINES,TRMLIN-4) .EQ. 0) 
     +                          .AND. (MUNIT .EQ. OUNIT)) THEN
         CALL PAUSE('to continue...',1)
         CALL CLEAR
         CALL CTITL(OUNIT)
         NLINES=2
      END IF
C
      IF (NCROSS .GT. 0) THEN 
         FORBDN='Yes'
      ELSE
         FORBDN='No '
      END IF
      WRITE (MUNIT,20) E,DE,F,NODES,FORBDN
      NLINES=NLINES+1
C                 
20    FORMAT (8X,F10.7,7X,1PE10.3,7X,1PE10.3,8X,I2,8X,A3)
C      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE EOUT(E,LTP,RTP,NODES,MUNIT)
C output information about eigenstate
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E3'
C Passed variables:
      INTEGER MUNIT                 !unit to which we are writing
      INTEGER NODES                 !number of nodes
      REAL E                        !trial eigenenergy 
      INTEGER LTP,RTP               !classical turning points (lattice)
      REAL XLTP,XRTP                !classical turning points
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      XLTP=VXMIN+(LTP-1)*DX         !turning points
      XRTP=VXMIN+(RTP-1)*DX
      WRITE (MUNIT,500) E, NODES
      WRITE (MUNIT,510) XLTP,XRTP
C
500   FORMAT ('0',20X,'Eigenvalue = ',F10.7,' with ',I2,' nodes.')
510   FORMAT (16X,'and classical turning points ',F7.4,' and ',F7.4)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFOUT(DEVICE,ENERGY,LEFTTP,RGHTTP,NFOUND,PSI)
C outputs potential and eigenvalues (plot 1) and eigenfunctions (plot 2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables                                                   
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E3'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE               !which device is being used?
      REAL ENERGY(MAXLEV)          !array of eigenvalues
      INTEGER LEFTTP(MAXLEV),RGHTTP(MAXLEV)   !array of turning points
      INTEGER NFOUND               !number of levels found
      REAL PSI(MAXLAT)             !wave function
C Local variables
      INTEGER IX                   !index of lattice
      INTEGER ILEVEL               !index of eigenvalue
      REAL XE(2),GRAPHE(2)         !arrays to graph eigenvalue
      CHARACTER*9 CGAMMA,CE        !GAMMA, ENERGY as a character string
      INTEGER LEN                  !length of string
      REAL NORM,MAXPSI             !norm and max value for PSI
      INTEGER SCREEN               !send to terminal
      INTEGER PAPER                !make a hardcopy
      INTEGER FILE                 !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
      NORM=0.
      DO 10 IX=1,NPTS               !normalize wave function 
        NORM=NORM+PSI(IX)**2 
10    CONTINUE
      NORM=SQRT(NORM*DX)
      MAXPSI=0.
      DO 30 IX=1,NPTS
         PSI(IX)=PSI(IX)/NORM
         IF (ABS(PSI(IX)) .GT. MAXPSI) MAXPSI = ABS(PSI(IX))
30    CONTINUE
C
      IF (DEVICE .NE. FILE) THEN
          NPLOT=2                        !how many plots?
C                                
C         parameters that are the same for both plots
          XMIN=VXMIN                     !axis parameters
          Y0VAL=VXMIN
          XMAX=VXMAX
          NPOINT=NPTS
          LABEL(1)= 'x (scaled units)'
C
          ILINE=1                        !line and symbol styles
          ISYM=1
          IFREQ=0
          NXTICK=5
          NYTICK=5
C                                                             
          CALL CONVRT(GAMMA,CGAMMA,LEN)
          TITLE='Solutions to Schroedinger''s Equation, Gamma='//CGAMMA
C
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
C                                                
C         first plot: potential and eigenvalues                            
          IPLOT=1
          YMIN=-1.
          YMAX=VMAX
          X0VAL=0.
          INFO=' '
          LABEL(2)= 'potential and eigenvalues (scaled units)'
          CALL LNLNAX                             
          CALL XYPLOT (X,V)                    !plot potential
C
          ILINE=3                              !plot eigenvalues
          NPOINT=2
          DO 20 ILEVEL=1,NFOUND
            XE(1)=VXMIN+(LEFTTP(ILEVEL)-1)*DX                          
            XE(2)=VXMIN+(RGHTTP(ILEVEL)-1)*DX
            GRAPHE(1)=ENERGY(ILEVEL)
            GRAPHE(2)=ENERGY(ILEVEL)
            CALL XYPLOT(XE,GRAPHE)
20        CONTINUE
C
C         second plot:  wave function
          IPLOT=2
          ILINE=1
          NPOINT=NPTS
          YMIN=-MAXPSI
          YMAX=+MAXPSI
          X0VAL=0.
          CALL CONVRT(ENERGY(NFOUND),CE,LEN)
          INFO='Eigenvalue = '//CE
          LABEL(2)= 'normalized wave function'
          CALL LNLNAX
          CALL XYPLOT (X,PSI)
C
      ELSE
C         output to file
          WRITE (GUNIT,80)
          WRITE (GUNIT,85)
          WRITE (GUNIT,70) (X(IX),V(IX),PSI(IX),IX=1,NPTS)
      END IF                 
C 
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)  !end graphing package
      IF (DEVICE .EQ. SCREEN) CALL TMODE        !switch to text mode
C
100   FORMAT (/,' Patience, please; output going to a file.')
80    FORMAT ('0',13X,'Position',14X,'Potential',13X,'Wave function')
85    FORMAT (13X,'--------',14X,'---------',13X,'------------')
70    FORMAT (3(11X,1PE12.5))
C
      RETURN
      END
                                                                        
