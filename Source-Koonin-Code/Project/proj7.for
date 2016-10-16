CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM PROJ7
C     Project 7: The Brusselator in two dimensions
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT          !display header screen, setup parameters
5     CONTINUE           !main loop/ execute once for each set of param
        CALL PARAM       !get input from screen
        CALL ARCHON      !calculate time evolution of the chem reactions
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON  
C calculates the time evolution of X and Y concentrations
C according to diffusion-reaction equations
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P7'
      INCLUDE 'IO.ALL'
C Local variables:
      REAL X(MAXX,MAXY)          !X species concentration
      REAL Y(MAXX,MAXY)          !Y species concentration
      REAL TIME                  !time
      REAL DT                    !time step
      REAL DTMIN,DTMAX           !limits on time step
      INTEGER IT                 !time index
      INTEGER IX,IY              !horiz and vert indices
      INTEGER NLINES             !number of lines printed to terminal
      INTEGER SCREEN             !send to terminal
      INTEGER PAPER              !make a hardcopy
      INTEGER FILE               !send to a file
      INTEGER NSTEP              !number of time steps to take
      LOGICAL MORE,NEWDT         !options for continuing
      INTEGER NFREQ              !graphing frequency (movies only)
      REAL XLO,XHI,XTOT          !min,max, and total conc of species X
      REAL YLO,YHI,YTOT          !min,max, and total conc of species Y
      REAL SGRF(MAXX,MAXY)       !species values for graphing
      INTEGER XSPEC,YSPEC        !flag to indicate species
C Functions:
      REAL GETFLT                !get floating point number from screen
      INTEGER YESNO              !get yes/no answer from screen
      LOGICAL LOGCVT             !change from 1,0 to true and false
      INTEGER GETINT             !get integer data from screen
      REAL RANNOS                !returns uniform random number
      DATA SCREEN,PAPER,FILE/1,2,3/
      DATA XSPEC,YSPEC/1,2/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     output summary of parameters
      IF (TTERM) CALL PRMOUT(OUNIT,NLINES)
      IF (TFILE) CALL PRMOUT(TUNIT,NLINES)
      IF (GFILE) CALL PRMOUT(GUNIT,NLINES)
C                 
C     initialize X and Y concentrations
      DO 5 IX=1,NX           !species concentrations fluctuate
         DO 6 IY=1,NY        !about equil values
            X(IX,IY)=A*(1.+XFLUCT*(RANNOS(DSEED)-.5))
            Y(IX,IY)=B/A*(1.+YFLUCT*(RANNOS(DSEED)-.5))
6        CONTINUE                                                  
5     CONTINUE
C
      TIME=0.             !initialize time            
      NSTEP=50            !default num of time steps until next prompt
      NFREQ=10            !default graphing freq
      DT=.07              !default, min, and max time steps
      DTMIN=1.E-05
      DTMAX=1.
C                                                                 
10    CONTINUE            !loop over different DT and/or NSTEP
        DT=GETFLT(DT,-DTMAX,DTMAX,'Enter time step')
        IF (ABS(DT) .LT. DTMIN) DT=SIGN(DTMIN,DT)    
        NSTEP=GETINT(NSTEP,1,1000,'Enter number of time steps')
        IF ((TTERM) .OR. (GTERM))
     +      NFREQ=GETINT(NFREQ,1,1000,
     +      'Enter display frequency for the concentration')
        NLINES=NLINES+6  
C           
C       calculate time step dependent constants
        IF (NX .GT. 1) THEN
          XXPLUS=-XDIFF*DT/HX/HX           
          YXPLUS=-YDIFF*DT/HX/HX
        ELSE 
          XXPLUS=0.
          YXPLUS=0.
        END IF
        IF (NY .GT. 1) THEN
          XYPLUS=-XDIFF*DT/HY/HY           
          YYPLUS=-YDIFF*DT/HY/HY
        ELSE
          XYPLUS=0.
          YYPLUS=0.
        END IF
        XXZERO=1-2*XXPLUS
        YXZERO=1-2*YXPLUS
        XYZERO=1-2*XYPLUS
        YYZERO=1-2*YYPLUS
        ADT=A*DT
        BDT=B*DT
        BP1DT=1.-(B+1.)*DT
        CALL TRDIAG              !calculate ALPHA and GAMMA for this DT
C        
        IF (TFILE) CALL TITLES(TUNIT,NLINES,DT)
C
15      CONTINUE                 !loop over sets of NSTEP time steps
         IF (TTERM) CALL TITLES(OUNIT,NLINES,DT)
C
          DO 20 IT=1,NSTEP       !time evolution
             TIME=TIME+DT
             CALL EVOLVE(X,Y,DT,XLO,XHI,XTOT,YLO,YHI,YTOT)
C
             IF (MOD(IT,NFREQ) .EQ. 0) THEN  !graphics output
               IF (GTERM) THEN
                 IF (TTERM) CALL PAUSE('to see concentrations...',1)
                 CALL GRFOUT(SCREEN,XSPEC,X,XLO,XHI,XTOT,SGRF,NX,NY)
                 CALL GRFOUT(SCREEN,YSPEC,Y,YLO,YHI,YTOT,SGRF,NX,NY)
                 IF ((TTERM) .AND. (IT .NE. NSTEP)) 
     +                    CALL TITLES(OUNIT,NLINES,DT)
               ELSE IF (TTERM) THEN
                 CALL PAUSE('to see X concentration...',1)
                 CALL DISPLY(X,OUNIT,XLO,XHI,XTOT,XSPEC)
                 CALL DISPLY(Y,OUNIT,YLO,YHI,YTOT,YSPEC)
                 CALL CLEAR
                 NLINES=0
                 IF (IT .NE. NSTEP) CALL TITLES(OUNIT,NLINES,DT)
               END IF
             ELSE                            !text output
               IF (TTERM) 
     +         CALL TXTOUT(OUNIT,NLINES,TIME,XLO,XHI,XTOT,YLO,YHI,YTOT)
               IF (TFILE) 
     +         CALL TXTOUT(TUNIT,NLINES,TIME,XLO,XHI,XTOT,YLO,YHI,YTOT)
             END IF
20        CONTINUE          
C 
          IF (TFILE) THEN                    !graphics output
             CALL DISPLY(X,TUNIT,XLO,XHI,XTOT,XSPEC)
             CALL DISPLY(Y,TUNIT,YLO,YHI,YTOT,YSPEC)
             CALL TITLES(TUNIT,NLINES,DT)
          END IF
          IF (GHRDCP) THEN
             CALL GRFOUT(PAPER,XSPEC,X,XLO,XHI,XTOT,SGRF,NX,NY)
             CALL GRFOUT(PAPER,YSPEC,Y,YLO,YHI,YTOT,SGRF,NX,NY)
          END IF
          IF (GFILE) CALL WRTOUT(X,Y,TIME,DT)
C
          MORE=LOGCVT(YESNO(1,'Continue iterating?'))
        IF (MORE) THEN
            NLINES=0
            IF (TTERM) CALL CLEAR
            GOTO 15
        END IF
C
        NEWDT=
     +  LOGCVT(YESNO(1,'Change time step or frequency and continue?'))
      IF (NEWDT) THEN
        NLINES=0
        IF (TTERM) CALL CLEAR
        GOTO 10
      END IF                                      
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE EVOLVE(X,Y,DT,XLO,XHI,XTOT,YLO,YHI,YTOT)
C evolves the species X and Y according to the diffusion-reaction
C equation and calculates minima, maxima, and average species values;
C all the reaction terms are in the horiz sweep of the lattice;
C see Eq 7.11-7.16 and VII.3a,b
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:              
      INCLUDE 'PARAM.P7'
C Input/Output variables:
      REAL X(MAXX,MAXY)      !X species concentration (I/O)
      REAL Y(MAXX,MAXY)      !Y species concentration (I/O)
      REAL DT                !time step (I)
      REAL XLO,XHI,XTOT      !min,max, and total conc of species X (out)
      REAL YLO,YHI,YTOT      !min,max, and total conc of species Y (out)
C Local variables:
      INTEGER IX,IY          !horiz and vert indices
                             !variables for matrix inversion
      REAL BETAXX(0:MAXX),BETAYX(0:MAXX),BETAXY(0:MAXY),BETAYY(0:MAXY) 
      REAL X2YDT,XRHS,YRHS   !temp variables
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (NX .GT. 1) THEN
        DO 10 IY=1,NY            !for each value of y, do an x sweep
           BETAXX(NX)=0.         !no flux bound. cond
           BETAYX(NX)=0.
           DO 20 IX=NX-1,0,-1    !backward recursion  
              X2YDT=Y(IX+1,IY)*DT*X(IX+1,IY)*X(IX+1,IY)
              XRHS=ADT+BP1DT*X(IX+1,IY)+X2YDT
              YRHS=-X2YDT+BDT*X(IX+1,IY)+Y(IX+1,IY)
              BETAXX(IX)=GAMXX(IX+1)*(XXPLUS*BETAXX(IX+1)-XRHS)
              BETAYX(IX)=GAMYX(IX+1)*(YXPLUS*BETAYX(IX+1)-YRHS)
20         CONTINUE
C
           X(1,IY)=BETAXX(0)/(1.-ALPHXX(0))    !no-flux b.c.
           Y(1,IY)=BETAYX(0)/(1.-ALPHYX(0))
           DO 30 IX=2,NX                       !forward recursion
              X(IX,IY)=ALPHXX(IX-1)*X(IX-1,IY)+BETAXX(IX-1)
              Y(IX,IY)=ALPHYX(IX-1)*Y(IX-1,IY)+BETAYX(IX-1)
30         CONTINUE
10      CONTINUE
      END IF
C
      IF (NY .GT. 1) THEN
        DO 40 IX=1,NX            !for each value of x, do a y sweep
           BETAXY(NY)=0.         !no flux bound. cond
           BETAYY(NY)=0.
           DO 50 IY=NY-1,0,-1    !backward recursion
              BETAXY(IY)=GAMXY(IY+1)*(XYPLUS*BETAXY(IY+1)-X(IX,IY+1))
              BETAYY(IY)=GAMYY(IY+1)*(YYPLUS*BETAYY(IY+1)-Y(IX,IY+1))
50         CONTINUE
C         
           X(IX,1)=BETAXY(0)/(1.-ALPHXY(0))      !no-flux b.c.
           Y(IX,1)=BETAYY(0)/(1.-ALPHYY(0))
           DO 60 IY=2,NY                         !forward recursion
              X(IX,IY)=ALPHXY(IY-1)*X(IX,IY-1)+BETAXY(IY-1)
              Y(IX,IY)=ALPHYY(IY-1)*Y(IX,IY-1)+BETAYY(IY-1)
60         CONTINUE
40      CONTINUE
      END IF
C
      XTOT=0.             !calculate min, max, and averages
      YTOT=0.
      XLO=X(1,1)
      YLO=Y(1,1)
      XHI=XLO
      YHI=YLO
      DO 70 IX=1,NX
         DO 80 IY=1,NY
            IF (X(IX,IY) .GT. XHI) XHI=X(IX,IY)
            IF (X(IX,IY) .LT. XLO) XLO=X(IX,IY)
            IF (Y(IX,IY) .GT. YHI) YHI=Y(IX,IY)
            IF (Y(IX,IY) .LT. YLO) YLO=Y(IX,IY)
            XTOT=XTOT+X(IX,IY)
            YTOT=YTOT+Y(IX,IY)
80       CONTINUE
70    CONTINUE
      XTOT=XTOT/NX/NY       !average value
      YTOT=YTOT/NX/NY
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TRDIAG
C calculate ALPHA and GAMMA for the inversion of the tridiagonal matrix
C see Eq. 7.11-7.16 and VII.3a,b
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:              
      INCLUDE 'PARAM.P7'
C Local variables:
      INTEGER IX,IY                    !lattice indices
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     horizontal recursion
      IF (NX .GT. 1) THEN
        ALPHXX(NX)=1.                 !no-flux boundary conditions
        ALPHYX(NX)=1.
        DO 10 IX=NX,1,-1              !backward recursion
           GAMXX(IX)=-1./(XXZERO+XXPLUS*ALPHXX(IX))
           ALPHXX(IX-1)=GAMXX(IX)*XXPLUS
           GAMYX(IX)=-1./(YXZERO+YXPLUS*ALPHYX(IX))
           ALPHYX(IX-1)=GAMYX(IX)*YXPLUS
10      CONTINUE                            
      END IF
C
C     vertical recursion
      IF (NY .GT. 1) THEN
        ALPHXY(NY)=1.                 !no-flux boundary conditions
        ALPHYY(NY)=1.
        DO 20 IY=NY,1,-1              !backward recursion
           GAMXY(IY)=-1./(XYZERO+XYPLUS*ALPHXY(IY))
           ALPHXY(IY-1)=GAMXY(IY)*XYPLUS
           GAMYY(IY)=-1./(YYZERO+YYPLUS*ALPHYY(IY))
           ALPHYY(IY-1)=GAMYY(IY)*YYPLUS
20      CONTINUE
      END IF
C
      RETURN
      END        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SUGGST
C suggests values for B based on linear stability analysis
C (see section VII.2), 
C also calculates graphing scales based on this analysis
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:                           
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'                                            
      INCLUDE 'PARAM.P7'
      INCLUDE 'MAP.P7'
C Local variables:
      REAL PI                 !3.14159
      REAL TEMP               !temp variable to find oscillations
      REAL DELTM,MPISQ        !functions of the mode M
      COMPLEX DISC            !discriminant
      REAL M,MOSC             !wave modes
      REAL BLO,BHI,BCRIT      !critical B values
      INTEGER IM              !mode index
      REAL AM,BM              !temp values to calculate frequencies
      REAL MAXFRQ             !maximum real part of the frequencies
      LOGICAL UNSTBL          !is the system probably unstable?
C Function:
      REAL GETFLT             !get floating point number from screen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C for 2-dimensions Meffective=sqrt(mx**2+my**2) and so needn't be
C an integer; if we replace M by Meffective in all the formulas
C we have the correct 2-dim analysis
C
      WRITE (OUNIT,5)
      WRITE (OUNIT,6)
      WRITE (OUNIT,*) ' '
5     FORMAT (' The following numbers from linear stability analysis',
     +        ' will guide you')
6     FORMAT (' in your choice of B for the chosen ',
     +        'values of A, XDIFF, and YDIFF:')
C
      MLO=0.  !lowest M value; =0 for no-flux bc; =SQRT(2) for fixed bc
C
C     what is the highest mode which gives rise to oscillations in time?
C     what values of B give oscillations? instability?
      PI=4.*ATAN(1.0)
      TEMP=1./(YDIFF-XDIFF)
C
      IF (TEMP .GT. 0.) THEN  !oscillations are possible
        MOSC=(SQRT(TEMP)/PI)                           
        WRITE (OUNIT,10) MOSC
10      FORMAT (' The highest mode to give oscillations is M = ',F7.3)
        DO 20 IM=1,2                     !loop over small and large M
          IF (IM .EQ. 1) M=MLO                !smallest M 
          IF (IM .EQ. 2) M=MOSC               !largest M  to give osc
          MPISQ=(M*PI)**2 
          DELTM=1+MPISQ*(XDIFF-YDIFF)
          BCRIT=1+A**2+MPISQ*(XDIFF+YDIFF)    !smallest B to give instb.
          DELTM=SQRT(DELTM)
          BLO=(A-DELTM)**2                    !smallest B to give osc
          BHI=(A+DELTM)**2                    !largest B to give osc
          WRITE (OUNIT,*) ' '
          WRITE (OUNIT,30) M,BCRIT
          WRITE (OUNIT,40) BLO,BHI
30        FORMAT (' for M = ',F7.3,
     +    ' the smallest B to give unstable oscillations is = ',F7.3)
40        FORMAT (' and oscillations arise only for B between ',F7.3,
     +    ' and ', F7.3)
20      CONTINUE
C
      ELSE IF (TEMP .LE. 0.) THEN  !oscillations are not possible
        MOSC=0.
        WRITE (OUNIT,50)                 
50      FORMAT (' There are no oscillations for these parameter values')
      END IF
C
C     what values of B and M give unstable, non-oscillatory behavior?
      M=SQRT(A/PI/PI/SQRT(XDIFF*YDIFF))
      MPISQ=M*M*PI*PI
      BCRIT=1+A**2*(XDIFF/YDIFF+1./YDIFF/MPISQ)+XDIFF*MPISQ
      WRITE (OUNIT,*) ' '
      WRITE (OUNIT,55)
      WRITE (OUNIT,60) M,BCRIT
55    FORMAT (' Instability w.r.t. exponentially growing behavior')
60    FORMAT (' first occurs for M = ',F7.3,' and B = ',F7.3)
C
C     allow user to adjust value of B using this information
      WRITE (OUNIT,*) ' '
      B=GETFLT(B,0.,20.,'Enter a value for B concentration')
      MREALS(IB)=B
C
C     calculate the frequencies for smallest M (depends on b.c.)
C     and largest M (which depends on the lattice spacing)
      MHI=SQRT(REAL(NX/2)**2+REAL(NY/2)**2)
      DO 70 IM=1,2
         IF (IM .EQ. 1) MPISQ=MLO*MLO*PI*PI
         IF (IM .EQ. 2) MPISQ=MHI*MHI*PI*PI
         AM=B-1-MPISQ*XDIFF
         BM=A**2+MPISQ*YDIFF
         DISC=(AM+BM)**2-4*A**2*B
         OMEGAP(IM)=(AM-BM+SQRT(DISC))/2
         OMEGAM(IM)=(AM-BM-SQRT(DISC))/2
70    CONTINUE
C
      UNSTBL=.FALSE.                !is the system unstable?
      MAXFRQ=REAL(OMEGAP(1))        !find the max freq
      IF (REAL(OMEGAP(2)) .GT. MAXFRQ) MAXFRQ=REAL(OMEGAP(2))
      IF (REAL(OMEGAM(1)) .GT. MAXFRQ) MAXFRQ=REAL(OMEGAM(1))
      IF (REAL(OMEGAM(2)) .GT. MAXFRQ) MAXFRQ=REAL(OMEGAM(2))
      IF (MAXFRQ .GE. 0) UNSTBL=.TRUE.
      IF (MAXFRQ .LT. .1) MAXFRQ=2*MAXFRQ  !adjust very small scales
      IF (UNSTBL) THEN              !scale for graphing
          YSCALE=MAXFRQ*B/A         !this scale is set for all times
          XSCALE=MAXFRQ*A           !for these values of the parameters
      ELSE
          YSCALE=2*YFLUCT
          XSCALE=2*XFLUCT
      END IF
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
      INCLUDE 'PARAM.P7'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)         
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get environment parameters
      CALL SETUP                
C 
C     display header screen     
      DESCRP(1)= 'PROJECT 7'
      DESCRP(2)= 'The Brusselator in two dimensions'
      NHEAD=2
C 
C     text output description
      DESCRP(3)= 'time, X and Y species minima, maxima and averages'
      NTEXT=1
C 
C     graphics output description
      DESCRP(4)= 'X and Y species concentrations'
      NGRAPH=1
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C 
C     setup menu arrays, beginning with constant part
      CALL MENU                      
C                
      MTYPE(13)=FLOAT
      MPRMPT(13)='Enter value for A species concentration'
      MTAG(13)='A concentration'
      MLOLIM(13)=0.001
      MHILIM(13)=20.
      MREALS(13)=2.
C
      MTYPE(14)=FLOAT
      MPRMPT(14)='Enter value for B species concentration'
      MTAG(14)='B concentration'
      MLOLIM(14)=0.
      MHILIM(14)=20.
      MREALS(14)=4.
C
      MTYPE(15)=FLOAT
      MPRMPT(15)='Enter value for XDIFF (X species diffusion constant)'
      MTAG(15)='X species diffusion constant'
      MLOLIM(15)=0.
      MHILIM(15)=10.
      MREALS(15)=0.001
C
      MTYPE(16)=FLOAT
      MPRMPT(16)='Enter value for YDIFF (Y species diffusion constant)'
      MTAG(16)='Y species diffusion constant'
      MLOLIM(16)=0.
      MHILIM(16)=10.
      MREALS(16)=0.003
C
      MTYPE(17)=FLOAT
      MPRMPT(17)='Enter value for XFLUCT (X species fluctuations)'
      MTAG(17)='X species fluctuation'
      MLOLIM(17)=-1.
      MHILIM(17)=1.
      MREALS(17)=.01
C
      MTYPE(18)=FLOAT
      MPRMPT(18)='Enter value for YFLUCT (Y species fluctuations)'
      MTAG(18)='Y species fluctuation'
      MLOLIM(18)=-1.
      MHILIM(18)=1.
      MREALS(18)=.01
C
      MTYPE(19)=SKIP
      MREALS(19)=35.
C 
      MTYPE(38)=NUM
      MPRMPT(38)= 'Enter number of X lattice points'
      MTAG(38)= 'Number of X lattice points'
      MLOLIM(38)=1.
      MHILIM(38)=MAXX
      MINTS(38)=20
C 
      MTYPE(39)=NUM
      MPRMPT(39)= 'Enter number of Y lattice points'
      MTAG(39)= 'Number of Y lattice points'
      MLOLIM(39)=1.
      MHILIM(39)=MAXY
      MINTS(39)=20
C 
      MTYPE(40)=NUM
      MPRMPT(40)= 'Integer random number seed for init fluctuations'
      MTAG(40)= 'Random number seed'
      MLOLIM(40)=1000.
      MHILIM(40)=99999.
      MINTS(40)=54765
C 
      MTYPE(41)=SKIP
      MREALS(41)=60.
C 
      MSTRNG(MINTS(75))= 'proj7.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C 
      MSTRNG(MINTS(86))= 'proj7.grf'
C                     
      MTYPE(87)=NUM
      MPRMPT(87)= 'Enter number of contour levels'
      MTAG(87)= 'Number of contour levels'
      MLOLIM(87)=2.
      MHILIM(87)=50.
      MINTS(87)=6
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
      INCLUDE 'PARAM.P7'
      INCLUDE 'MAP.P7'
C Functions:
      LOGICAL LOGCVT          !converts 1 and 0 to true and false 
      REAL GETFLT             !get floating point number from screen
      INTEGER GETINT          !get integer from screen
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
      A=MREALS(IA)
      B=MREALS(IB)
      XDIFF=MREALS(IXDIFF)
      YDIFF=MREALS(IYDIFF)
      XFLUCT=MREALS(IXFL)
      YFLUCT=MREALS(IYFL)
      NX=MINTS(INX)
      NY=MINTS(INY)
      DSEED=DBLE(MINTS(IDSEED))
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
      NLEV=MINTS(INLEV)
C 
C     open files 
      IF (TFILE) CALL FLOPEN(TNAME,TUNIT)
      IF (GFILE) CALL FLOPEN(GNAME,GUNIT)
      !files may have been renamed
      MSTRNG(MINTS(ITNAME))=TNAME
      MSTRNG(MINTS(IGNAME))=GNAME
C
C     make sure that problem is at least one-dimensional
60    IF ((NX .EQ. 1) .AND. (NY .EQ. 1)) THEN
         WRITE (OUNIT,50)
50       FORMAT (' Both NX and NY can''t be = 1')
         NX=GETINT(20,1,MAXX,'Enter new value for NX')
         NY=GETINT(20,1,MAXY,'Enter new value for NY')
         MINTS(INX)=NX
         MINTS(INY)=NY
         GOTO 60
      END IF
C
      CALL SUGGST                   !suggest values for B
      CALL CLEAR
C
C     derivative parameters
      HX=0.                         !allow for one-dimensional systems
      HY=0.                         
      IF (NX .GT. 1) HX=1./(NX-1)   !calculate space steps 
      IF (NY .GT. 1) HY=1./(NY-1)
C     calculate parameters for best looking text display
      IF (2*NX .LE. TRMWID) THEN
           XSKIP=.TRUE.              !skip spaces in x
           XCNTR=(TRMWID-2*NX)/2     !how to center display
      ELSE
           XSKIP=.FALSE.
           XCNTR=(TRMWID-NX)/2
      END IF
      IF (XCNTR .LT. 1) XCNTR=1
C
      IF (2*NY .LE. TRMLIN-4) THEN
          YSKIP=.TRUE.               !skip lines in y
          YCNTR=(TRMLIN-2*NY)/2-2    !how to center display
      ELSE
          YSKIP=.FALSE.
          YCNTR=(TRMLIN-NY)/2-2
      END IF
      IF (YCNTR .LT. 0) YCNTR=0
C
      RETURN               
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DISPLY(S,MUNIT,SLO,SHI,STOT,SPEC)
C display species (S) as ascii characters
C values above the equilibrium value are displayed as capital letters,
C those less than the equilibrium value are displayed as small letters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:          
      INCLUDE 'PARAM.P7'
      INCLUDE 'IO.ALL'
C Input variables:
      REAL S(MAXX,MAXY)          !species values
      INTEGER MUNIT              !unit we're writing to
      REAL SLO,SHI,STOT          !min,max, and total conc of species S
      INTEGER SPEC               !which species
C Local variables:
      INTEGER IX,IY              !lattice indices
      INTEGER STEMP              !S as an integer
      CHARACTER*1 SPECS(MAXX)    !species as character data
      REAL SZERO                 !'zero' value for species
      REAL SCALE                 !scale for species
      CHARACTER*80 BLNK          !blanks for centering in X
      INTEGER XSPEC,YSPEC        !flag to indicate species
      INTEGER CSCALE             !scale of ascii data
      CHARACTER*1 ASKII(0:25),NEGASK(0:25)!charac data for display
      DATA BLNK /' '/
      DATA ASKII/'A','B','C','D','E','F','G','H','I','J','K','L','M',
     +           'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      DATA NEGASK/'a','b','c','d','e','f','g','h','i','j','k','l','m',
     +           'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA XSPEC,YSPEC/1,2/
      DATA CSCALE/25/                     !since there are 26 charac
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) THEN
        CALL CLEAR
        DO 20 IY=1,YCNTR                   !center output
           WRITE (OUNIT,*) ' '
20      CONTINUE                           !print out data
        IF (SPEC .EQ. XSPEC) WRITE (OUNIT,25) SLO,SHI,STOT
        IF (SPEC .EQ. YSPEC) WRITE (OUNIT,26) SLO,SHI,STOT
25      FORMAT(15X,'X min =',F7.3,'  X max =',F7.3,'  X average =',F7.3)
26      FORMAT(15X,'Y min =',F7.3,'  Y max =',F7.3,'  Y average =',F7.3)
      ELSE
        WRITE (MUNIT,*) ' '
        IF (SPEC .EQ. XSPEC) WRITE(MUNIT,*)
     +               ' X species concentration - A:'
        IF (SPEC .EQ. YSPEC) WRITE(MUNIT,*)
     +               ' Y species concentration - A/B:'
      END IF
C
      IF (SPEC .EQ. XSPEC) THEN            !set scales
            SZERO=A                        !X equil value
            SCALE=XSCALE/CSCALE
      ELSE IF (SPEC .EQ. YSPEC) THEN
            SZERO=B/A                      !Y equil value
            SCALE=YSCALE/CSCALE
      END IF
C
      DO 100 IY=NY,1,-1                    
         DO 50 IX=1,NX       !STEMP is scaled deviation from SZERO
            STEMP=NINT((S(IX,IY)-SZERO)/SCALE)   
            IF (STEMP .GT. CSCALE)   STEMP=CSCALE
            IF (STEMP .LT. -CSCALE)  STEMP=-CSCALE
            IF (STEMP .GT. 0) THEN         !convert species to ASCII
                SPECS(IX)=ASKII(STEMP) 
            ELSE IF (STEMP .LT. 0) THEN
                SPECS(IX)=NEGASK(-STEMP) 
            ELSE IF (STEMP .EQ. 0) THEN
                IF (S(IX,IY) .GT. SZERO) THEN
                   SPECS(IX)=ASKII(0)
                ELSE
                   SPECS(IX)=ASKII(0)
                END IF
            END IF
50       CONTINUE    
C
C        write out a line at a time (no centering done for TUNIT)
         IF (MUNIT .EQ. TUNIT) THEN
                  WRITE (TUNIT,16) (SPECS(IX),IX=1,NX)
         ELSE
           IF (XSKIP) THEN
             WRITE (OUNIT,10) BLNK(1:XCNTR),(SPECS(IX),IX=1,NX)
           ELSE                   
             WRITE (OUNIT,15) BLNK(1:XCNTR),(SPECS(IX),IX=1,NX)
           END IF
           IF (YSKIP) WRITE (OUNIT,*) ' '
         END IF
10       FORMAT (1X,A,100(A1,1X))                    
15       FORMAT (1X,A,100(A1))
16       FORMAT (1X,100(A1))
100   CONTINUE
      IF (MUNIT .EQ. OUNIT) THEN
          IF (SPEC .EQ. XSPEC) THEN
              CALL PAUSE('to see Y concentration...',0)
          ELSE
              CALL PAUSE('to continue...',0)
          END IF
      END IF
C
      RETURN
      END  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE PRMOUT(MUNIT,NLINES)
C outputs parameter summary to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:                                                     
       INCLUDE 'IO.ALL'
       INCLUDE 'PARAM.P7'
C Input/Output variables:               
       INTEGER MUNIT            !unit number for output (input)
       INTEGER NLINES           !number of lines written so far (output)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (MUNIT .EQ. OUNIT) CALL CLEAR
C 
       WRITE (MUNIT,2)
       WRITE (MUNIT,4)
       WRITE (MUNIT,6) A,B
       WRITE (MUNIT,8) XDIFF,YDIFF
       WRITE (MUNIT,10) XFLUCT,YFLUCT
       WRITE (MUNIT,12) NX,NY
       WRITE (MUNIT,16) MLO,OMEGAP(1),OMEGAM(1)
       WRITE (MUNIT,16) MHI,OMEGAP(2),OMEGAM(2)
       WRITE (MUNIT,14)
       WRITE (MUNIT,2)
C
       NLINES=10
C                                 
2      FORMAT (' ')
4      FORMAT (' Project 7: Brusselator in two dimensions')
6      FORMAT (' A concentration = ',F8.3,5X,'B concentration =',F8.3)
8      FORMAT (' X species diffusion = ',1PE12.5,5X,
     +         '  Y species diffusion = ',1PE12.5)
10     FORMAT (' X species init fluct = '1PE12.5,5X,
     +         ' Y species init fluct = '1PE12.5)
12     FORMAT (' NX = ',I4,5X,' NY = ',I4)
14     FORMAT (' all values are in scaled units')
16     FORMAT (' for M=',F7.3,' the frequencies are ', 
     +         2( 2X,'(',F7.3,',',F7.3,')' ) )
C
       RETURN                                             
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,NLINES,TIME,XLO,XHI,XTOT,YLO,YHI,YTOT)
C writes results for one time step to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
C Input variables:
      INTEGER MUNIT              !output unit specifier
      REAL TIME                  !time
      INTEGER NLINES             !number of lines printed so far (I/O)
      REAL XLO,XHI,XTOT          !min,max, and total conc of species X
      REAL YLO,YHI,YTOT          !min,max, and total conc of species Y
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     if screen is full, clear screen and retype headings          
      IF ((MOD(NLINES,TRMLIN-4) .EQ. 0) 
     +    .AND. (MUNIT .EQ. OUNIT) .AND. (.NOT. GTERM)) THEN
         CALL PAUSE('to continue...',1)
         CALL CLEAR
         WRITE (MUNIT,14)
         WRITE (MUNIT,16)
         NLINES=NLINES+2
      END IF
C
      WRITE (MUNIT,40) TIME,XLO,XHI,XTOT,YLO,YHI,YTOT
C     keep track of printed lines only for terminal output
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+1
C
40    FORMAT (3X,1PE12.5,5X,6(0PF8.3,2X))
14    FORMAT (8X,'time',11X,'Xlo',7X,'Xhi',7X,'Xave',7X,'Ylo',
     +         7X,'Yhi',6X,'Yave')
16    FORMAT (8X,'----',11X,'---',7X,'---',7X,'----',7X,'---',
     +         7X,'---',6X,'----')
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE TITLES(MUNIT,NLINES,DT)
C write out time step and column titles to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'IO.ALL'
C Passed variables:               
       INTEGER MUNIT            !unit number for output (input)
       INTEGER NLINES           !number of lines written so far (output)
       REAL DT                  !time step (input)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF ((MOD(NLINES,TRMLIN-7) .EQ. 0) 
     +    .AND. (MUNIT .EQ. OUNIT) .AND. (.NOT. GTERM)) 
     +    CALL CLEAR
C
       WRITE (MUNIT,*)  ' '
       WRITE (MUNIT,20) DT
       WRITE (MUNIT,14)
       WRITE (MUNIT,16)
C 
       IF (MUNIT .EQ. OUNIT) NLINES=NLINES+4
C
20     FORMAT (' time step = ',1PE12.5)
14     FORMAT (8X,'time',11X,'Xlo',7X,'Xhi',7X,'Xave',7X,'Ylo',
     +         7X,'Yhi',6X,'Yave')
16     FORMAT (8X,'----',11X,'---',7X,'---',7X,'----',7X,'---',
     +         7X,'---',6X,'----')
C
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFOUT(DEVICE,SPEC,S,SLO,SHI,STOT,SGRF,MX,MY)
C display contours of the species-equilibrium value
C this routine will also do 1-dim plots if either NX or NY =1
C     the field values must be in an array that is exactly NX by NY;
C     this can be accomplished with implicit dimensioning which 
C     requires that SGRF and its dimensions be passed to this routine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'   
      INCLUDE 'PARAM.P7'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE                !which device
      INTEGER SPEC                  !which species
      REAL S(MAXX,MAXY)             !species values
      INTEGER MX,MY                 !NX and NY in disguise
      REAL SGRF(MX,MY)              !species values for graphing
      REAL SLO,SHI,STOT             !min,max, average species values
C Local variables:
      REAL SZERO                    !'zero' value of species
      INTEGER IX,IY                 !lattice indices
      INTEGER SCREEN                !send to terminal
      INTEGER PAPER                 !make a hardcopy
      INTEGER FILE                  !send to a file
      INTEGER NCONT                 !number of contours
      REAL LATTIC(1:MAXX)           !array for 1-dim plots
      CHARACTER*9 CLO,CHI,CTOT      !data as characters
      CHARACTER*9 CA,CB             !data as characters
      INTEGER XSPEC,YSPEC           !flag to indicate species
      INTEGER LLO,LHI,LTOT,LA,LB    !length of char strings
      REAL SMIN,SMAX        
      DATA XSPEC,YSPEC/1,2/
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (SPEC .EQ. XSPEC) SZERO=A       !set equil values
      IF (SPEC .EQ. YSPEC) SZERO=B/A
      DO 10 IX=1,MX             
         DO 11 IY=1,MY
           SGRF(IX,IY)=S(IX,IY)-SZERO    !load field into SGRF
11       CONTINUE                        !display deviations from Szero
10    CONTINUE
      IF (MX .EQ. 1) THEN
          DO 15 IY=1,MY
             LATTIC(IY)=REAL(IY)         !need to reload data for  
             SGRF(IY,1)=S(1,IY)-SZERO    !one-dim plots
15        CONTINUE                          
      END IF
      IF (MY .EQ. 1) THEN
          DO 16 IX=1,MX
             LATTIC(IX)=REAL(IX)
16        CONTINUE
      END IF
C
C     messages for the impatient
      IF ((DEVICE .NE. SCREEN) .AND. (SPEC .EQ. XSPEC)) WRITE(OUNIT,100) 
C
C     calculate parameters for graphing
      NPLOT=2                       !how many plots
C
      YMAX=MY                       !axis data
      YMIN=1.
      XMIN=1.
      XMAX=MX
      Y0VAL=XMIN
      X0VAL=YMIN
      NXTICK=5
      NYTICK=5
      NPOINT=5
C
      LABEL(1)='NX'                 !graph description
      LABEL(2)='NY'
      CALL CONVRT(SLO,CLO,LLO)
      CALL CONVRT(SHI,CHI,LHI)
      CALL CONVRT(STOT,CTOT,LTOT)
      CALL CONVRT(A,CA,LA)
      CALL CONVRT(B,CB,LB)
      TITLE='Brusselator with A='//CA(1:LA)//' and B='//CB(1:LB)
C
      IF (MX .EQ. 1) THEN           !labels and axes are different for
         XMIN=1.                    !one-dim plots
         XMAX=MY
         Y0VAL=XMIN
         LABEL(1)='NY'
         IF (SPEC .EQ. XSPEC) THEN
             LABEL(2)='X Concentration - A'
             YMIN=MAX(-XSCALE,-A)                 !scale is the same
             YMAX=XSCALE                          !for all times
         ELSE IF (SPEC .EQ. YSPEC) THEN
             LABEL(2)='Y Concentration - B/A'
             YMIN=MAX(-YSCALE,-B/A)
             YMAX=YSCALE
         END IF
         X0VAL=YMIN
         NPOINT=NY
      ELSE IF (MY .EQ. 1) THEN
         IF (SPEC .EQ. XSPEC) THEN
             LABEL(2)='X Concentration - A'
             YMIN=MAX(-XSCALE,-A)
             YMAX=XSCALE
         ELSE IF (SPEC .EQ. YSPEC) THEN
             LABEL(2)='Y Concentration - B/A'
             YMIN=MAX(-YSCALE,-B/A)
             YMAX=YSCALE
         END IF
         X0VAL=YMIN
         NPOINT=NX
      END IF
C
      IF (SPEC .EQ. XSPEC) THEN   !species-dependent parameters
         IPLOT=1
         SMIN=MAX(-XSCALE,-A)     !scale for species
         SMAX=XSCALE
         INFO='Xmin='//CLO(1:LLO)//' Xmax='
     +             //CHI(1:LHI)//' Xave='//CTOT(1:LTOT)
         CALL GTDEV(DEVICE)                   !device nomination
         IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
      ELSE IF (SPEC .EQ. YSPEC) THEN
         IPLOT=2
         SMIN=MAX(-YSCALE,-B/A)
         SMAX=YSCALE
         INFO='Ymin='//CLO(1:LLO)//' Ymax='
     +        //CHI(1:LHI)//' Yave='//CTOT(1:LTOT)
      END IF
C 
      CALL LNLNAX                          !draw axes
C
      IF ((MX .EQ. 1) .OR. (MY .EQ. 1)) THEN
          CALL XYPLOT(LATTIC,SGRF)
      ELSE
          CALL CONTOR(SGRF,MX,MY,SMIN,SMAX,NLEV)
      END IF
C 
C     end graphing session
      IF (SPEC .EQ. YSPEC) THEN
       IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)  !close graphics package
       IF (DEVICE .EQ. SCREEN) CALL TMODE        !switch to text mode
      END IF
C 
100   FORMAT (/,' Patience, please; output going to a file.')
C 
      RETURN
      END                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE WRTOUT(X,Y,TIME,DT)
C write out X and Y species concentrations to GUNIT
C for graphing with an external package
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P7'
C Input variables:
      REAL TIME                  !time
      REAL DT                    !time step
      REAL X(MAXX,MAXY)          !X species concentration
      REAL Y(MAXX,MAXY)          !Y species concentration
C Local variables:
      INTEGER IX,IY              !lattice indices
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE (GUNIT,20) TIME,DT
20    FORMAT (' time = ',1PE12.5,'  time step=',1PE12.5)
      WRITE (GUNIT,*) ' X concentration'
      DO 100 IX=1,NX                            
        WRITE (GUNIT,10) (X(IX,IY),IY=1,NY)
100   CONTINUE
      WRITE (GUNIT,*) ' Y concentration'
      DO 200 IX=1,NX                            
        WRITE (GUNIT,10) (Y(IX,IY),IY=1,NY)
200   CONTINUE
10    FORMAT (5(2X,1PE14.7))
C
      RETURN
      END                                                               
