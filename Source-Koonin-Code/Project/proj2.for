CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM PROJ2
C     Project 2: The structure of white dwarf stars
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop; execute once for each set of param
         CALL PARAM       !get input parameters
         CALL ARCHON      !calculate density, mass and radius of stars
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON
C calculate the mass(r), density(r) for a range of central densities
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P2'
C Local variables:
      REAL MSTOR(MAXMOD)            !mass of star
      REAL RADIUS(MAXMOD)           !radius of star
      REAL STEP(MAXMOD)             !step size
      REAL RMIN(MAXMOD)             !starting radius
      INTEGER NSTEP(MAXMOD)         !number of radial steps
      REAL RHOSTR(MAXMOD,0:MAXSTP)  !density(r) of star
      INTEGER IMODEL                !index of current model
      INTEGER NLINES                !number of lines printed out
      REAL M,RHO                    !current mass and density
      REAL R                        !current radius
      REAL DR                       !radial step
      REAL EGRAV                    !gravitational energy of star
      REAL EKINT                    !kinetic energy of star
      REAL EXPONT                   !temporary exponent
      REAL RHOCEN                   !central density for this model
      INTEGER IR                    !number of radial steps
      INTEGER DEVICE                !current graphing device
      INTEGER SCREEN                !send to terminal
      INTEGER PAPER                 !make a hardcopy
      INTEGER FILE                  !send to a file
C Functions:
      REAL GAMMA                    !dPressure/dDensity
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     output summary of parameters
      IF (TTERM) CALL PRMOUT(OUNIT,NLINES)
      IF (TFILE) CALL PRMOUT(TUNIT,NLINES)
      IF (GFILE) CALL PRMOUT(GUNIT,NLINES)
C
      DO 20 IMODEL=1,NMODEL                      !loop over central dens
C        
C       calculate central density
        EXPONT=FLOAT(IMODEL-1)/FLOAT(NMODEL-1)   !spacing between dens
        RHOCEN=(RHO1/RHO0)*((RHO2/RHO1)**EXPONT) !equal on log scale
C
        DR=((3.0*0.001/RHOCEN)**(1.0/3.0))/3.0   !radial step
        R=DR/10.                                 !begin at finite radius
        RMIN(IMODEL)=R*R0/10.E+08                !starting radius
        STEP(IMODEL)=DR*R0/10.E+08               !step size
C
C       use Taylor expansion to obtain initial conditions
        M=(R**3.0)*(RHOCEN/3.0)                  !initial mass
C       initial density
        RHO=RHOCEN*(1-RHOCEN*R**2.0/6.0/GAMMA(RHOCEN**(1./3.))) 
        RHOSTR(IMODEL,0)=RHO*RHO0
C
C       integrate equations and find energies, mass, and radius of star
        CALL INTGRT(M,RHO,R,DR,EGRAV,EKINT,IMODEL,IR,RHOSTR)
C
        MSTOR(IMODEL)=M*M0/10.E+33               !save values for graph
        RADIUS(IMODEL)=R*R0/10.E+08
        NSTEP(IMODEL)=IR                          
C
        IF (TTERM)                               !text output
     +      CALL TXTOUT(OUNIT,NLINES,R,M,DR,RHOCEN,EKINT,EGRAV,IR)
        IF (TFILE) 
     +      CALL TXTOUT(TUNIT,NLINES,R,M,DR,RHOCEN,EKINT,EGRAV,IR)
C                                 
20    CONTINUE
C
      IF (TTERM) CALL PAUSE('to continue...',1)   
      IF (TTERM) CALL CLEAR
C
C     graphics output
      IF (GTERM) CALL GRFOUT(SCREEN,MSTOR,RADIUS,STEP,RMIN,NSTEP,RHOSTR)
      IF (GFILE) CALL GRFOUT(FILE,MSTOR,RADIUS,STEP,RMIN,NSTEP,RHOSTR)
      IF (GHRDCP) CALL GRFOUT(PAPER,MSTOR,RADIUS,STEP,RMIN,NSTEP,RHOSTR)
C      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTGRT(M,RHO,R,DR,EGRAV,EKINT,IMODEL,IR,RHOSTR)
C integrates mass and density beginning with M and RHO at radius R;
C calculates gravitational and kinetic energies
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P2'
C Passed variables:
      REAL RHOSTR(MAXMOD,0:MAXSTP)!density(r) of star (output)
      INTEGER IMODEL              !index of current model (input)
      REAL M,RHO                  !current mass and density (I/O)
      REAL R                      !current radius (I/O)
      REAL DR                     !radial step (input)
      REAL EGRAV                  !gravitational energy of star (output)
      REAL EKINT                  !kinetic energy of star (output)
      INTEGER IR                  !current lattice point (output)
C Function:
      REAL EPSRHO                 !function in energy density
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      EGRAV=0.0                   !zero sums
      EKINT=0.0
      IR=0
10    IF ((IR .LT. MAXSTP) .AND. (RHO .GT. RHOCUT)) THEN
         IR=IR+1                          !update radial index
         CALL RUNGE(M,RHO,R,DR)           !take a Runge-Kutta step
         IF (RHO .LT. RHOCUT) RHO=RHOCUT  !avoid small densities
         RHOSTR(IMODEL,IR)=RHO*RHO0       !save values for graphing
         EGRAV=EGRAV+M*RHO*R              !contribution to energy integ.
         EKINT=EKINT+R**2*EPSRHO(RHO**(1./3.))
      GOTO 10
      END IF
      EGRAV=-EGRAV*DR*E0                  !unscaled energies
      EKINT=EKINT*DR*E0
      RETURN     
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RUNGE(M,RHO,R,DR)
C take a Runge-Kutta step to integrate mass and density
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
      REAL M,RHO                    !current mass and density (I/O)
      REAL R                        !current radius (input)
      REAL DR                       !radial step (input)
C Local variables:
      REAL DMDR,DRHODR              !radial derivatives of mass and dens
      REAL K1M,K1RHO                !intermediate increments of mass and
      REAL K2M,K2RHO                !                           density
      REAL K3M,K3RHO
      REAL K4M,K4RHO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL DERS(M,RHO,R,DMDR,DRHODR) !calculate k1's
      K1M=DR*DMDR
      K1RHO=DR*DRHODR
C
      R=R+DR/2.0                     !calculate k2's
      M=M+.5*K1M
      RHO=RHO+.5*K1RHO
      CALL DERS(M,RHO,R,DMDR,DRHODR)
      K2M=DR*DMDR
      K2RHO=DR*DRHODR
C
      M=M+.5*(K2M-K1M)               !calculate k3's
      RHO=RHO+.5*(K2RHO-K1RHO)
      CALL DERS(M,RHO,R,DMDR,DRHODR)
      K3M=DR*DMDR
      K3RHO=DR*DRHODR
C
      R=R+DR/2.0                     !calculate k4's
      M=M+K3M-.5*K2M
      RHO=RHO+K3RHO-.5*K2RHO
      CALL DERS(M,RHO,R,DMDR,DRHODR)
      K4M=DR*DMDR
      K4RHO=DR*DRHODR
C
C     values of mass and density at new radius
      M=M-K3M+(K1M+2.0*K2M+2.0*K3M+K4M)/6.0
      RHO=RHO-K3RHO+(K1RHO+2.0*K2RHO+2.0*K3RHO+K4RHO)/6.
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DERS(M,RHO,R,DMDR,DRHODR)
C calculate the derivatives of mass and density with respect to radius
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
      REAL M,RHO           !current mass and density (input)
      REAL R               !current radius (input)
      REAL DMDR,DRHODR     !radial derivatives of mass and dens (output)
C Functions
      REAL GAMMA           !dPressure/dDensity
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (RHO .GT. 0.0) THEN   !avoid division by zero
          DRHODR=-M*RHO/((R**2)*GAMMA(RHO**(1.0/3.0)))
          DMDR=(R**2)*RHO
      END IF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION GAMMA(X)
C derivative of pressure with respect to density
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
      REAL X                           !cubed root of density
      GAMMA=X*X/3.0/SQRT(X*X+1)
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION EPSRHO(X)
C function needed in the calculation of the energy
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
      REAL X                          !cubed root of density
C Local variables:
      REAL PART1,PART2                !terms in the function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PART1=X*(1.0+2.0*X**2)*SQRT(1.0+X**2)
      PART2=LOG(X+SQRT(1.0+X**2))
      EPSRHO=0.375*(PART1-PART2)
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INIT
C initializes constants, displays header screen,
C initializes arrays for input parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'MENU.ALL'
      INCLUDE 'PARAM.P2'
C Local parameters:
      CHARACTER*80 DESCRP        !program description
      DIMENSION DESCRP(22)
      INTEGER NHEAD,NTEXT,NGRF   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL SETUP                !get environment parameters
C 
C     display header screen     
      DESCRP(1)= 'PROJECT 2'
      DESCRP(2)= 'The structure of white dwarf stars'
      NHEAD=2
C
C     text output description
      DESCRP(3)= 'mass, radius, and energy of each star'
      NTEXT=1
C
C     graphics output description
      DESCRP(4)= 'total mass vs. final radius for all models'
      DESCRP(5)= 'density vs. radius for each model'
      NGRF=2
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRF)
C
C     calculate constants                
      MEORMP = 1./1836.
C
C     setup menu
      CALL MENU                       !setup constant part of menu
C
      MPRMPT(4)='2) (not used)'
C                
      MTYPE(13)=FLOAT
      MPRMPT(13)= 'Enter electron fraction, YE'
      MTAG(13)= 'Electron fraction, YE'
      MLOLIM(13)=.001
      MHILIM(13)=10.
      MREALS(13)=1.0
C                
      MTYPE(14)=FLOAT
      MPRMPT(14)= 'Enter central density (gm/cm^3) for the first model'
      MTAG(14)= 'Initial central density (gm/cm^3)'
      MLOLIM(14)=1.E5
      MHILIM(14)=1.E15
      MREALS(14)=1.E5
C 
      MTYPE(15)=FLOAT
      MPRMPT(15)= 'Enter central density (gm/cm^3) for the last model'
      MTAG(15)= 'Final central density (gm/cm^3)'
      MLOLIM(15)=1.E5
      MHILIM(15)=1.E15
      MREALS(15)=1.E11
C               
      MTYPE(16)=NUM
      MPRMPT(16)= 'Enter the number of models to calculate'
      MTAG(16)= 'Number of models to calculate'
      MLOLIM(16)=1.
      MHILIM(16)=MAXMOD
      MINTS(16)=4.
C                                                 
      MTYPE(17)=SKIP
      MREALS(17)=35.
C
      MTYPE(36)=SKIP
      MREALS(36)=60.
C
      MSTRNG(MINTS(75))= 'proj2.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C            
      MSTRNG(MINTS(86))= 'proj2.grf'
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
      INCLUDE 'PARAM.P2'
      INCLUDE 'MAP.P2'
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
C     set new parameter values 
C     physical and numerical:
      YE=MREALS(IYE)
      RHO1=MREALS(IRHO1)
      RHO2=MREALS(IRHO2)
      NMODEL=MINTS(INMODL)
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
      CALL PCHECK
      CALL CLEAR
C
C     calculate parameter dependent quantities
      R0=7.72E+08*YE             !scaling radius
      M0=5.67E+33*(YE**2)        !scaling mass
      RHO0=979000./YE            !scaling density
      E0=YE*MEORMP*9*(M0/1E+31)  !scaling energy
      RHOCUT=1000./RHO0          !nearly zero density
C      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE PCHECK
C ensure that RHO1 < RHO2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global parameters:
      INCLUDE 'PARAM.P2'
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'
      INCLUDE 'MAP.P2'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
10    IF (RHO1 .GE. RHO2) THEN
          WRITE (OUNIT,15) 
          CALL ASK(IRHO1,IRHO2)
          RHO1=MREALS(IRHO1)
          RHO2=MREALS(IRHO2)
      GOTO 10
      END IF
15    FORMAT 
     + (' The initial density must be smaller than the final density')
      RETURN
      END 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRMOUT(MUNIT,NLINES)
C outputs parameters to the specified unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables
      INCLUDE 'PARAM.P2'
      INCLUDE 'IO.ALL'
C Passed variables
      INTEGER MUNIT              !unit number for output (input)
      INTEGER NLINES             !number of lines written so far (I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) CALL CLEAR
C
      WRITE (MUNIT,19)
      WRITE (MUNIT,21)
      WRITE (MUNIT,25)  NMODEL
      WRITE (MUNIT,27)  YE
      WRITE (MUNIT,29)  RHO1,RHO2
      WRITE (MUNIT,19)
C
      IF (MUNIT .NE. GUNIT) THEN       !headers for text output only
        WRITE (MUNIT,31)
        WRITE (MUNIT,32)
        WRITE (MUNIT,33)
      END IF
C
      NLINES = 9
C 
19    FORMAT  (' ')
21    FORMAT  (' Output from project 2: ',
     +         'The structure of white dwarf stars ')
25    FORMAT  (' Number of models = ',I3)
27    FORMAT  (' Electron fraction YE= ',1PE9.3)
29    FORMAT  (' First central density=',1PE9.3, 
     +         5X,' Final central density=',1PE9.3)
31    FORMAT  (3X,'Step',4X,'Cent Den',1X,'Nstep',3X,'Radius',
     +         5X,'Mass', 7X,'Kin E',6X,'Grav E',7X,'Tot E')
32    FORMAT  (3X,' cm ',4X,'gm/cm^-3',1X,'     ',3X,'  cm  ',
     +         5X,' gm ', 7X,9X,'10^51 ergs ',4X,'     ')
33    FORMAT  (3X,'----',4X,'--------',1X,'-----',3X,'------',
     +         5X,'----', 7X,'-----',6X,'------',7X,'-----')
C
      RETURN                          
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,NLINES,R,M,DR,RHOCEN,EKINT,EGRAV,IR)
C writes results for one model to the requested unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P2'
C Input variables:
      INTEGER MUNIT                 !unit to which we are writing
      INTEGER NLINES                !number of lines printed out (I/O)
      REAL M,RHO                    !current mass and density
      REAL RHOCEN                   !central density
      REAL R                        !current radius
      REAL DR                       !radial step
      REAL EGRAV                    !gravitational energy of star
      REAL EKINT                    !kinetic energy of star
      INTEGER IR                    !number of radial steps
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     if screen is full, clear screen and retype headings          
      IF ((MOD(NLINES,TRMLIN-4) .EQ. 0) 
     +                          .AND. (MUNIT .EQ. OUNIT)) THEN
         CALL PAUSE('to continue',1)
         CALL CLEAR
         WRITE (MUNIT,31)
         WRITE (MUNIT,32)
         WRITE (MUNIT,33)
      END IF
C
C     write unscaled values to screen
      WRITE (MUNIT,35) DR*R0,RHOCEN*RHO0,IR,
     +                 R*R0,M*M0,EKINT,EGRAV,EKINT+EGRAV
C
C     keep track of printed lines only for terminal output
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+1
C
31    FORMAT  (3X,'Step',4X,'Cent Den',1X,'Nstep',3X,'Radius',
     +         5X,'Mass', 7X,'Kin E',6X,'Grav E',7X,'Tot E')
32    FORMAT  (3X,' cm ',4X,'gm/cm^-3',1X,'     ',3X,'  cm  ',
     +         5X,' gm ', 7X,9X,'10^51 ergs ',4X,'     ')
33    FORMAT  (3X,'----',4X,'--------',1X,'-----',3X,'------',
     +         5X,'----', 7X,'-----',6X,'------',7X,'-----')
35    FORMAT  (2(1X,1PE8.2),1X,I4,3(2X,1PE9.3),2(2X,1PE10.3))
C      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFOUT(DEVICE,MSTOR,RADIUS,STEP,RMIN,NSTEP,RHOSTR)
C outputs two plots:  1) total mass vs. final radius of star
C                     2) density vs. radius
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables                                                   
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P2'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE                !which device is being used?
      REAL MSTOR(MAXMOD)            !mass of stars
      REAL RADIUS(MAXMOD)           !radius of stars
      REAL STEP(MAXMOD)             !step size
      REAL RMIN(MAXMOD)             !starting radius
      INTEGER NSTEP(MAXMOD)         !number of radial steps
      REAL RHOSTR(MAXMOD,0:MAXSTP)  !density(r) of star
C Local variables
      REAL R,DEN                    !radius and density for one model
      INTEGER IM,IR                 !current model and radius
      DIMENSION R(0:MAXSTP),DEN(0:MAXSTP) !radius and density 
      CHARACTER*9 CYE               !electron fraction as charac string
      INTEGER LEN                   !length of string
      INTEGER SCREEN                !send to terminal
      INTEGER PAPER                 !make a hardcopy
      INTEGER FILE                  !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     define parameters which are the same for both plots
      IF (DEVICE .NE. FILE) THEN
          NPLOT=2                        !how many plots?
C
          XMIN=0.                     
          Y0VAL=0.
C
          ILINE=1                      !line and symbol styles
          ISYM=1
          IFREQ=1
          NXTICK=5
          NYTICK=5
C                                                             
          CALL CONVRT(YE,CYE,LEN)
          TITLE = 'White Dwarf Stars with electron fraction='//CYE
C
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
      END IF
C                                      
      DO 20 IPLOT=1,NPLOT
C       define parameters which are different for each plot
        IF (DEVICE .NE. FILE) THEN
          IF (IPLOT .EQ. 1) THEN
             YMIN=MSTOR(1)
             X0VAL=MSTOR(1)
             XMAX=RADIUS(1)
             YMAX=MSTOR(NMODEL)
             NPOINT=NMODEL
             LABEL(1)= 'radius (10**8 cm)'
             LABEL(2)= 'mass (10**33 gm)'
             CALL LNLNAX                          
          ELSE
             XMAX=RADIUS(1) 
             YMAX=RHO2
             YMIN=RHOCUT*RHO0
             X0VAL=RHOCUT*RHO0
             IFREQ=0
             LABEL(1)= 'radius (10**8 cm)'
             LABEL(2)= 'density (gm*cm**(-3))'
             CALL LGLNAX
          END IF
        END IF                                                    
C 
C       output results
        IF (IPLOT .EQ. 1) THEN
C         total mass vs. final radius
          IF (DEVICE .EQ. FILE) THEN
              WRITE (GUNIT,80)
              WRITE (GUNIT,70) (RADIUS(IM),MSTOR(IM),IM=1,NMODEL)
          ELSE
             CALL XYPLOT (RADIUS,MSTOR)
          END IF
C
        ELSE
C         density(r) vs. r for each model
          DO 30 IM=1,NMODEL
             NPOINT=NSTEP(IM)+1
             DO 40 IR=0,NPOINT-1
                R(IR)=RMIN(IM)+IR*STEP(IM)
                DEN(IR)=RHOSTR(IM,IR)
40           CONTINUE            
             IF (DEVICE .EQ. FILE) THEN
                 WRITE (GUNIT,85) IM
                 WRITE (GUNIT,70) (R(IR),DEN(IR),IR=0,NPOINT)
             ELSE
                CALL XYPLOT(R,DEN)
             END IF
30        CONTINUE                             
        END IF
20    CONTINUE                                                      
C
      IPLOT=NPLOT
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)  !end graphing package
      IF (DEVICE .EQ. SCREEN) CALL TMODE        !switch to text mode
C
70    FORMAT (2(5X,E11.3))
80    FORMAT (/,10X,'radius',10X,'mass')
85    FORMAT (/,'  Density vs. radius for model ',I2)
100   FORMAT (/,' Patience, please; output going to a file.')
C
      RETURN
      END
                                                              
