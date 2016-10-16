CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM PROJ8
C     Project 8: Monte Carlo solution of the H2 molecule
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop/ execute once for each set of param
        CALL PARAM        !get input from screen
        CALL ARCHON       !calculate the eigenvalue for this value of S
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON  
C calculates the electronic eigenvalue or energy auto-correlation 
C for a given separation of the protons 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P8'
      INCLUDE 'IO.ALL'
C Local variables:
C     energy has two indices           
C     first index is the level: sweep, group, or total
C     second index is the value: quantity, quant**2, or sigma**2
      REAL ENERGY(3,3)           !energy
      REAL CONFIG(NCOORD)        !configuration
      REAL W                     !weight for single variational config
      REAL WEIGHT(MAXENS)        !weight of ensemble members
      REAL ENSMBL(NCOORD,MAXENS) !ensemble of configurations
      REAL ESAVE(MAXCRR)         !array of local energies for corr
      REAL EPSILN                !local energy of CONFIG
      REAL ACCPT                 !acceptance ratio 
      INTEGER ITHERM             !thermalization index
      INTEGER ISWP,ISMPL         !sweep and sample index
      INTEGER IQUANT             !quantity index
      INTEGER IGRP               !group index
      INTEGER NLINES             !number of lines printed to terminal
      INTEGER MORE               !how many more groups
      INTEGER SWEEP,GROUP,TOTAL  !which level of calculation
      INTEGER VALUE,SQUARE,SIGSQ !which quantity
      INTEGER CORR,EPS           !what is being calculated?
      INTEGER PIMC,VARY          !which method?
C Functions:
      REAL ELOCAL                !local energy
      INTEGER GETINT             !get integer data from screen
      DATA SWEEP,GROUP,TOTAL/1,2,3/
      DATA VALUE,SQUARE,SIGSQ/1,2,3/
      DATA EPS,CORR /1,2/
      DATA VARY,PIMC /1,2/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     output summary of parameters
      IF (TTERM) CALL PRMOUT(OUNIT,NLINES)         
      IF (TFILE) CALL PRMOUT(TUNIT,NLINES)
      IF (GFILE) CALL PRMOUT(GUNIT,NLINES)
C                  
C     generate initial configuration or ensemble of configurations
      IF (METHOD .EQ. VARY) THEN
         CALL INTCFG(CONFIG,W)
      ELSE IF (METHOD .EQ. PIMC) THEN
         CALL INTENS(ENSMBL,WEIGHT,CONFIG)
      END IF
C
C     take thermalization steps
      DO 10 ITHERM=1,NTHERM
        IF (ITHERM .EQ. 1) WRITE (OUNIT,*) ' Thermalizing...'
        IF (ITHERM .EQ. NTHERM) WRITE (OUNIT,*) ' '
        IF (METHOD .EQ. VARY) THEN
           CALL METROP(CONFIG,W,ACCPT)
        ELSE IF (METHOD .EQ. PIMC) THEN
           CALL TSTEP(ENSMBL,WEIGHT,EPSILN)
        END IF
10    CONTINUE
C
      DO 11 IQUANT=1,3                       !zero total sums
         ENERGY(TOTAL,IQUANT)=0.
11    CONTINUE
      ACCPT=0                                !zero acceptance
      MORE=NGROUP                            !initial number of groups
C
15    CONTINUE                               !allow for more groups
        DO 20 IGRP=NGROUP-MORE+1,NGROUP      !loop over groups
C
           DO 21 IQUANT=1,3                  !zero group sums
             ENERGY(GROUP,IQUANT)=0.
21         CONTINUE
C          
           DO 30 ISWP=1,NFREQ*NSMPL          !loop over sweeps
C
             IF (METHOD .EQ. VARY) THEN      !take a Metrop step
                CALL METROP(CONFIG,W,ACCPT)
             ELSE IF (METHOD .EQ. PIMC) THEN !or a time step
                CALL TSTEP(ENSMBL,WEIGHT,EPSILN)
             END IF
C
             IF (MOD(ISWP,NFREQ) .EQ. 0) THEN !sometimes save the energy
                ISMPL=ISWP/NFREQ             
                IF (METHOD .EQ. VARY) EPSILN=ELOCAL(CONFIG)
                ENERGY(GROUP,VALUE)=ENERGY(GROUP,VALUE)+EPSILN
                ENERGY(GROUP,SQUARE)=ENERGY(GROUP,SQUARE)+EPSILN**2
                IF (.NOT. TERSE) THEN
                  IF (TTERM) WRITE (OUNIT,100) IGRP,ISMPL,NSMPL,EPSILN
                  IF (TFILE) WRITE (TUNIT,100) IGRP,ISMPL,NSMPL,EPSILN
100               FORMAT (5X,' Group ',I4, ',  sample ',I4,' of ',I4,5X,
     +                   'Energy =',F9.4)
                END IF
             END IF
C
30        CONTINUE                            !this group is done
          IF (CALC .EQ. CORR) THEN            !save energy for corr
              ESAVE(IGRP)=ENERGY(GROUP,VALUE)
          ELSE                                !or calc averages
              CALL AVERAG(ENERGY,ACCPT,IGRP)
          END IF
C
20      CONTINUE                                   
        IF (CALC .EQ. CORR) CALL CRLTNS(ESAVE)   !calc corr
C
C       allow for more groups, taking care not to exceed array bounds
        MORE=GETINT(10,0,1000,'How many more groups?')
        IF ((CALC .EQ. CORR) .AND. (NGROUP+MORE .GT. MAXCRR)) THEN
            WRITE (OUNIT,200) MAXCRR-NGROUP
200         FORMAT(' You will run out of storage space for '
     +      'corr if you do more than ',I3,' more groups')
            MORE=GETINT(MAXCRR-NGROUP,0,MAXCRR-NGROUP,
     +                    'How many more groups?')
        END IF
      IF (MORE .GT. 0) THEN
           NGROUP=NGROUP+MORE
           NLINES=0                                                     
           IF (TTERM) CALL CLEAR
           GOTO 15
      END IF
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TSTEP(ENSMBL,WEIGHT,EPSILN)
C take a time step using the Path Integral Monte Carlo method
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P8'
C Input/Output variables
      REAL WEIGHT(MAXENS)       !weight of ensemble members (I/O)
      REAL ENSMBL(NCOORD,MAXENS)!ensemble of configurations (I/O)
      REAL EPSILN               !local energy of CONFIG (output)
C Local variables:
      REAL CONFIG(NCOORD)       !configuration
      REAL W                    !weight for single config
      REAL EBAR,WBAR            !ensmble average local energy and weight
      INTEGER IENSEM            !ensemble index
      INTEGER ICOORD            !coordinate index
      REAL NORM                 !normalization of weights
      REAL SHIFT(NCOORD)        !array containing drift vector
C Functions:
      REAL GAUSS                !Gaussian random number   
      REAL ELOCAL               !local energy of the configuration
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      EBAR=0.                    !zero sums
      WBAR=0.
      DO 10 IENSEM=1,NENSEM                        !loop over ensemble
         DO 20 ICOORD=1,NCOORD
            CONFIG(ICOORD)=ENSMBL(ICOORD,IENSEM)   !get a configuration
20       CONTINUE
         CALL DRIFT(CONFIG,SHIFT)                  !calc shifts
         DO 30 ICOORD=1,NCOORD
            CONFIG(ICOORD)=CONFIG(ICOORD)+         !shift configuration
     +                     GAUSS(DSEED)*SQHBDT+SHIFT(ICOORD)
30       CONTINUE
C
         EPSILN=ELOCAL(CONFIG)                     !calculate energy
         WEIGHT(IENSEM)=WEIGHT(IENSEM)*EXP(-EPSILN*DT)  !calc weight
         EBAR=EBAR+WEIGHT(IENSEM)*EPSILN           !update sums
         WBAR=WBAR+WEIGHT(IENSEM)
C       
         DO 40 ICOORD=1,NCOORD
            ENSMBL(ICOORD,IENSEM)=CONFIG(ICOORD)   !save configuration
40       CONTINUE
10    CONTINUE
C
      EPSILN=EBAR/WBAR                    !weighted average energy
      NORM=NENSEM/WBAR
      DO 50 IENSEM=1,NENSEM               !renormalize weights
         WEIGHT(IENSEM)=NORM*WEIGHT(IENSEM)
50    CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE METROP(CONFIG,W,ACCPT)
C take a Metropolis step
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P8'
C Input/Output variables:
      REAL CONFIG(NCOORD)        !configuration 
      REAL W                     !weight for single config 
      REAL ACCPT                 !acceptance ratio 
C Local variables:
      INTEGER ICOORD             !coordinate index
      REAL CSAVE(NCOORD)         !temp storage for last config
      REAL WTRY                  !weight for trial config
C Function:
      REAL PHI                   !total wave function
      REAL RANNOS                !uniform random number
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 10 ICOORD=1,NCOORD
       CSAVE(ICOORD)=CONFIG(ICOORD)        !save previous values
       CONFIG(ICOORD)=CONFIG(ICOORD)+DELTA*(RANNOS(DSEED)-.5)!trial step
10    CONTINUE
      WTRY=PHI(CONFIG)**2                  !trial weight
C
      IF (WTRY/W .GT. RANNOS(DSEED)) THEN  !sometimes accept the step
         W=WTRY                            !save new weight
         ACCPT=ACCPT+1                     !update accpt ratio
      ELSE
        DO 20 ICOORD=1,NCOORD
          CONFIG(ICOORD)=CSAVE(ICOORD)     !or else restore old config
20      CONTINUE
      END IF
      RETURN
      END     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION ELOCAL(CONFIG)
C calculate the local energy for CONFIG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variable:
      INCLUDE 'PARAM.P8'
C Input variables:
      REAL CONFIG(NCOORD)        !configuration
C Local variables:
      REAL TPOP,VPOP             !kinetic and potential contributions
      REAL EECORR                !elec-elec correlation
      REAL CROSS,CROSS1,CROSS2   !cross terms
      REAL ONEE1,ONEE2           !one electron terms
      REAL X1,X2,Y1,Y2,Z1,Z2     !coordinates of 2 electrons
      REAL R1L,R1R,R2L,R2R,R12   !relative distances
      REAL CHI1,CHI2,F           !parts of the wave function
      REAL DOTR1L,DOTR2L,DOTR1R,DOTR2R !dot products with R12
      REAL R12DR1,SR12Z          !temp vars for dot products
      REAL CHI,FDCHI,SDCHI,LAPCHI !atomic orbitals
      REAL FEE,FDFEE,SDFEE,LAPFEE !elec-elec correlations
      REAL DIST                   !Euclidean distance
      REAL R,X,Y,Z                !dummy variables
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     define functions
      CHI(R)=EXP(-R/A)                        !atomic orbital
      FDCHI(R)=-CHI(R)/A                      !its first derivative,
      SDCHI(R)=CHI(R)/A/A                     !second derivative,
      LAPCHI(R)=SDCHI(R)+2*FDCHI(R)/R         !and Laplacian
C
      FEE(R)=EXP(R/(ALPHA*(1+BETA*R)))        !elec-elec correlation
      FDFEE(R)=FEE(R)/(ALPHA*(1.+BETA*R)**2)  !its first,second deriv,
      SDFEE(R)=FDFEE(R)**2/FEE(R)-2.*BETA*FEE(R)/ALPHA/(1+BETA*R)**3
      LAPFEE(R)=SDFEE(R)+2*FDFEE(R)/R         !and Laplacian
C
      DIST(X,Y,Z)=SQRT(X**2+Y**2+Z**2)        !Euclidean distance
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get coordinates and radii
      CALL RADII(X1,X2,Y1,Y2,Z1,Z2,R1L,R2L,R1R,R2R,R12,CONFIG)
C                                                             
C     calculate dot products with R12
      R12DR1=X1*(X1-X2)+Y1*(Y1-Y2)+Z1*(Z1-Z2) !convenient starting place
      SR12Z=S*(Z1-Z2)/2                !useful constant
      DOTR1L=R12DR1+SR12Z              !dot products with R12
      DOTR1R=R12DR1-SR12Z             
      DOTR2L=DOTR1L-R12**2
      DOTR2R=DOTR1R-R12**2
      DOTR1L=DOTR1L/R12/R1L            !dot products of unit vectors
      DOTR2L=DOTR2L/R12/R2L
      DOTR1R=DOTR1R/R12/R1R
      DOTR2R=DOTR2R/R12/R2R
C
      CHI1=CHI(R1R)+CHI(R1L)           !pieces of the total wave function
      CHI2=CHI(R2R)+CHI(R2L)
      F=FEE(R12)
C
      EECORR=2*LAPFEE(R12)/F           !correlation contribution
      ONEE1=(LAPCHI(R1L)+LAPCHI(R1R))/CHI1               !electron one
      ONEE2=(LAPCHI(R2L)+LAPCHI(R2R))/CHI2               !electron two
      CROSS1=(FDCHI(R1L)*DOTR1L+FDCHI(R1R)*DOTR1R)/CHI1  !cross terms
      CROSS2=(FDCHI(R2L)*DOTR2L+FDCHI(R2R)*DOTR2R)/CHI2
      CROSS=2*FDFEE(R12)*(CROSS1-CROSS2)/F
C
      TPOP=-HBM*(EECORR+ONEE1+ONEE2+CROSS)/2                !kinetic
      VPOP=-E2*(1./R1L + 1./R1R + 1./R2L + 1./R2R - 1./R12) !potential
      ELOCAL=TPOP+VPOP                                      !total
C
      RETURN                                                            
      END 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DRIFT(CONFIG,SHIFT)
C calculate the drift vector (SHIFT) for CONFIG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P8'
C Input/Output variables:
      REAL CONFIG(NCOORD)        !configuration (input)
      REAL SHIFT(NCOORD)         !array containing drift vector (output)
C Local variables:
      INTEGER ICOORD              !coordinate index
      REAL X1,X2,Y1,Y2,Z1,Z2      !coordinates of 2 electrons
      REAL R1L,R1R,R2L,R2R,R12    !relative distances
      REAL CHI1,CHI2,F            !parts of the wave function
      REAL CHI,FDCHI,SDCHI,LAPCHI !atomic orbital
      REAL FEE,FDFEE,SDFEE,LAPFEE !elec-elec correlations
      REAL R                      !dummy variables
      REAL FACTA,FACTB,FACTE      !useful factors
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     define functions
      CHI(R)=EXP(-R/A)                        !atomic orbital
      FDCHI(R)=-CHI(R)/A                      !its first derivative,
      SDCHI(R)=CHI(R)/A/A                     !second derivative,
      LAPCHI(R)=SDCHI(R)+2*FDCHI(R)/R         !and Laplacian
C
      FEE(R)=EXP(R/(ALPHA*(1+BETA*R)))        !elec-elec correlation
      FDFEE(R)=FEE(R)/(ALPHA*(1.+BETA*R)**2)  !its first, second deriv,
      SDFEE(R)=FDFEE(R)**2/FEE(R)-2.*BETA*FEE(R)/ALPHA/(1+BETA*R)**3
      LAPFEE(R)=SDFEE(R)+2*FDFEE(R)/R         !and Laplacian
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get coordinates and radii
      CALL RADII(X1,X2,Y1,Y2,Z1,Z2,R1L,R2L,R1R,R2R,R12,CONFIG)
C
      CHI1=CHI(R1R)+CHI(R1L)           !pieces of the total wave function
      CHI2=CHI(R2R)+CHI(R2L)
      F=FEE(R12)
C
      FACTA=HBMDT*(FDCHI(R1L)/R1L+FDCHI(R1R)/R1R)/CHI1  !useful factors
      FACTB=HBMDT*(FDCHI(R1L)/R1L-FDCHI(R1R)/R1R)/CHI1
      FACTE=HBMDT*FDFEE(R12)/F/R12
C
      SHIFT(1)=FACTA*X1+FACTE*(X1-X2)  !shift for electron one
      SHIFT(2)=FACTA*Y1+FACTE*(Y1-Y2)
      SHIFT(3)=FACTA*Z1+FACTE*(Z1-Z2)+FACTB*S2
C
      FACTA=HBMDT*(FDCHI(R2L)/R2L+FDCHI(R2R)/R2R)/CHI2
      FACTB=HBMDT*(FDCHI(R2L)/R2L-FDCHI(R2R)/R2R)/CHI2
C
      SHIFT(4)=FACTA*X2-FACTE*(X1-X2)  !shift for electron two
      SHIFT(5)=FACTA*Y2-FACTE*(Y1-Y2)
      SHIFT(6)=FACTA*Z2-FACTE*(Z1-Z2)+FACTB*S2
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION GAUSS(DSEED)
C returns a Gaussian random number with zero mean and unit variance
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER IGAUSS                 !sum index
      DOUBLE PRECISION DSEED         !random number seed
      REAL RANNOS                    !uniform random number
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      GAUSS=0.                       !sum 12 uniform random numbers
      DO 10 IGAUSS=1,12
         GAUSS=GAUSS+RANNOS(DSEED)
10    CONTINUE
      GAUSS=GAUSS-6.                 !subtract six so that mean=0
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL FUNCTION PHI(CONFIG)
C calculates the total variational wave function for CONFIG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P8'
C Input variables:
      REAL CONFIG(NCOORD)        !configuration
C Local variables:                                               
      REAL X1,X2,Y1,Y2,Z1,Z2     !coordinates of 2 electrons
      REAL R1L,R1R,R2L,R2R,R12   !relative distances
      REAL CHI1R,CHI1L,CHI2R,CHI2L,F  !parts of the wave function
      REAL CHI,FEE               !terms in the wave function
      REAL R                     !radius
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CHI(R)=EXP(-R/A)                   !atomic orbital
      FEE(R)=EXP(R/(ALPHA*(1+BETA*R)))   !electron-electron correlation
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     calculate the radii
      CALL RADII(X1,X2,Y1,Y2,Z1,Z2,R1L,R2L,R1R,R2R,R12,CONFIG)
C      
      CHI1R=CHI(R1R)          !pieces of the total wave function
      CHI1L=CHI(R1L)
      CHI2R=CHI(R2R)
      CHI2L=CHI(R2L)
      F=FEE(R12)
      PHI=(CHI1L +CHI1R)*(CHI2L+CHI2R)*F   !the whole thing
C
      RETURN
      END  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTENS(ENSMBL,WEIGHT,CONFIG)
C generate the ENSMBL at t=0 for PIMC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P8'
C Output variables:
      REAL CONFIG(NCOORD)        !configuration
      REAL WEIGHT(MAXENS)        !weight of ensemble members
      REAL ENSMBL(NCOORD,MAXENS) !ensemble of configurations
C Local variables:
      INTEGER ISTEP              !step index
      INTEGER ICOORD             !coordinate index
      REAL W                     !weight for single config
      REAL ACCPT                 !acceptance ratio
      INTEGER IENSEM             !ensemble index
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INTCFG(CONFIG,W)      !generate a single intial configuration
C
      DO 10 ISTEP=1,20                 !do 20 thermalization steps
         CALL METROP(CONFIG,W,ACCPT)   !using Metropolis algorithm
10    CONTINUE
C
      DO 30 ISTEP=1,10*NENSEM          !generate the ensemble
         CALL METROP(CONFIG,W,ACCPT)   !take a Metrop step
         IF (MOD(ISTEP,10) .EQ. 0) THEN
             IENSEM=ISTEP/10           !save every 10th config
             DO 20 ICOORD=1,NCOORD     
                ENSMBL(ICOORD,IENSEM)=CONFIG(ICOORD)
20           CONTINUE
             WEIGHT(IENSEM)=1.         !set all weights=1
         END IF
30    CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTCFG(CONFIG,W)
C generate a configuration (CONFIG) and calculate its weight (W)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variable:
      INCLUDE 'PARAM.P8'
C Output variables:
      REAL CONFIG(NCOORD)        !configuration
      REAL W                     !weight for single config
C Local variables:
      INTEGER ICOORD             !coordinate index
C Function:
      REAL PHI                   !total wave function
      REAL RANNOS                !uniform random number
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 10 ICOORD=1,NCOORD      !pick configuration at random
         CONFIG(ICOORD)=A*(RANNOS(DSEED)-.5)
10    CONTINUE
      CONFIG(3)=CONFIG(3)+S2    !center elec 1. at right
      CONFIG(6)=CONFIG(6)-S2    !center elec 2. at left 
      W=PHI(CONFIG)**2          !weight=phi**2
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CRLTNS(ESAVE)
C calculate the energy auto-correlations
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM .P8'
      INCLUDE 'IO.ALL'
C Input variables:
      REAL ESAVE(MAXCRR)           !array of local energies for corr.
C Local variables:
      REAL EI,EIK,ESQI,ESQIK,EIEK  !sums
      INTEGER I,K                  !index of ESAVE
      REAL ECORR(0:MAXCRR)         !energy auto-correlations            
      INTEGER NI                   !number of energies in sum
      INTEGER SCREEN               !send to terminal
      INTEGER PAPER                !make a hardcopy
      INTEGER FILE                 !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 10 K=0,NCORR              !loop over correlation lengths
         EI=0.                     !zero sums
         EIK=0.
         ESQI=0.
         ESQIK=0.
         EIEK=0.
         NI=NGROUP-K
         DO 20 I=1,NI
            EI=EI+ESAVE(I)         !calculate sums
            EIK=EIK+ESAVE(I+K)
            ESQI=ESQI+ESAVE(I)**2
            ESQIK=ESQIK+ESAVE(I+K)**2
            EIEK=EIEK+ESAVE(I)*ESAVE(I+K)
20       CONTINUE   
         EI=EI/NI                  !calculate averages
         EIK=EIK/NI
         ESQI=ESQI/NI
         ESQIK=ESQIK/NI
         EIEK=EIEK/NI
         ECORR(K)=(EIEK-EI*EIK)/(SQRT(ESQI-EI**2))/(SQRT(ESQIK-EIK**2))
10    CONTINUE
C
      IF (GTERM) THEN                   !display results
          CALL PAUSE ('to see the energy auto-correlations...',1)
          CALL GRFOUT(SCREEN,ECORR)
      ELSE IF (TTERM) THEN 
          CALL PAUSE ('to see the energy auto-correlations...',1)
          CALL CRROUT(OUNIT,ECORR)
      END IF
      IF (TFILE) CALL CRROUT(TUNIT,ECORR)
      IF (GHRDCP) CALL GRFOUT(PAPER,ECORR)
      IF (GFILE) CALL GRFOUT(FILE,ECORR)
C
      RETURN
      END    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE AVERAG(ENERGY,ACCPT,IGRP)
C calculate group averages, add to totals, print out
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P8'
      INCLUDE 'IO.ALL'
C Input variables:
C     energy has two indices
C     first index is the level: sweep, group, or total
C     second index is the value: quantity, quant**2, or sigma**2
      REAL ENERGY(3,3)           !energy
      INTEGER IGRP               !group index
      REAL ACCPT                 !acceptance ratio 
C Local variables:
      REAL EVALUE                !current average energy
      REAL SIG1,SIG2             !uncertainties in energy
      REAL U                     !total pot energy of the system 
      INTEGER NLINES             !number of lines printed to terminal
      INTEGER SWEEP,GROUP,TOTAL  !which level of calculation
      INTEGER VALUE,SQUARE,SIGSQ !which quantity
      DATA SWEEP,GROUP,TOTAL/1,2,3/
      DATA VALUE,SQUARE,SIGSQ/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     calculate group averages and uncertainties
      ENERGY(GROUP,VALUE)=ENERGY(GROUP,VALUE)/NSMPL
      ENERGY(GROUP,SQUARE)=ENERGY(GROUP,SQUARE)/NSMPL
      ENERGY(GROUP,SIGSQ)= 
     +  (ENERGY(GROUP,SQUARE)-ENERGY(GROUP,VALUE)**2)/NSMPL
      IF (ENERGY(GROUP,SIGSQ) .LT. 0.) ENERGY(GROUP,SIGSQ)=0.
C
C     add to totals
      ENERGY(TOTAL,VALUE)=ENERGY(TOTAL,VALUE)+ENERGY(GROUP,VALUE)
      ENERGY(TOTAL,SQUARE)=ENERGY(TOTAL,SQUARE)+ENERGY(GROUP,SQUARE)
      ENERGY(TOTAL,SIGSQ)=ENERGY(TOTAL,SIGSQ)+ENERGY(GROUP,SIGSQ)
C
C     calculate current grand averages
      EVALUE=ENERGY(TOTAL,VALUE)/IGRP
      SIG1=(ENERGY(TOTAL,SQUARE)/IGRP-EVALUE**2)/IGRP/NSMPL
      IF (SIG1 .LT. 0.) SIG1=0.
      SIG1=SQRT(SIG1)
      SIG2=SQRT(ENERGY(TOTAL,SIGSQ))/IGRP
C
C     calculate total energy of the system
      IF (S .GT. .01) THEN
          U=EVALUE+E2/S+E2/ABOHR
      ELSE
          U=0.
      END IF         
C
      IF (TTERM) CALL TXTOUT(IGRP,ENERGY,EVALUE,SIG1,SIG2,U,ACCPT,OUNIT)
      IF (TFILE) CALL TXTOUT(IGRP,ENERGY,EVALUE,SIG1,SIG2,U,ACCPT,TUNIT)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RADII(X1,X2,Y1,Y2,Z1,Z2,R1L,R2L,R1R,R2R,R12,CONFIG)
C calculates cartesian coordinates and radii given CONFIG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variable:
      INCLUDE 'PARAM.P8'
C Input variables:
      REAL CONFIG(NCOORD)        !configuration
C Output variables:
      REAL X1,X2,Y1,Y2,Z1,Z2     !coordinates of 2 electrons
      REAL R1L,R1R,R2L,R2R,R12   !relative distances
      REAL DIST                  !Euclidean distance
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DIST(X1,Y1,Z1)=SQRT(X1**2+Y1**2+Z1**2)        !Euclidean distance
C
      X1=CONFIG(1)            !give the CONFIG elements their real names
      X2=CONFIG(4)
      Y1=CONFIG(2)
      Y2=CONFIG(5)
      Z1=CONFIG(3)
      Z2=CONFIG(6)
C
      R1L=DIST(X1,Y1,Z1+S2)   !calculate separations
      R1R=DIST(X1,Y1,Z1-S2)
      R2L=DIST(X2,Y2,Z2+S2)
      R2R=DIST(X2,Y2,Z2-S2)
      R12=DIST(X1-X2,Y1-Y2,Z1-Z2)
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
      INCLUDE 'PARAM.P8'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)         
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get environment parameters
      CALL SETUP                
C 
C     display header screen     
      DESCRP(1)= 'PROJECT 8'
      DESCRP(2)= 'Monte Carlo solution of the H2 molecule'
      NHEAD=2
C 
C     text output description
      DESCRP(3)= 'electronic eigenvalue and its uncertainty'
      DESCRP(4)= 'or energy auto-correlation'
      NTEXT=2
C             
C     graphics output description
      DESCRP(5)= 'energy auto-correlation vs. correlation length'
      NGRAPH=1
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C 
C     calculate constants 
      HBM=7.6359               !hbar*2/(mass)
      E2=14.409                !charge of the electron
      ABOHR=HBM/E2
C
C     setup menu arrays, beginning with constant part
      CALL MENU                      
C                
      MTYPE(12)=TITLE
      MPRMPT(12)= 'PHYSICAL PARAMETERS'
      MLOLIM(12)=0.
      MHILIM(12)=1.
C                
      MTYPE(13)=FLOAT                                        
      MPRMPT(13)='Enter the interproton separation S (Angstroms)'
      MTAG(13)='Inter proton separation (Angstroms)'
      MLOLIM(13)=0.
      MHILIM(13)=10.
      MREALS(13)=0.
C
      MTYPE(14)=FLOAT     
      MPRMPT(14)=
     + 'Enter value for variational parameter Beta (Angstroms**-1)'
      MTAG(14)='variational parameter Beta (Angstroms**-1)'
      MLOLIM(14)=0.
      MHILIM(14)=10.
      MREALS(14)=.25
C
      MTYPE(15)=SKIP
      MREALS(15)=35.
C 
      MTYPE(37)=TITLE
      MPRMPT(37)= 'NUMERICAL PARAMETERS'
      MLOLIM(37)=1.
      MHILIM(37)=1.
C                    
      MTYPE(38)=TITLE
      MPRMPT(38)='Methods of calculation:'
      MLOLIM(38)=0.
      MHILIM(38)=0.
C
      MTYPE(39)=MTITLE
      MPRMPT(39)='1) Variational'
      MLOLIM(39)=0.
      MHILIM(39)=0.
C
      MTYPE(40)=MTITLE
      MPRMPT(40)='2) Path Integral Monte Carlo'
      MLOLIM(40)=0.
      MHILIM(40)=1.
C
      MTYPE(41)=MCHOIC
      MPRMPT(41)='Make a menu choice and press return'
      MTAG(41)='44 42'
      MLOLIM(41)=1.
      MHILIM(41)=2.
      MINTS(41)=1
      MREALS(41)=1.
C
      MTYPE(42)=NUM
      MPRMPT(42)= 'Enter size of the ensemble'
      MTAG(42)= 'Ensemble size'
      MLOLIM(42)=1.
      MHILIM(42)=MAXENS
      MINTS(42)=20.
C            
      MTYPE(43)=FLOAT
      MPRMPT(43)='Enter time step (units of 1E-16 sec/hbar)'
      MTAG(43)='Time step (units of 1E-16 sec/hbar)'
      MLOLIM(43)=0.
      MHILIM(43)=10.
      MREALS(43)=.01
C
      MTYPE(44)=FLOAT
      MPRMPT(44)= 'Enter step size for sampling PHI (Angstroms)'
      MTAG(44)= 'Sampling step size (Angstroms)'
      MLOLIM(44)=.01
      MHILIM(44)=10.
      MREALS(44)=.4
C 
      MTYPE(45)=NUM
      MPRMPT(45)= 'Number of thermalization sweeps'
      MTAG(45)= 'Thermalization sweeps'
      MLOLIM(45)=0
      MHILIM(45)=1000
      MINTS(45)=20
C 
      MTYPE(46)=TITLE
      MPRMPT(46)='Quantity to calculate:'
      MLOLIM(46)=1.
      MHILIM(46)=0.
C
      MTYPE(47)=MTITLE
      MPRMPT(47)='1) Energy'
      MLOLIM(47)=0.
      MHILIM(47)=0.
C
      MTYPE(48)=MTITLE
      MPRMPT(48)='2) Correlations'
      MLOLIM(48)=0.
      MHILIM(48)=1.
C               
      MTYPE(49)=MCHOIC
      MPRMPT(49)='Make a menu choice and press return'
      MTAG(49)='50 53'
      MLOLIM(49)=1.
      MHILIM(49)=2.
      MINTS(49)=1
      MREALS(49)=1.
C                
      MTYPE(50)=NUM
      MPRMPT(50)=  'Enter sampling frequency (to avoid correlations)'
      MTAG(50)= 'Sampling frequency'
      MLOLIM(50)=1
      MHILIM(50)=100
      MINTS(50)=6
C 
      MTYPE(51)=NUM
      MPRMPT(51)= 'Enter number of samples in a group'
      MTAG(51)= 'Group sample size'
      MLOLIM(51)=1
      MHILIM(51)=1000
      MINTS(51)=10
C
      MTYPE(52)=SKIP
      MREALS(52)=54.
C 
      MTYPE(53)=NUM
      MPRMPT(53)= 'Enter maximum correlation length'
      MTAG(53)= 'Maximum correlation length'
      MLOLIM(53)=1
      MHILIM(53)=100.
      MINTS(53)=40
C                  
      MTYPE(54)=NUM
      MPRMPT(54)= 'Enter number of groups'
      MTAG(54)= 'Number of groups'
      MLOLIM(54)=1
      MHILIM(54)=1000
      MINTS(54)=10
C 
      MTYPE(55)=NUM
      MPRMPT(55)= 'Integer random number seed for init fluctuations'
      MTAG(55)= 'Random number seed'
      MLOLIM(55)=1000.
      MHILIM(55)=99999.
      MINTS(55)=34767
C 
      MTYPE(56)=SKIP 
      MREALS(56)=60.
C 
      MSTRNG(MINTS(75))= 'proj8.txt'
C
      MTYPE(76)=BOOLEN
      MPRMPT(76)='Do you want the SHORT version of the output?'
      MTAG(76)='Short version of output'
      MINTS(76)=0
C                      
      MTYPE(77)=SKIP
      MREALS(77)=80.
C                                                 
      MSTRNG(MINTS(86))= 'proj8.grf'
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
      INCLUDE 'PARAM.P8'
C Local variables:
      REAL AOLD           !temp variable to search for A
      INTEGER CORR,EPS    !what is being calculated?
      INTEGER PIMC,VARY   !which method?
C gives the map between menu indices and parameters
      INTEGER IS,IBETA,IMETHD,INENSM,IDT,IDELTA,ITHERM,ICALC,
     +        INFREQ,INSMPL,INCORR,IGROUP,IDSEED,ITERSE
      PARAMETER (IS     = 13)
      PARAMETER (IBETA  = 14)
      PARAMETER (IMETHD = 41)
      PARAMETER (INENSM = 42)
      PARAMETER (IDT    = 43)
      PARAMETER (IDELTA = 44)
      PARAMETER (ITHERM = 45)
      PARAMETER (ICALC  = 49)
      PARAMETER (INFREQ = 50)
      PARAMETER (INSMPL = 51)
      PARAMETER (INCORR = 53)
      PARAMETER (IGROUP = 54)
      PARAMETER (IDSEED = 55)
      PARAMETER (ITERSE = 76)
C Functions:
      LOGICAL LOGCVT      !converts 1 and 0 to true and false 
      INTEGER GETINT      !get integer from screen   
      DATA VARY,PIMC /1,2/
      DATA EPS,CORR /1,2/
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
      S=MREALS(IS)
      BETA=MREALS(IBETA)
      METHOD=MINTS(IMETHD)
      NENSEM=MINTS(INENSM)
      DT=MREALS(IDT)
      DELTA=MREALS(IDELTA)
      NTHERM=MINTS(ITHERM)
      CALC=MINTS(ICALC)
      NFREQ=MINTS(INFREQ)
      NSMPL=MINTS(INSMPL)
      NCORR=MINTS(INCORR)
      NGROUP=MINTS(IGROUP)
      DSEED=DBLE(MINTS(IDSEED))
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
C     check parameters for correlations, fix NFREQ, NSMPL
      IF (CALC .EQ. CORR) THEN
         NFREQ=1            !fixed for correlations
         NSMPL=1
         IF ((NGROUP .GT. MAXCRR) .OR. ((NGROUP-NCORR) .LE. 20)) THEN
           WRITE (OUNIT,*) ' '
           WRITE (OUNIT,20)
           WRITE (OUNIT,30) NGROUP,NCORR+20,MAXCRR
20         FORMAT (5X,' For reasonable values of the correlations ')
30         FORMAT (5X,' NGROUP (',I4,') must be between NCORR+20 (', 
     +             I4,') and MAXCRR (',I4,')')
           WRITE (OUNIT,*) ' '
           NCORR=GETINT(NCORR,1,100,'Reenter NCORR')
           NGROUP=GETINT(NCORR+100,NCORR+20,MAXCRR,'Re-enter NGROUP')
           MINTS(INCORR)=NCORR
           MINTS(IGROUP)=NGROUP
         END IF
      END IF
C
      CALL CLEAR
C
C     calculate derivative parameters
      A=ABOHR
      AOLD=0.
10    IF (ABS(A-AOLD) .GT. 1.E-6) THEN
         AOLD=A
         A=ABOHR/(1+EXP(-S/AOLD))
      GOTO 10
      END IF
      S2=S/2
      HBMDT=HBM*DT
      SQHBDT=SQRT(HBMDT)
      ALPHA=2*ABOHR
      IF (METHOD .EQ. PIMC) DELTA=1.5*A      
C
      RETURN               
      END    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRMOUT(MUNIT,NLINES)
C write out parameter summary to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P8'
C Input variables:
      INTEGER MUNIT                   !fortran unit number
      INTEGER NLINES                  !number of lines sent to terminal
C Local variables:
      INTEGER CORR,EPS                !what is being calculated?
      INTEGER PIMC,VARY               !which method?
      DATA EPS,CORR /1,2/
      DATA VARY,PIMC /1,2/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) THEN
         CALL CLEAR
      ELSE
         WRITE (MUNIT,*) ' '
         WRITE (MUNIT,*) ' '
      END IF
C
      WRITE (MUNIT,5)
      WRITE (MUNIT,7) S
      WRITE (MUNIT,8) BETA
      WRITE (MUNIT,9) A
      WRITE (MUNIT,*) ' '
      IF (METHOD .EQ. PIMC) THEN
         WRITE (MUNIT,10) NENSEM, DT
      ELSE 
         WRITE (MUNIT,11)
      END IF
      IF (CALC .EQ. CORR) WRITE (MUNIT,12) NCORR
      WRITE (MUNIT,13)DELTA
      WRITE (MUNIT,15) NTHERM
      WRITE (MUNIT,20) NFREQ,NSMPL
      WRITE (MUNIT,*) ' '
C
      NLINES=11
C      
5     FORMAT (' Output from project 8:',
     +        ' Monte Carlo solution of the H2 molecule')
7     FORMAT (' Proton separation (Angstroms) = ',F7.4)
8     FORMAT (' Variational parameter Beta (Angstroms**-1) = ',F7.4)
9     FORMAT (' Wave function parameter A (Angstroms) = ',F7.4)

10    FORMAT (' Path Integral Monte Carlo with ensemble size = ',I4,
     +        ' and time step = ',1PE12.5)
11    FORMAT (' Variational Monte Carlo method')
12    FORMAT (' correlations will be calculated up to K = ', I4)
13    FORMAT (' Metropolis step in coordinate space (Angstroms)=',F7.4)
15    FORMAT (' number of thermalization sweeps =',I4)
20    FORMAT (' sweep frequency = ',I4,' group size =',I4)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(IGRP,ENERGY,EVALUE,SIG1,SIG2,U,ACCPT,MUNIT)
C write out results to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P8'
      INCLUDE 'IO.ALL'
C Input variables:
C     energy has two indices
C     first index is the level: sweep, group, or total
C     second index is the value: quantity, quant**2, or sigma**2
      REAL ENERGY(3,3)           !energy
      INTEGER IGRP               !group index
      REAL EVALUE                !current average energy
      REAL SIG1,SIG2             !uncertainties in energy
      REAL U                     !total energy of the system at this S
      REAL ACCPT                 !acceptance ratio 
      INTEGER MUNIT              !unit to write to
C Local variables:
      INTEGER SWEEP,GROUP,TOTAL  !which level of calculation
      INTEGER VALUE,SQUARE,SIGSQ !which quantity
      INTEGER PIMC,VARY          !which method?
      DATA SWEEP,GROUP,TOTAL/1,2,3/
      DATA VALUE,SQUARE,SIGSQ/1,2,3/
      DATA VARY,PIMC /1,2/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE (MUNIT,10) IGRP,NGROUP,
     +        ENERGY(GROUP,VALUE),SQRT(ENERGY(GROUP,SIGSQ))
      IF (METHOD .EQ. VARY) THEN
        WRITE (MUNIT,20) EVALUE,SIG1,SIG2,U,ACCPT/IGRP/NFREQ/NSMPL
      ELSE
        WRITE (MUNIT,30) EVALUE,SIG1,SIG2,U
      END IF
      IF (MUNIT .EQ. TUNIT) WRITE (MUNIT,*) ' ' 
C
      IF ((MUNIT .EQ. OUNIT) .AND. (.NOT. TERSE)) 
     +      CALL PAUSE('to continue...',1)

10    FORMAT (2X,'Group ', I4,' of ', I4,5X,'Eigenvalue = ',F9.4,
     +       ' +- ',F8.4)
20    FORMAT (2X,'Grand average E =',F9.4,'+-',F8.4,'/',F8.4,
     +        '   U=',F9.4,'  acceptance=',F6.4)
30    FORMAT (2X,'Grand average E =',F9.4,'+-',F8.4,'/',F8.4,
     +        '   U=',F9.4)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFOUT(DEVICE,ECORR)
C outputs energy auto-correlation vs. correlation length
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables                                                   
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P8'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      REAL ECORR(0:MAXCRR)         !energy auto-correlations            
      INTEGER DEVICE               !which device is being used?
C Local variables                  
      REAL K(0:MAXCRR)             !correlation length
      INTEGER IK                   !correlation length index
      CHARACTER*9 CB,CS,CG         !Beta, S, NGROUP as character data
      INTEGER SCREEN               !send to terminal
      INTEGER PAPER                !make a hardcopy
      INTEGER FILE                 !send to a file
      INTEGER LB,LS,LG             !true length of character data
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     calculate parameters for graphing
      IF (DEVICE .NE. FILE) THEN
          NPLOT=1                      !how many plots?
          IPLOT=1
C
          YMIN=-1.                     !limits on plot
          YMAX=1.
          XMIN=0.
          XMAX=NCORR
          X0VAL=0.
          Y0VAL=XMIN
C 
          NPOINT=NCORR+1
C 
          ILINE=1                      !line and symbol styles
          ISYM=1
          IFREQ=1
          NXTICK=5
          NYTICK=5
C
          CALL CONVRT(BETA,CB,LB)             !titles and labels
          CALL CONVRT(S,CS,LS)
          CALL ICNVRT(NGROUP,CG,LG)
          INFO='NGROUP = '//CG(1:LG)
          TITLE = 'H2 molecule, S='//CS(1:LS)//',  Beta='//CB(1:LB)
          LABEL(1)= 'Correlation length'
          LABEL(2)= 'Energy auto-correlation'
C             
          CALL GTDEV(DEVICE)                  !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE  !change to graphics mode
          CALL LNLNAX                         !draw axes
      END IF
C
      DO 10 IK=0,NCORR                        !fill array of corr length
         K(IK)=REAL(IK)
10    CONTINUE
C                                                      
C     output results
      IF (DEVICE .EQ. FILE) THEN
          WRITE (GUNIT,*) ' '
          WRITE (GUNIT,25) NGROUP
          WRITE (GUNIT,70) (K(IK),ECORR(IK),IK=0,NCORR)
      ELSE
          CALL XYPLOT (K,ECORR)
      END IF
C
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)!close graphics package
      IF (DEVICE .EQ. SCREEN) CALL TMODE      !switch to text mode
C
70    FORMAT (2(5X,E11.3))
25    FORMAT (6X,'corr length',5X,
     +    'energy auto-correlation for NGROUP=',I5)
100   FORMAT (/,' Patience, please; output going to a file.')

      RETURN                                                    
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CRROUT(MUNIT,ECORR)
C write out correlations to MUNIT                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM .P8'
      INCLUDE 'IO.ALL'
C Input variables:
      REAL ECORR(0:MAXCRR)         !energy auto-correlations            
      INTEGER MUNIT                !unit to write to
C Local variables:
      INTEGER K                    !correlation length
      INTEGER NLINES               !number of lines written to screen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) THEN
         CALL CLEAR
      ELSE
         WRITE (MUNIT,*) ' '
      END IF
C
      NLINES=1                
      WRITE (MUNIT,30) NGROUP
30    FORMAT(' Correlations with NGROUP = ',I5)
      DO 10 K=0,NCORR
         NLINES=NLINES+1
         WRITE (MUNIT,20) K,ECORR(K)
         IF ((MUNIT .EQ. OUNIT) .AND. (MOD(NLINES,TRMLIN-3).EQ. 0)) THEN
            CALL PAUSE('to continue...',0)
            NLINES=0
         END IF
10    CONTINUE
      IF (MUNIT .NE. OUNIT) WRITE (MUNIT,*) ' '
20    FORMAT (5X,' Correlation length = ', I3, 5X, 
     +        'Energy auto-correlation = ', F12.5)
      RETURN
      END
                                                              
