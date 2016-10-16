CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM PROJ3  
C     Project 3:  Hartree-Fock solutions of small atomic systems in the
C                 filling approximation
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop/ execute once for each set of param
        CALL PARAM        !get input from screen
        CALL ARCHON       !find the Hartree-Fock wave functions
      GOTO 5                          
      END                                                     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON
C find the Hartree-Fock wave functions for the specified atom
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P3'
C Local variables:
      REAL E(MAXSTT+1,8)              !all energies of all states
      REAL FOCK(0:MAXSTP,MAXSTT)      !Fock terms
      REAL RHO(0:MAXSTP)              !density
      REAL PSTOR(0:MAXSTP,MAXSTT)     !radial wave function
      REAL PHI(0:MAXSTP)              !electron potential
      REAL ESP                        !single particle energy of state
      INTEGER ITER                    !iteration index
      INTEGER STATE                   !single particle state index
      REAL ZSTAR                      !optimal effective nuclear charge
      INTEGER DEVICE                  !current graphing device
      INTEGER ISTOP,ISTART            !current limits on iteration
      INTEGER NLINES                  !number of lines written to screen
      INTEGER SCREEN                  !send to terminal
      INTEGER PAPER                   !make a hardcopy
      INTEGER FILE                    !send to a file
C Functions
      INTEGER GETINT                  !integer screen input
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     begin iterations with a good guess
      MIX=1.                    !no old density to mix with new
      ZSTAR=Z                   
      CALL HYDRGN(ZSTAR,PSTOR)  !find hydrogenic wave functions
      CALL ENERGY(E,FOCK,RHO,PHI,PSTOR)  !and energy
C     optimal ZSTAR using Virial theorem
      ZSTAR=-Z*(E(NSTATE+1,IVTOT)/(2*E(NSTATE+1,IKTOT)) )
      CALL HYDRGN(ZSTAR,PSTOR)  !find new hydrogenic wave functions 
      CALL ENERGY(E,FOCK,RHO,PHI,PSTOR) ! and energies
C
      !output summary of parameters
      IF (TTERM) CALL PRMOUT(OUNIT,ZSTAR,NLINES)
      IF (TFILE) CALL PRMOUT(TUNIT,ZSTAR,NLINES)
      IF (GFILE) CALL PRMOUT(GUNIT,ZSTAR,NLINES)
C 
      MIX=.5                    !mix old and new density for stability
      ITER=0                    !zero iteration counters
      ISTART=0                                        
      ISTOP=0
C 
C     output initial energies
      IF (TTERM) CALL TXTOUT(OUNIT,ITER,E,NLINES)
      IF (TFILE) CALL TXTOUT(TUNIT,ITER,E,NLINES)
C
10    CONTINUE                  !loop over iterations
        ISTART=ISTOP+1          !update iteration counters
        ISTOP=ISTOP+NITER
        DO 100 ITER=ISTART,ISTOP
           DO 50 STATE=1,NSTATE
              IF (NOCC(STATE) .NE. 0) THEN  !loop over all occpd states
                  ESP=E(STATE,ITOT)         !single particle energy
                  IF (ESP .GT. 0) ESP=-10   !keep particle bound
C                 find single part wave function
                  CALL SNGWFN(STATE,ESP,PSTOR,FOCK,PHI)
              END IF
50         CONTINUE
C          calculate new energies and output
           CALL ENERGY(E,FOCK,RHO,PHI,PSTOR)
           IF (TTERM) CALL TXTOUT(OUNIT,ITER,E,NLINES)
           IF (TFILE) CALL TXTOUT(TUNIT,ITER,E,NLINES)
 100     CONTINUE
C
C        allow for more iterations
         NITER=GETINT(0,0,20,'How many more iterations?')
      IF (NITER .NE. 0 ) GOTO 10
C
      IF (TTERM) CALL CLEAR
C
C     graphics output
      IF (GTERM) CALL GRFOUT(SCREEN,PSTOR,RHO,E)
      IF (GFILE) CALL GRFOUT(FILE,PSTOR,RHO,E)
      IF (GHRDCP) CALL GRFOUT(PAPER,PSTOR,RHO,E)
C               
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HYDRGN(ZSTAR,PSTOR)
C creates hydrogenic wave functions PSTOR with nuclear charge=ZSTAR      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P3'
C Input variables:
      REAL ZSTAR                             !effective nuclear charge
      REAL PSTOR(0:MAXSTP,MAXSTT)            !radial wave function
C Local variables:
      INTEGER IR                             !radial index
      REAL RSTAR                             !scaled radius
      REAL ERSTAR                            !useful exponential
      INTEGER STATE                          !state index
      REAL NORM                              !norm of wave function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     store radial parts of hydrogenic wave functions
      DO 20 IR=0,NR
         RSTAR=IR*DR*ZSTAR/ABOHR              !scaled radius
         ERSTAR=EXP(-RSTAR/2.)                !useful exponential
         IF (NOCC(1) .NE. 0) PSTOR(IR,1)=RSTAR*ERSTAR**2
         IF (NOCC(2) .NE. 0) PSTOR(IR,2)=(2-RSTAR)*RSTAR*ERSTAR
         IF (NOCC(3) .NE. 0) PSTOR(IR,3)=(RSTAR**2)*ERSTAR
20     CONTINUE
C            
C      normalize wave functions
       DO 40 STATE=1,NSTATE
          IF (NOCC(STATE) .NE. 0) THEN
              NORM=0.
              DO 30 IR=0,NR
                  NORM=NORM+PSTOR(IR,STATE)**2
30            CONTINUE
              NORM=1./(SQRT(NORM*DR))
              DO 35 IR=0,NR
                  PSTOR(IR,STATE)=PSTOR(IR,STATE)*NORM
35            CONTINUE
          END IF
40     CONTINUE
C 
       RETURN
       END                                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ENERGY(E,FOCK,RHO,PHI,PSTOR)
C subroutine to calculate the energies of a normalized
C set of single-particle wave functions (PSTOR); 
C also calculates Fock terms, density, and electron potential
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P3'
C Passed variables:
      REAL PSTOR(0:MAXSTP,MAXSTT)     !radial wave function (input)
      REAL FOCK(0:MAXSTP,MAXSTT)      !Fock terms (output)
      REAL PHI(0:MAXSTP)              !electron potential (output)
      REAL RHO(0:MAXSTP)              !density (output)
      REAL E(MAXSTT+1,8)       !all energies of all states (output)
C Local variables:
      INTEGER STATE            !state index
      INTEGER IR               !radial index
      INTEGER IE               !energy index
      REAL R                   !current radius
      REAL LL1                 !square of angular momentum
      REAL PM,PZ,PZ2           !values to calc d(rho)/dr
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL SOURCE(PSTOR,FOCK,RHO)  !calculate Fock terms and density
      CALL POISSN(PHI,RHO)         !calc potntl due to electron charge
C
      DO 20 IE=IKEN,ITOT           !zero total energies
         E(NSTATE+1,IE)=0
20    CONTINUE
C 
      DO 38 STATE=1,NSTATE
         IF (NOCC(STATE) .NE. 0) THEN    
C
            DO 10 IE=IKEN,ITOT     !zero energy for this state
               E(STATE,IE)=0
10          CONTINUE
C       
            LL1=ANGMOM(STATE)*(ANGMOM(STATE)+1)
            PM=0
            DO 48 IR=1,NR          !integrate the energy densities
              R=IR*DR
              PZ=PSTOR(IR,STATE)
              PZ2=PZ**2
C
              E(STATE,IKEN)=E(STATE,IKEN)+(PZ-PM)**2        !kinetic
              E(STATE,ICEN)=E(STATE,ICEN)+PZ2*LL1/R**2      !centrifugal
              E(STATE,IVEN)=E(STATE,IVEN)-PZ2/R             !elec-nucl
              E(STATE,IVEE)=E(STATE,IVEE)+PHI(IR)*PZ2       !elec-elec
              E(STATE,IVEX)=E(STATE,IVEX)+FOCK(IR,STATE)*PZ !exchange
C
               PM=PZ     !roll values for derivative
48          CONTINUE
C
C           put in constant factors
            E(STATE,IKEN)=E(STATE,IKEN)*HBARM/(2*DR)
            E(STATE,ICEN)=E(STATE,ICEN)*DR*HBARM/2       
            E(STATE,IVEN)=E(STATE,IVEN)*ZCHARG*DR
            E(STATE,IVEE)=E(STATE,IVEE)*DR
            E(STATE,IVEX)=E(STATE,IVEX)*DR  
C 
C           calculate totals for this state
            E(STATE,IKTOT)=E(STATE,IKEN)+E(STATE,ICEN)
            E(STATE,IVTOT)=E(STATE,IVEN)+E(STATE,IVEE)+E(STATE,IVEX)
            E(STATE,ITOT)=E(STATE,IVTOT)+E(STATE,IKTOT)
C 
C           add this state's contribution to total energy
            DO 30 IE=IKEN,IVEX
               E(NSTATE+1,IE)=E(NSTATE+1,IE)+E(STATE,IE)*NOCC(STATE)
30          CONTINUE
C   
         END IF
38    CONTINUE
C 
      STATE=NSTATE+1              !calculate total energies
C     don't double count electron-electron and exchange energies
      E(STATE,IVEE)=E(STATE,IVEE)/2
      E(STATE,IVEX)=E(STATE,IVEX)/2
C          
C     find total kinetic, potential and total energies
      E(STATE,IKTOT)=E(STATE,IKEN)+E(STATE,ICEN)
      E(STATE,IVTOT)=E(STATE,IVEN)+E(STATE,IVEE)+E(STATE,IVEX)
      E(STATE,ITOT)=E(STATE,IVTOT)+E(STATE,IKTOT)
C
      RETURN
      END                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SOURCE(PSTOR,FOCK,RHO)
C subroutine to compute the density and the Fock terms given PSTOR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P3'
C Passed variables:
      REAL PSTOR(0:MAXSTP,MAXSTT)!radial wave function (input)
      REAL RHO(0:MAXSTP)         !density (output)
      REAL FOCK(0:MAXSTP,MAXSTT) !Fock terms (output)
C Local variables:
      REAL DF                    !increment in Fock term
      REAL FAC                   !constant factor in Fock term
      INTEGER IR                 !radial index
      REAL R                     !current radius
      REAL RLAM                  !r**lambda
      REAL RLAM1                 !r**(lambda+1)
      INTEGER STATE              !indexes state
      INTEGER STATE2             !indexes second state in Fock term
      REAL SUM                   !sum for Fock integrals
      REAL TERM                  !temporary term for sum
      INTEGER L1,L2              !angular momentum for two states
      INTEGER LSTART,LSTOP       !limits on sum of L1+L2
      INTEGER LAM                !current value of ang mom sum
      REAL THREEJ                !three j coefficient
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 20 IR=1,NR              !include fraction of old density
           RHO(IR)=(1.-MIX)*RHO(IR)
20    CONTINUE
C
      DO 30 STATE=1,NSTATE
         IF (NOCC(STATE) .NE. 0) THEN     !loop over occupied states
C
            DO 40 IR=1,NR     !contribution of this state to density
               RHO(IR)=RHO(IR)+MIX*NOCC(STATE)*(PSTOR(IR,STATE)**2)
               FOCK(IR,STATE)=0           !zero FOCK term
40          CONTINUE
C
C           begin calculation of Fock term
            DO 50 STATE2=1,NSTATE
               IF (NOCC(STATE2) .NE. 0) THEN  !loop over occupied states
C
                  L1=ANGMOM(STATE)
                  L2=ANGMOM(STATE2)
                  LSTART=IABS(L1-L2)          !limits on ang mom sum
                  LSTOP=L1+L2                 
C
                  DO 60 LAM=LSTART,LSTOP,2    !loop over ang mom values
C
                     CALL SQR3J(L1,L2,LAM,THREEJ)
                     FAC=-CHARGE/2*NOCC(STATE2)*THREEJ
C
                     SUM=0
                     DO 80 IR=1,NR          !outward integral for Fock
                        R=IR*DR
                        RLAM=R**LAM
                        TERM=PSTOR(IR,STATE2)*PSTOR(IR,STATE)*RLAM/2
                        SUM=SUM+TERM
                        DF=PSTOR(IR,STATE2)*FAC*SUM*DR/(RLAM*R)
                        FOCK(IR,STATE)=FOCK(IR,STATE)+DF
                        SUM=SUM+TERM
   80                CONTINUE
C
                     SUM=0
                     DO 90 IR=NR,1,-1       !inward integral for Fock
                        R=IR*DR
                        RLAM1=R**(LAM+1)
                        TERM=PSTOR(IR,STATE2)*PSTOR(IR,STATE)/RLAM1/2
                        SUM=SUM+TERM
                        DF=PSTOR(IR,STATE2)*FAC*SUM*DR*RLAM1/R
                        FOCK(IR,STATE)=FOCK(IR,STATE)+DF
                        SUM=SUM+TERM
90                   CONTINUE
C
60                CONTINUE  !end loop over lam
                END IF      
50          CONTINUE        !end loop over state2
         END IF
30    CONTINUE              !end loop over state
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE POISSN(PHI,RHO)
C subroutine to solve Poisson's equation for the direct potential
C given the electron density RHO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P3'
C Passed variables:
      REAL PHI(0:MAXSTP)               !electron potential (output)
      REAL RHO(0:MAXSTP)               !density (input)
C Local variables:
      REAL CON                         !useful constant
      INTEGER IR                       !radial index
      REAL M                           !linear behavior to subtract
      REAL R                           !radius
      REAL SM,SP,SZ                    !source terms
      REAL SUM                         !total charge
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUM = 0.                         !quadrature of density to get 
      DO 19 IR = 1,NR                  ! initial value for PHI(1)
          SUM = SUM+RHO(IR)/REAL(IR)
19    CONTINUE
C 
      CON = DR**2/12                   !initial values for outward integ
      SM = 0
      SZ = -CHARGE*RHO(1)/DR
      PHI(0) = 0
      PHI(1) = SUM*CHARGE*DR           
C
      DO 29 IR = 1,NR-1                !Numerov algorithm
            SP=-CHARGE*RHO(IR+1)/((IR+1)*DR)
            PHI(IR+1)=2*PHI(IR)-PHI(IR-1)+CON*(10*SZ+SP+SM)
            SM=SZ
            SZ=SP
29    CONTINUE
C 
      M=(PHI(NR)-PHI(NR-10))/(10*DR)   !subtract off linear behavior
      DO 39 IR=1,NR
            R=IR*DR
            PHI(IR)=PHI(IR)/R-M        !factor of 1/r for true potl
39    CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SNGWFN(STATE,E,PSTOR,FOCK,PHI)
C subroutine to solve the single-particle wave function as an
C inhomogeneous boundary-value problem for a given state, energy E,
C source term (FOCK), and potential (PHI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P3'
C Passed variables:
      INTEGER STATE                 !which single particle state (input)
      REAL E                        !single particle energy (input)
      REAL FOCK(0:MAXSTP,MAXSTT)    !Fock terms (input)
      REAL PHI(0:MAXSTP)            !electron potential (input)
      REAL PSTOR(0:MAXSTP,MAXSTT)   !radial wave function (output)
C Local variables:
      REAL PSIIN(0:MAXSTP),PSIOUT(0:MAXSTP)  !homogeneous solutions
      INTEGER NR2                   !midpoint for the lattice
      REAL LL1                      !angular momentum
      REAL DRHBM                    !useful constant
      REAL K2M,K2Z,K2P              !local wave numbers
      REAL NORM                     !normalization
      REAL R                        !current radius
      REAL SUM                      !sum for integration
      REAL TERM                     !temp term in integration
      REAL WRON                     !Wronskian at middle of lattice
      INTEGER IR                    !radial index
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DRHBM=DR**2/HBARM/6           !useful constants
      LL1=ANGMOM(STATE)*(ANGMOM(STATE)+1)*HBARM/2                
C     
      K2M=0                         !integrate outward homogeneous soln
      K2Z=DRHBM*(E-PHI(1)+(ZCHARG-LL1/DR)/DR)
      PSIOUT(0)=0
      PSIOUT(1)=1.0E-10
      DO 49 IR=2,NR
            R=DR*IR                 !Numerov algorithm
            K2P=DRHBM*(E-PHI(IR)+(ZCHARG-LL1/R)/R)
            PSIOUT(IR)=(PSIOUT(IR-1)*(2-10*K2Z)-
     +                 PSIOUT(IR-2)*(1+K2M))/(1+K2P)
            K2M=K2Z                 !roll values
            K2Z=K2P
49    CONTINUE
C 
      K2P=0                         !integrate inward homogeneous soln
      R=(NR-1)*DR
      K2Z=DRHBM*(E-PHI(NR-1)+(ZCHARG-LL1/R)/R)
      PSIIN(NR)=0
      PSIIN(NR-1)=1.0E-10
      DO 59 IR=NR-2,1,-1
            R=DR*IR                 !Numerov algorithm
            K2M=DRHBM*(E-PHI(IR)+(ZCHARG-LL1/R)/R)
            PSIIN(IR)=(PSIIN(IR+1)*(2-10*K2Z)-
     +                PSIIN(IR+2)*(1+K2P))/(1+K2M)
            K2P=K2Z                 !roll values
            K2Z=K2M
59    CONTINUE
C 
      NR2=NR/2                       !Wronskian at middle of mesh
      WRON=(PSIIN(NR2+1)-PSIIN(NR2-1))/(2*DR)*PSIOUT(NR2)
      WRON=WRON-(PSIOUT(NR2+1)-PSIOUT(NR2-1))/(2*DR)*PSIIN(NR2)
C 
      SUM=0                           !outward integral in Green's soln
      DO 69 IR=1,NR
            TERM=-PSIOUT(IR)*FOCK(IR,STATE)/2
            SUM=SUM+TERM
            PSTOR(IR,STATE)=PSIIN(IR)*SUM*DR
            SUM=SUM+TERM
69    CONTINUE
C                     
      SUM=0                           !inward integral in Green's soln
      DO 79 IR=NR,1,-1
            TERM=-PSIIN(IR)*FOCK(IR,STATE)/2
            SUM=SUM+TERM
            PSTOR(IR,STATE)=(PSTOR(IR,STATE)+PSIOUT(IR)*SUM*DR)/WRON
            SUM=SUM+TERM
79    CONTINUE
C
      IF (STATE .EQ. 2) THEN          !keep 1s and 2s states orthogonal
             SUM=0
             DO 89 IR=1,NR
                   SUM=SUM+PSTOR(IR,1)*PSTOR(IR,2)
89           CONTINUE
             SUM=SUM*DR
             DO 99 IR=1,NR
                   PSTOR(IR,2)=PSTOR(IR,2)-SUM*PSTOR(IR,1)
99           CONTINUE
      END IF
C
      NORM=0                          !normalize the soln
      DO 18 IR=1,NR
         NORM=NORM+PSTOR(IR,STATE)**2
18    CONTINUE
      NORM=1/SQRT(NORM*DR)
      DO 28 IR=1,NR
         PSTOR(IR,STATE)=PSTOR(IR,STATE)*NORM
28    CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SQR3J(L1,L2,LAM,THREEJ)
C subroutine to calculate square of the 3-j coefficient appearing
C                in the exchange energy
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P3'
C Passed variables:
      INTEGER L1,L2,LAM          !angular momentum (input)
      REAL THREEJ                !three j coefficient (output)
C Local variables:
      REAL DELTA                 !intermediate term for calc of 3J
      INTEGER IMAX,P             !useful combinations of ang mom
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMAX=L1+L2+LAM+1           !useful combinations
      P=(L1+L2+LAM)/2
C
      DELTA=FACTRL(L1+L2-LAM)*FACTRL(-L1+L2+LAM)   !temp values
      DELTA=DELTA*FACTRL(L1-L2+LAM)/FACTRL(IMAX)
      THREEJ=DELTA*(FACTRL(P)**2)                  !calculate 3J squared
      THREEJ=THREEJ/(FACTRL(P-L1)**2)
      THREEJ=THREEJ/(FACTRL(P-L2)**2)
      THREEJ=THREEJ/(FACTRL(P-LAM)**2)
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
      INCLUDE 'PARAM.P3'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
      INTEGER I                    !index for factorial loop
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL SETUP                   !get environment parameters
C 
C     display header screen     
      DESCRP(1)= 'PROJECT 3'
      DESCRP(2)= 'Hartree-Fock solutions of small atomic systems'
      DESCRP(3)= 'in the filling approximation'  
      NHEAD=3
C 
C     text output description
      DESCRP(4)='kinetic, potential '
     +   //'(electron-nucleus, electron-electron,'
      DESCRP(5)= 'exchange and total) and total energy for each state' 
      NTEXT=2
C 
C     graphics output description
      DESCRP(6)= 'electron probability density for each state'        
      NGRAPH=1
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C 
C     calculate constants                
      !atomic constants
      HBARM=7.6359             !hbar*2/(mass)
      CHARGE=14.409            !charge of the electron
      ABOHR=HBARM/CHARGE       !Bohr radius
C 
      NSTATE=MAXSTT            !descriptions of states
      ANGMOM(1)=0
      ANGMOM(2)=0
      ANGMOM(3)=1
      ID(1)=' 1S  '
      ID(2)=' 2S  '
      ID(3)=' 2P  '
      ID(4)='TOTAL'
C      
      FACTRL(0)=1              !factorials
      DO 10 I=1,10
         FACTRL(I)=I*FACTRL(I-1)
10    CONTINUE      
C 
C     setup menu arrays
      CALL MENU                 !setup constant part of menu
C                
      MTYPE(13)=FLOAT
      MPRMPT(13)= 'Enter nuclear charge' 
      MTAG(13)= 'Nuclear charge'
      MLOLIM(13)=1.
      MHILIM(13)=20.
      MREALS(13)=6.
C 
      MTYPE(14)=NUM
      MPRMPT(14)= 'Enter number of electrons in the 1s state' 
      MTAG(14)= 'Occupation of 1s state'       
      MLOLIM(14)=0.
      MHILIM(14)=2.
      MINTS(14)=2.
C              
      MTYPE(15)=NUM
      MPRMPT(15)= 'Enter number of electrons in the 2s state' 
      MTAG(15)= 'Occupation of 2s state'       
      MLOLIM(15)=0.
      MHILIM(15)=2.
      MINTS(15)=2.
C              
      MTYPE(16)=NUM
      MPRMPT(16)= 'Enter number of electrons in the 2p state' 
      MTAG(16)= 'Occupation of 2p state'       
      MLOLIM(16)=0.
      MHILIM(16)=6.
      MINTS(16)=2.
C              
      MTYPE(17)=SKIP
      MREALS(17)=35.
C               
      MTYPE(38)=FLOAT
      MPRMPT(38)= 'Enter radial step size (Angstroms)'
      MTAG(38)= 'Radial step size (Angstroms)'
      MLOLIM(38)=.0001
      MHILIM(38)=.5
      MREALS(38)=.01
C             
      MTYPE(39)=FLOAT
      MPRMPT(39)= 'Enter outer radius of the lattice (Angstroms)' 
      MTAG(39)= 'Outer radius of lattice (Angstroms)'
      MLOLIM(39)=.01
      MHILIM(39)=10.
      MREALS(39)=2.
C 
      MTYPE(40)=NUM
      MPRMPT(40)= 'Enter number of iterations' 
      MTAG(40)= 'Number of iterations'         
      MLOLIM(40)=1.
      MHILIM(40)=50.
      MINTS(40)=4.
C                   
      MTYPE(41)=SKIP
      MREALS(41)=60.
C 
      MSTRNG(MINTS(75))= 'proj3.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C 
      MSTRNG(MINTS(86))= 'proj3.grf'
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
      INCLUDE 'PARAM.P3'
      INCLUDE 'MAP.P3'
C Local variables:
      INTEGER I                     !index of current state
C Function:                      
      LOGICAL LOGCVT                !converts 1 to TRUE, others to FALSE
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
      Z=MREALS(IZ)
      NOCC(1)=MINTS(IONE)
      NOCC(2)=MINTS(ITWO)
      NOCC(3)=MINTS(ITHREE)
      DR=MREALS(IDR)
      RMAX=MREALS(IRMAX)
      NITER=MINTS(INITER)
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
C     calculate derivative parameters:
      ZCHARG=Z*CHARGE      !nuclear charge
      NR=INT(RMAX/DR)      !number of radial steps
      NOCC(NSTATE+1)=0     !total number of electrons
      DO 10 I=1,NSTATE
         NOCC(NSTATE+1)=NOCC(NSTATE+1)+NOCC(I)
10    CONTINUE      
C 
      CALL PCHECK        !check input
      CALL CLEAR
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE PCHECK
C ensure that the number of radial steps is not greater than the
C size of the arrays
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global parameters:
      INCLUDE 'PARAM.P3'
      INCLUDE 'MENU.ALL'
      INCLUDE 'IO.ALL'
      INCLUDE 'MAP.P3'
C Functions:
      REAL GETFLT                       !get float from screen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
10    IF ( NR .GT. MAXSTP) THEN
          WRITE (OUNIT,15) REAL(NR),MAXSTP
          MLOLIM(IDR)=RMAX/MAXSTP   !revise lower limit
          MREALS(IDR)=MLOLIM(IDR)
          MREALS(IDR) =  GETFLT(MREALS(IDR),MLOLIM(IDR),
     +                   MHILIM(IDR),'Enter a larger step')
          DR=MREALS(IDR)                                
          NR=INT(RMAX/DR)      
      GOTO 10
      END IF
C 
15    FORMAT (' Total number of radial steps (=',1PE9.2, 
     +          ') is larger than maxstp (=',i5,')')
C 
      RETURN
      END 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRMOUT(MUNIT,ZSTAR,NLINES)
C outputs text results to the specified unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables
      INCLUDE 'IO.ALL'                 
      INCLUDE 'PARAM.P3'
C Passed variables
      INTEGER MUNIT             !current unit number (input)
      REAL ZSTAR                !optimal charge (input)
      INTEGER IS                !indexes states (input)
      INTEGER NLINES            !number of lines written to screen (I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) CALL CLEAR
C
      WRITE (MUNIT,19)
      WRITE (MUNIT,21)
      WRITE (MUNIT,27)  Z,ZSTAR
      WRITE (MUNIT,23)  RMAX,DR
      WRITE (MUNIT,30) (NOCC(IS),IS=1,NSTATE)
      WRITE (MUNIT,33)
C
      IF (MUNIT .EQ. GUNIT) THEN   !special heading for graphics file
          WRITE (MUNIT,19)
          WRITE (MUNIT,35)
          WRITE (MUNIT,40)
          WRITE (MUNIT,19)
      ELSE                                                         
          WRITE (MUNIT,19)
      END IF
C
      NLINES=7
C
19    FORMAT  (' ')
21    FORMAT  (' Output from project 3:',
     +         ' Hartree-Fock solutions for small atomic systems')
27    FORMAT  (' nuclear charge=',F6.3,5X,' zstar =',F6.3)
23    FORMAT  (' Rmax (Angstroms)=',F6.3,5X,' radial step (Angstroms)=',
     +         1PE12.5)                                                 
30    FORMAT  (' Occupation of the states are:',3(2X,I2))
33    FORMAT  (' All energies are in eV')
35    FORMAT  (4X,'   radius  ',4X,'1s density',5X,'2s density',
     +         5X,'2p density',4X,'total density')
40    FORMAT  (4X,'-----------',4X,'----------',5X,'----------',
     +         5X,'----------',4X,'-------------')              
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,ITER,E,NLINES)
C writes energies for one iteration to the requested unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'   
      INCLUDE 'PARAM.P3'
C Passed variables:                   
      REAL E(MAXSTT+1,8)        !all energies of all states (input)
      INTEGER MUNIT             !current unit (input)
      INTEGER ITER              !current iteration (input)
      INTEGER NLINES            !number of lines written to screen (I/O)
C Local variables:
      INTEGER IS                !state index
      INTEGER DELINE            !number of lines written per call
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DELINE=3+NSTATE           !number of lines written per call
C     clear the screen at a convenient place
      IF ((NLINES+DELINE  .GT. TRMLIN-6) .AND. (MUNIT .EQ. OUNIT)) THEN
          CALL PAUSE('to continue...',1)
          CALL CLEAR
          NLINES=0
      END IF
C             
      WRITE (MUNIT,20) ITER 
      WRITE (MUNIT,30) 
C
      DO 50 IS=1,NSTATE+1
        IF (NOCC(IS) .NE. 0) THEN   !occupied states only
          WRITE (MUNIT,35) ID(IS),NOCC(IS),E(IS,IKTOT),E(IS,IVEN),
     +    E(IS,IVEE),E(IS,IVEX),E(IS,IVTOT),E(IS,ITOT)
        END IF
50    CONTINUE
C
20    FORMAT  (26X,'------- Iteration ',I2,' -------')
30    FORMAT  (2X,'State',2X,'Nocc',4X,'Ktot',7X,'Ven',7X,'Vee',7X,
     +         'Vex',6X,'Vtot',6X,'Etot')
35    FORMAT  (2X,A5,3X,I2,1X,6(1X,F9.3))
C
      IF (MUNIT .EQ. OUNIT) NLINES=NLINES+DELINE
C
      RETURN  
      END                             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFOUT(DEVICE,PSTOR,RHO,E)
C graph the densities of each state and the total density
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables                                                   
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P3'
      INCLUDE 'GRFDAT.ALL'
C Input variables:              
      INTEGER DEVICE               !which device is being used?
      REAL PSTOR(0:MAXSTP,MAXSTT)  !radial wave function
      REAL RHO(0:MAXSTP)           !density
      REAL E(MAXSTT+1,8)           !all energies of all states
C Local variables
      REAL X,Y                     !radius and density
      REAL R                       !radius for GUNIT
      REAL DEN                     !density for gunit
      INTEGER IS                   !array of occupied states
      INTEGER OCC                  !number of occupied states
      INTEGER IR                   !indexes radius
      INTEGER STATE                !current state
      CHARACTER*9 CZ,CE,COCC       !Z,E,NOCC  as a character strings
      INTEGER ZLEN,ELEN,OLEN       !character length
      DIMENSION X(0:MAXSTP),Y(0:MAXSTP)
      DIMENSION IS(MAXSTT+1),DEN(MAXSTT+1)
      INTEGER SCREEN               !send to terminal
      INTEGER PAPER                !make a hardcopy
      INTEGER FILE                 !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     how many occupied states are there, and which are they?
      OCC=0                                      
      DO 60 STATE=1,NSTATE
         IF (NOCC(STATE) .NE. 0) THEN
             OCC=OCC+1
             IS(OCC)=STATE
         END IF
60    CONTINUE
      IF (OCC .GT. 1) THEN                   !include total density
         OCC=OCC+1                           ! if more than one sp state
         IS(OCC)=NSTATE+1
      END IF
C
C     define parameters which are the same for all plots
      IF (DEVICE .NE. FILE) THEN
          NPLOT=OCC                          !how many plots?
C
          XMIN=0.                            !x value is radius
          XMAX=RMAX
          Y0VAL=0. 
C                                                 
          YMAX=0.                            !y values are density
          DO 50 IR=1,NR
             IF (RHO(IR) .GT. YMAX) YMAX=RHO(IR)
50        CONTINUE
          YMIN=0.
          X0VAL=0.
C
          NPOINT=NR+1
          ILINE=1                      !line and symbol styles
          ISYM=1
          IFREQ=1
          NXTICK=5
          NYTICK=5
C                                                             
          CALL CONVRT(Z,CZ,ZLEN)       !titles and labels
          CALL CONVRT(E(OCC,8),CE,ELEN)
          TITLE = 'Atomic Hartree Fock, z='//CZ(1:ZLEN)//
     +             ', Total Energy='//CE
          LABEL(1)= 'radius (Angstroms)'
C
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
      END IF
C                                      
      IF (DEVICE .NE. FILE) THEN                           
C       for graphics pack, do one plot for each occupied state
        DO 20 IPLOT=1,NPLOT
C      
          STATE=IS(IPLOT)
          CALL CONVRT(E(IPLOT,8),CE,ELEN)
          CALL ICNVRT(NOCC(STATE),COCC,OLEN)
          INFO='Energy = '//CE(1:ELEN)//', occp='//COCC
          LABEL(2)=ID(STATE)//' probability density'
          CALL LNLNAX
C
          DO 30 IR=0,NR
            X(IR)=IR*DR
            IF (STATE .EQ. NSTATE+1) THEN
                Y(IR)=RHO(IR)
            ELSE
                Y(IR)=NOCC(STATE)*PSTOR(IR,STATE)**2
            END IF
30        CONTINUE
          CALL XYPLOT(X,Y)
20      CONTINUE                                                      
C
      ELSE
C       for gunit, write one only line for each radius 
        DO 120 IR=0,NR
          R=IR*DR
          DO 110 STATE=1,NSTATE+1
             IF (STATE .EQ. NSTATE+1) THEN
                DEN(STATE)=RHO(IR)
             ELSE
                DEN(STATE)=NOCC(STATE)*PSTOR(IR,STATE)**2
             END IF
110       CONTINUE
          WRITE (GUNIT,70) R,(DEN(STATE),STATE=1,NSTATE+1)
120     CONTINUE
C
      END IF
C
C     end graphing session
      IPLOT=NPLOT                               !reset index 
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)  !end graphics package
      IF (DEVICE .EQ. SCREEN) CALL TMODE        !switch to text mode
C
70    FORMAT (5(4X,1PE11.3))
100   FORMAT (/,' Patience, please; output going to a file.')
C
      RETURN
      END
