CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM EXMPL5
C     Example 5: Determining Nuclear Charge Density
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company, Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop/ execute once for each set of param
        CALL PARAM        !get input from screen
        CALL ARCHON       !calculate nuclear charge density
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON
C calculates nuclear charge density from best fit to electron scattering
C data; the density is expanded in a finite Fourier series
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E5'
      INCLUDE 'IO.ALL'
C Local variables:
      INTEGER ITER                      !number of iterations
      INTEGER INCRS,CONT                !increase nsine? continue iter?
      INTEGER NSINE                     !number of coefficients
      REAL CZERO(CMAX)                  !Fourier coefficients
      REAL SIGT(DATMAX)                 !total cross section
      COMPLEX FTOTAL(DATMAX)            !total scattering amplitude
      REAL QBASIS                       !max num momentum transfer
      DOUBLE PRECISION A(CMAX+1,CMAX+1) !matrix for inverting
      REAL CHI(NLEG)                    !profile functions
      REAL RHO(NGRF)                    !density for graphing
      REAL DRHO(NGRF)                   !error in density 
      INTEGER SCREEN                    !send to terminal
      INTEGER PAPER                     !make a hardcopy
      INTEGER FILE                      !send to a file
C Function:
      INTEGER GETINT,YESNO              !get input from screen
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     read in nuclear data and calculate derivative parameters
      CALL NDATA
C     set up arrays that do not depend on Cn's
      CALL TABLES(A)
C     prompt for value of N to begin with
      NSINE=GETINT(NBASIS,1,NBASIS,
     +         'Enter number of sine functions to begin with')
      CALL CLEAR
C
C     output summary of parameters
      IF (TTERM) CALL PRMOUT(OUNIT)
      IF (TFILE) CALL PRMOUT(TUNIT)
      IF (GFILE) CALL PRMOUT(GUNIT)
C
C     obtain initial Cn's from Fermi distribution
      ITER=0
      QBASIS=NSINE*PI/RMAX                 !largest Q described by basis
      CALL FERMI(NSINE,CZERO)
C     calculate initial fit and display
      CALL FIT(NSINE,CZERO,ITER,SIGT,FTOTAL,QBASIS,CHI,A,RHO,DRHO)
C
100   CONTINUE                             !loop over iterations
          ITER=ITER+1   
C         find changes in CZERO's to minimize the error 
          CALL MINIMZ(CZERO,NSINE,SIGT,FTOTAL,QBASIS,CHI,A)
C         calculate quality of fit, errors; display results
          CALL FIT(NSINE,CZERO,ITER,SIGT,FTOTAL,QBASIS,CHI,A,RHO,DRHO)
C
C         various options for continuing    
          IF (NSINE .LT. NBASIS) THEN
            INCRS=
     +      YESNO(0,'Do you want to increment number of sines by one?')
            IF (INCRS .EQ. 1) THEN
               NSINE=NSINE+1
               CZERO(NSINE)=0.      
               CONT=1
               QBASIS=NSINE*PI/RMAX       !largest Q described by basis
            ELSE 
               CONT=YESNO(1,'Do you want to continue iterating?')
            END IF
          ELSE
            CONT=YESNO(1,'Do you want to continue iterating?')
          END IF
C
       IF (CONT .EQ. 1) GOTO 100
C
C      write out only the final version to a file
       IF (GHRDCP) CALL GRFOUT(PAPER,RHO,SIGT,NSINE,DRHO)
       IF (GFILE) CALL GRFOUT(FILE,RHO,SIGT,NSINE,DRHO)
C
       RETURN
       END         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NDATA
C read in nuclear data
C define variables that are nucleus-dependent
C take recoil out of data; calculate effective momentum transfer
C prompt for NBASIS and RMAX
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E5'
      INCLUDE 'IO.ALL'
C Local variables:
      INTEGER I,IR                      !index of input data, radius
      REAL QLAB,RECOIL                  !variables to adjust input data
      REAL MTARGT                       !mass of the target in MeV
      INTEGER NZERO                     !suggestion for NBASIS
C Function:
      INTEGER GETINT                    !get integer input from terminal
      REAL GETFLT                       !get real input from terminal
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     open data files
      IF (NUCL .EQ. CA) THEN
          OPEN(UNIT=DUNIT,FILE='CA.DAT',STATUS='OLD')
      ELSE IF (NUCL .EQ. NI) THEN
          OPEN(UNIT=DUNIT,FILE='NI.DAT',STATUS='OLD')
      ELSE IF (NUCL .EQ. PB) THEN
          OPEN(UNIT=DUNIT,FILE='PB.DAT',STATUS='OLD')
      END IF
C
C     read in data description
      READ (DUNIT,5) TARGET
      READ (DUNIT,*) ZTARGT,ATARGT,EBEAM,NPTS
5     FORMAT (2X,A10)
C             
C     calculate target parameters
      RZERO=1.07*ATARGT**(1./3.)
      MTARGT=940*ATARGT
      ZA=ZTARGT*ALPHA
      KBEAM=EBEAM/HBARC
      VC1=1.+4./3.*ZA*HBARC/(EBEAM*RZERO)
C
C     read in cross sections; close file
      READ (DUNIT,20) (THETA(I),SIGE(I),DSIGE(I),
     +              THETA(I+1),SIGE(I+1),DSIGE(I+1),I=1,NPTS,2)
20    FORMAT (2(1X,F7.3,1X,E10.3,1X,E10.3))
      CLOSE (UNIT=DUNIT)
C
C     do some preliminary adjustments to the data
      QMAX=0.
      DO 30 I=1,NPTS
C        keep the error small
         IF (DSIGE(I) .GT. SIGE(I)) DSIGE(I)=.8*SIGE(I)
C        angle in radians
         THETA(I)=THETA(I)*PI/180.
C        correction to momentum transfer
         QLAB=2*KBEAM*SIN(THETA(I)/2.)
         QEFF(I)=QLAB*VC1
         IF (QEFF(I) .GT. QMAX) QMAX=QEFF(I)
C        take out recoil         
         RECOIL=1.+2.*EBEAM*SIN(THETA(I)/2)**2/MTARGT
         SIGE(I)=SIGE(I)*RECOIL
         DSIGE(I)=DSIGE(I)*RECOIL
30    CONTINUE
C
C     prompt for input of RMAX and NBASIS
      WRITE (OUNIT,40) TARGET
      WRITE (OUNIT,50) EBEAM 
      WRITE (OUNIT,60) QMAX
      WRITE (OUNIT,70) RZERO
      WRITE (OUNIT,*) '  '
C
      RMAX=GETFLT(RZERO+4.,1.,25.,
     +       ' Enter the boundary radius to use (in fermis)')
      NZERO=INT(QMAX*RMAX/PI)
      WRITE (OUNIT,80) NZERO
      NBASIS=GETINT(NZERO,1,15,' Enter NBASIS')
C
C     fill array RGRF that is used for charge density with radial values
      DRGRF=RMAX/NGRF
      DO 90 IR=1,NGRF
         RGRF(IR)=IR*DRGRF
90    CONTINUE                                                          
C
40    FORMAT (' For the target ',A10)
50    FORMAT (' the data are at a beam energy of ',F7.3,' MeV')
60    FORMAT (' the maximum momentum transfer covered is ',
     +        F7.3,' fm**-1')
70    FORMAT (' the nuclear radius is about',F7.3,' fm')
80    FORMAT (' to keep QMAX less than Qexperimental, NBASIS'
     +  ' must be less than ',I2)
C
      RETURN
      END                      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TABLES(A)
C calculates CHI(N), sets up table of B*J_0(QR), calculates FOUTER
C and zeroes the matrix A
C (none of these variables depend on Cn's, and therefore need 
C to be done only once per calculation)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E5'
      INCLUDE 'IO.ALL'
C Output variables:
      DOUBLE PRECISION A(CMAX+1,CMAX+1)!matrix for inverting
C Local variables:
      REAL B                           !impact parameter/RMAX
      INTEGER N,ILEG,JLEG,I,J          !indices
      REAL Z,X,Y                       !integration variables
      REAL SUM                         !sum for integration
      REAL PHI,BESSJ0,BESSJ1           !real functions
      REAL J0,J1                       !temp values for Fouter
      COMPLEX MU                       !arguments of the Lommel function
      REAL NU
      COMPLEX FACTOR                   !common factor in FOUTER
      COMPLEX LOMMEL                   !Lommel function * X**(1-MU)
      COMPLEX TEMP1,TEMP2              !temp storage for Lommel function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     function for CHIN 
      PHI(Y)=ALOG((1.+SQRT(1.-Y**2))/Y)-SQRT(1.-Y**2)
C     Lommel function * X**(1-MU)
      LOMMEL(X,MU,NU)=1.- ((MU-1.)**2-NU**2)/X**2 +
     +                ((MU-1)**2-NU**2)*((MU-3.)**2-NU**2)/X**4
C
C     zero the matrix A
      DO 2 I=1,CMAX+1
         DO 1 J=1,CMAX+1
            A(I,J)=0.
1        CONTINUE
2     CONTINUE
C
C     CHIN = part of profile function independent of Cn's
      WRITE (OUNIT,*) ' Calculating chi(n)...'
      DO 20 N=1,NBASIS
         DO 10 ILEG=1,NLEG
            B=(1+XLEG(ILEG))/2.           !all B's are scaled by RMAX
            SUM=0.
            DO 5 JLEG=1,NLEG
               Z=B+(1.-B)/2.*(1.+XLEG(JLEG))              !Eq. 5.66b
               SUM=SUM+WLEG(JLEG)*Z*SIN(N*PI*Z)*PHI(B/Z)  !Gaussian quad
5           CONTINUE
            CHIN(ILEG,N)=-4*PI*ALPHA*RMAX**2*SUM*(1.-B)/2.
10       CONTINUE                                            
20    CONTINUE
C
C     part of F_inner integrand independent of Cn's
      WRITE (OUNIT,*) ' Calculating Bessel functions * b ...'
      DO 60 ILEG=1,NLEG
         B=(1+XLEG(ILEG))/2.               !all B's are scaled by RMAX
         DO 50 I=1,NPTS                    !Eq. 5.62
            X=QEFF(I)*RMAX*B
            JTABLE(ILEG,I)=BESSJ0(X)*B*WLEG(ILEG)
50       CONTINUE
60    CONTINUE
C
C     calculate F_outer (which only sees a point charge)
      WRITE (OUNIT,*) ' Calculating outer scattering amplitude...'
      WRITE (OUNIT,*) '  '
      FACTOR=-2*SQRTM1*ZA                
      DO 80 I=1,NPTS
         X=QEFF(I)*RMAX            
         J0=BESSJ0(X)
         J1=BESSJ1(X)                      !Eq. 5.63
         TEMP1=LOMMEL(X,FACTOR,-1.)
         TEMP2=LOMMEL(X,1.+FACTOR,0.)
         FOUTER(I)=FACTOR*J0*TEMP1+X*J1*TEMP2-X*J1
         FOUTER(I)=SQRTM1*KBEAM*FOUTER(I)/(QEFF(I)**2)
80    CONTINUE
C                         
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE FERMI(NSINE,CZERO)
C calculate initial CZERO by finding Fourier coefficients 
C for the Fermi Function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'PARAM.E5'
C Input variables:
       INTEGER NSINE                     !number of coefficients
C Output variables
       REAL CZERO(CMAX)                  !Fourier coefficients
C Local variables:
       INTEGER ILEG,N                    !indices
       REAL RADIUS,R                     !radial step,radius
       REAL FERMIF,RRHO                  !Fermi Function
       REAL SUM                          !sum for normalization
       REAL THICK                        !thickness of Fermi function
       DATA THICK / 2.4 /
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       FERMIF(R)=1./(1.+EXP(4.4*(R-RZERO)/THICK))
C
       DO 30 N=1,NSINE                   !zero sums
          CZERO(N)=0.
30     CONTINUE
C      radial integrals to find C's
       DO 50 ILEG=1,NLEG
          RADIUS=RMAX*(1+XLEG(ILEG))/2
          RRHO=RADIUS*FERMIF(RADIUS)     !Fourier transform of Eq. 5.57
          DO 40 N=1,NSINE                !Gaussian quad 
             CZERO(N)=CZERO(N)+RRHO*SIN(N*PI*RADIUS/RMAX)*WLEG(ILEG)
40        CONTINUE
50     CONTINUE
C                       
C      normalize the C's to charge ZTARGT (Eq. 5.58)
       SUM=0.   
       DO 70 N=1,NSINE
            IF (MOD(N,2) .EQ. 1) SUM=SUM+CZERO(N)/N 
            IF (MOD(N,2) .EQ. 0) SUM=SUM-CZERO(N)/N
70     CONTINUE
       SUM=ZTARGT/(4*RMAX**2*SUM)
       DO 80 N=1,NBASIS                
          CZERO(N)=CZERO(N)*SUM
80     CONTINUE
C
       RETURN
       END                     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FIT(NSINE,CZERO,ITER,SIGT,FTOTAL,QBASIS,CHI,A,RHO,DRHO)
C calculates sigma (given CZERO), chi squared, the charge density, 
C and error in charge density; displays the results
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E5'
      INCLUDE 'IO.ALL'
C Input variables:
      REAL CZERO(CMAX)                     !Fourier coefficients
      INTEGER NSINE                        !number of coefficients
      INTEGER ITER                         !number of iterations
C Output variables:      
      REAL SIGT(DATMAX)                    !total cross section
      COMPLEX FTOTAL(DATMAX)               !total scattering amplitude
      REAL QBASIS                          !max num momentum transfer
      REAL RHO(NGRF)                       !charge density 
      REAL DRHO(NGRF)                      !error in density 
      DOUBLE PRECISION A(CMAX+1,CMAX+1)    !matrix for inverting
C Local variables:
      REAL CHI(NLEG)                       !profile functions
      COMPLEX FINNER                       !inner scattering amplitude
      INTEGER ILEG,IQ,N,IR,M               !indices
      REAL CHISQ                           !goodness of fit
      INTEGER NDOF                         !degrees of freedom in fit
      REAL SINES(CMAX)                     !table of sines
      REAL SUM,FAC                         !temp variable for DRHO
      REAL Z                               !total nuclear charge
      INTEGER SCREEN                       !send to terminal
      INTEGER PAPER                        !make a hardcopy
      INTEGER FILE                         !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     calculate total profile function at each b (Eq. 5.66a)
      DO 20 ILEG=1,NLEG                       
         CHI(ILEG)=-ZA*LOG((1+XLEG(ILEG))/2)
         DO 10 N=1,NSINE
            CHI(ILEG)=CHI(ILEG)+CZERO(N)*CHIN(ILEG,N)
10       CONTINUE
20    CONTINUE
C
C     calculate fit sigma, chisquare
      CHISQ=0
      NDOF=0                           
      DO 50 IQ=1,NPTS                          !loop over data points
C        calculate the inner and total scattering amplitude
         FINNER=0.
         DO 40 ILEG=1,NLEG                     !Eq. 5.62
           FINNER=FINNER+JTABLE(ILEG,IQ)*(EXP(2*SQRTM1*CHI(ILEG))-1.)
40       CONTINUE
         FINNER=-SQRTM1*KBEAM*RMAX**2*FINNER/2.
C        scattering amplitude 
         FTOTAL(IQ)=(FINNER+FOUTER(IQ)) 
C
C        sigma in mbarnes/sr, with momentum transfer corr  (Eq. 5.60)
         SIGT(IQ)=10*COS(THETA(IQ)/2)**2*VC1**2*ABS(FTOTAL(IQ))**2
         IF (QEFF(IQ) .LE. QBASIS) THEN 
           CHISQ=CHISQ+((SIGE(IQ)-SIGT(IQ))/DSIGE(IQ))**2  !Eq. 5.42
           NDOF=NDOF+1
         END IF
50    CONTINUE 
      NDOF=NDOF-NSINE                        !degrees of freedom
C
C     calculate the density and its error 
      Z=0.
      DO 70 IR=1,NGRF
         RHO(IR)=0.                            !zero sums
         DRHO(IR)=0.
         DO 60 N=1,NSINE
            SINES(N)=SIN(N*PI*RGRF(IR)/RMAX)   !common term
            RHO(IR)=RHO(IR)+CZERO(N)*SINES(N)  !Eq. 5.57; density*radius
            SUM=0.
            DO 80 M=1,N                        !sums to find err in dens
               FAC=2.
               IF (M .EQ. N) FAC=1.
               SUM=SUM+FAC*A(N,M)*SINES(M)     !Eq. 5.55, sum over M
80          CONTINUE
            DRHO(IR)=DRHO(IR)+SUM*SINES(N)     !Eq. 5.55, sum over N
60       CONTINUE
         DRHO(IR)=SQRT(ABS(DRHO(IR)))/RGRF(IR) !take out radius from RHO
         Z=Z+RHO(IR)*RGRF(IR)                  !integ dens for tot charge 
         RHO(IR)=RHO(IR)/RGRF(IR)              !now RHO is just density
70    CONTINUE
      Z=Z*4*PI*RMAX/NGRF
C
      IF (TTERM) 
     +  CALL TXTOUT(OUNIT,ITER,CHISQ,NDOF,NSINE,QBASIS,CZERO,A,Z)
      IF (TFILE)  
     +  CALL TXTOUT(TUNIT,ITER,CHISQ,NDOF,NSINE,QBASIS,CZERO,A,Z)
      IF (GTERM) THEN
         CALL GRFOUT(SCREEN,RHO,SIGT,NSINE,DRHO)
         CALL CLEAR
      END IF
C                        
      RETURN 
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MINIMZ(CZERO,NSINE,SIGT,FTOTAL,QBASIS,CHI,A)
C finds the corrections to CZERO to minimize the fit 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E5'
C Input variables:
      REAL CZERO(CMAX)                    !Fourier coefficients (I/O)
      INTEGER NSINE                       !number of coefficients
      REAL SIGT(DATMAX)                   !total cross section
      REAL CHI(NLEG)                      !profile functions
      COMPLEX FTOTAL(DATMAX)              !total scattering amplitude
      REAL QBASIS                         !max num momentum transfer
      DOUBLE PRECISION A(CMAX+1,CMAX+1)   !matrix for inverting (I/O)
C Local variables:
      REAL BVEC(CMAX+1)                   !lhs of linear equation
      INTEGER NDIM                        !dimension of matrix          
      COMPLEX DFDC                        !diff of F w.r.t. CZERO
      REAL W(CMAX)                        !diff of SIGT w.r.t. CZERO
      INTEGER IQ,N,ILEG,M                 !indices
      REAL SUM                            !temp value for matrix mulplct
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE (OUNIT,*) ' Calculating and inverting A...'
      NDIM=NSINE+1
C
C     zero the vector B and matrix A
      DO 10 N=1,NSINE
         BVEC(N)=0
         DO 20 M=1,N
            A(N,M)=0.
20       CONTINUE
10    CONTINUE
C
C     calculate the matrix A and vector B
      DO 70 IQ=1,NPTS
        IF (QEFF(IQ) .LT. QBASIS) THEN  !leave off Q's not desc by basis
          DO 40 N=1,NSINE               !first find W terms
            DFDC=0
            DO 30 ILEG=1,NLEG
               DFDC=DFDC+               !terms in W (Eq. 5.67b)
     +          JTABLE(ILEG,IQ)*EXP(2*SQRTM1*CHI(ILEG))*CHIN(ILEG,N)
30          CONTINUE
            DFDC=DFDC*KBEAM*RMAX**2
            W(N)=20.*VC1**4*COS(THETA(IQ)/2)**2   !Eq. 5.67a
     +           *REAL(CONJG(FTOTAL(IQ))*DFDC)    !W in mbarnes/sr
40        CONTINUE
C
          DO 60 N=1,NSINE                         !calculate matrix A
             DO 50 M=1,N                          !Eq. 5.48
                A(N,M)=A(N,M)+DBLE(W(N)*W(M)/DSIGE(IQ)**2)
50           CONTINUE                             !and vector B(Eq 5.48)
             BVEC(N)=BVEC(N)+(SIGE(IQ)-SIGT(IQ))*W(N)/DSIGE(IQ)**2
60        CONTINUE
         END IF
70     CONTINUE
C
C      include dZ/dCn terms, and fill out A using symmetry
       DO 100 N=1,NSINE                             !Eq 5.49a, 5.58
          IF (MOD(N,2) .EQ. 1) A(N,NDIM)=-A(1,1)/N  !scale by A(1,1)
          IF (MOD(N,2) .EQ. 0) A(N,NDIM)=A(1,1)/N
          A(NDIM,N)=A(N,NDIM)
          DO 90 M=1,N          !fill out A matrix using symmetry
             A(M,N)=A(N,M)
90        CONTINUE
100    CONTINUE
       A(NSINE+1,NSINE+1)=0.
C
       CALL MATINV(A,NDIM)     !invert the matrix
C
       DO 120 N=1,NSINE        !multiply A**-1 * B to solve for 
          SUM=0.               !corrections to Czero
          DO 110 M=1,NSINE
             SUM=SUM+A(N,M)*BVEC(M)
110       CONTINUE
          CZERO(N)=CZERO(N)+SUM
120    CONTINUE
C
       RETURN
       END    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE MATINV(A,NDIM)
C inverts the matrix A of dimension NDIM using Gauss-Jordan elimination;
C A is replaced by its inverse
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E5'
C Input/output variables:
      DOUBLE PRECISION A(CMAX+1,CMAX+1) !matrix for inverting
      INTEGER NDIM                 !dimension of matrix          
C Local variables:                 
      DOUBLE PRECISION P(CMAX+1),Q(CMAX+1) !temporary storage
      LOGICAL ZEROED(CMAX+1)       !keeps track of who's zeroed
      INTEGER I,J,K                !indices
      INTEGER KBIG                 !index of row being zeroed
      DOUBLE PRECISION PIVOT       !largest diag element in A
      REAL BIG                     !largest diag elem of non zeroed rows
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     which rows are zeroed? none so far
      DO 10 I=1,NDIM
         ZEROED(I)=.FALSE.
10    CONTINUE
C 
C     loop over all rows
      DO 100 I=1,NDIM
         KBIG=0                    !search for largest diag element
         BIG=0                     !in rows not yet zeroed
         DO 20 J=1,NDIM
           IF (.NOT. ZEROED(J)) THEN
               IF (ABS(A(J,J)) .GT. BIG) THEN
                  BIG=ABS(A(J,J))
                  KBIG=J
               END IF
           END IF
20       CONTINUE
C
C        store the largest diagonal element
         IF (I .EQ. 1) PIVOT=A(KBIG,KBIG)
C
C        if all diagonals are zero, then the matrix is singular
         IF (KBIG .EQ. 0) THEN
            WRITE (OUNIT,*) ' Matrix is singular'
            RETURN
         END IF
C        matrix is ill conditioned if size of elements varies greatly
         IF (ABS(A(KBIG,KBIG)/PIVOT) .LT. 1.E-14) THEN
            WRITE (OUNIT,*) ' Matrix is ill conditioned'
            RETURN
         END IF
C
C        begin zeroing row KBIG
         ZEROED(KBIG)=.TRUE.
         Q(KBIG)=1/A(KBIG,KBIG)
         P(KBIG)=1.
         A(KBIG,KBIG)=0.
C
C        elements above the diagonal
         IF (KBIG .GT. 1) THEN
            DO 30 J=1,KBIG-1
               P(J)=A(J,KBIG)
               Q(J)=A(J,KBIG)*Q(KBIG)
               IF (.NOT. ZEROED(J)) Q(J)=-Q(J)
               A(J,KBIG)=0.
30          CONTINUE
        END IF
C
C       elements to the right of the diagonal
        IF (KBIG .LT. NDIM) THEN
            DO 40 J=KBIG+1,NDIM
               P(J)=A(KBIG,J)
               Q(J)=-A(KBIG,J)*Q(KBIG)
               IF (ZEROED(J)) P(J)=-P(J)
               A(KBIG,J)=0.
40          CONTINUE
        END IF
C
C      transform all of A
       DO 60 J=1,NDIM
          DO 50 K=J,NDIM
             A(J,K)=A(J,K)+P(J)*Q(K)
50        CONTINUE
60     CONTINUE         
100   CONTINUE
C
C     symmetrize A**-1
      DO 110 J=2,NDIM
         DO 120 K=1,J-1
            A(J,K)=A(K,J)
120      CONTINUE
110   CONTINUE
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
      INCLUDE 'PARAM.E5'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     get environment parameters
      CALL SETUP                
C             
C     display header screen     
      DESCRP(1)= 'EXAMPLE 5'
      DESCRP(2)= 'Determining Nuclear Charge Densities'
      NHEAD=2
C 
C     text output description
      DESCRP(3)= 'Density parameters and quality of fit'
      NTEXT=1
C 
C     graphics output description
      DESCRP(4)= 'Charge density and scattering cross section'
      NGRAPH=1
C                
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C 
C     calculate constants                
      PI=4*ATAN(1.0)
      ALPHA=1./137.036
      HBARC=197.329
      SQRTM1=(0.,1.)
C     get data for Gauss-Legendre integration
      CALL GAUSS
C 
C     setup menu arrays, beginning with constant part
      CALL MENU                      
C                
      MPRMPT(4)='2) (not used)'
C                
      MTYPE(13)=MTITLE
      MPRMPT(13)='Choice of nuclei:'
      MLOLIM(13)=2
      MHILIM(13)=1
C 
      MTYPE(14)=MTITLE
      MPRMPT(14)='1) Calcium 40 '
      MLOLIM(14)=0
      MHILIM(14)=0
C 
      MTYPE(15)=MTITLE
      MPRMPT(15)='2) Nickel 58 '
      MLOLIM(15)=0            
      MHILIM(15)=0
C             
      MTYPE(16)=MTITLE
      MPRMPT(16)='3) Lead 208 '
      MLOLIM(16)=0
      MHILIM(16)=1
C 
      MTYPE(17)=MCHOIC
      MPRMPT(17)='Enter your choice'
      MTAG(17)='18 18 18'
      MLOLIM(17)=1
      MHILIM(17)=3
      MINTS(17)=1
      MREALS(17)=1.
C
      MTYPE(18)=SKIP
      MREALS(18)=35
C 
      MTYPE(36)=SKIP
      MREALS(36)=60.
C
      MSTRNG(MINTS(75))= 'exmpl5.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C 
      MSTRNG(MINTS(86))= 'exmpl5.grf'
C                     
      MTYPE(87)=SKIP
      MREALS(87)=90.
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
      INCLUDE 'PARAM.E5'
C Local variables:
      INTEGER INUCL                !map menu arrays to parameters
      PARAMETER (INUCL  = 17 )
C Function:
      LOGICAL LOGCVT               !converts 1 to true, others to false
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
      NUCL=MINTS(INUCL)
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
      SUBROUTINE GAUSS
C establish Gauss-Legendre weights and abscissae for 20 points
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.E5'
C Local variables:
      INTEGER ILEG                   !index for weights and abscissae
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      XLEG(1)=.993128599185094 
      WLEG(1)=.017614007139152
      XLEG(2)=.963971927277913
      WLEG(2)=.040601429800386
      XLEG(3)=.912234428251325
      WLEG(3)=.062672048334109
      XLEG(4)=.839116971822218 
      WLEG(4)=.083276741576704
      XLEG(5)=.74633190646015
      WLEG(5)=.10193011981724
      XLEG(6)=.636053680726515
      WLEG(6)=.118194531961518
      XLEG(7)=.510867001950827
      WLEG(7)=.131688638449176
      XLEG(8)=.373706088715419
      WLEG(8)=.142096109318382
      XLEG(9)=.227785851141645
      WLEG(9)=.149172986472603
      XLEG(10)=.076526521133497
      WLEG(10)=.152753387130725
C
      DO 10 ILEG=1,NLEG/2        !weights and abscissae are even and odd
         XLEG(21-ILEG)=-XLEG(ILEG)
         WLEG(21-ILEG)=WLEG(ILEG)
10    CONTINUE        
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       REAL FUNCTION BESSJ0(X)
C calculates zeroth order Bessel function at X
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
       REAL X
C Local variables:
	REAL TEMP,TEMP2,Y2,Y
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        Y=X/3.	
        IF(ABS(X) .LT. 3.0) THEN
             Y2=Y*Y	
             TEMP=-.0039444+.00021*Y2
             TEMP=.0444479+Y2*TEMP   
             TEMP=-.3163866+Y2*TEMP   
             TEMP=1.2656208+Y2*TEMP
             TEMP=-2.2499997+Y2*TEMP
             BESSJ0=1+Y2*TEMP
        ELSE
             Y=1/Y
             TEMP=-7.2805E-04+1.4476E-04*Y
             TEMP=1.37237E-03+TEMP*Y 
             TEMP=-9.512E-05+TEMP*Y         
             TEMP=-.0055274+TEMP*Y
             TEMP=-7.7E-07+TEMP*Y           
             TEMP=.79788456+TEMP*Y
             TEMP2=-2.9333E-04+1.3558E-04*Y 
             TEMP2=-5.4125E-04+TEMP2*Y
             TEMP2=2.62573E-03+TEMP2*Y      
             TEMP2=-3.954E-05+TEMP2*Y
             TEMP2=-4.166397E-02+TEMP2*Y
             TEMP2=X-.78539816+TEMP2*Y
             BESSJ0=TEMP*COS(TEMP2)/SQRT(X)
        END IF
        RETURN
        END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       REAL FUNCTION BESSJ1(X)
C calculates first order Bessel function at X
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Passed variables:
       REAL X
C Local variables:
	REAL TEMP,TEMP2,Y2,Y
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        Y=X/3.	
        IF(ABS(X) .LT. 3.0) THEN
             Y2=Y*Y	
             TEMP=-3.1761E-04+1.109E-05*Y2 
             TEMP=4.43319E-03+Y2*TEMP    
             TEMP=-3.954289E-02+Y2*TEMP    
             TEMP=.21093573+Y2*TEMP
             TEMP=-.56249985+Y2*TEMP
             BESSJ1=X*(.5+Y2*TEMP)
        ELSE
           Y=1/Y
             TEMP=1.13653E-03-2.0033E-04*Y
             TEMP=-2.49511E-03+TEMP*Y     
             TEMP=1.7105E-04+TEMP*Y
             TEMP=1.659667E-02+TEMP*Y
             TEMP=1.56E-06+TEMP*Y 
             TEMP=.79788456+TEMP*Y
             TEMP2=7.9824E-04-2.9166E-04*Y
             TEMP2=7.4348E-04+TEMP2*Y
             TEMP2=-6.37879E-03+TEMP2*Y 
             TEMP2=.0000565+TEMP2*Y
             TEMP2=.12499612+TEMP2*Y  
             TEMP2=X-2.35619449+TEMP2*Y
             BESSJ1=TEMP*COS(TEMP2)/SQRT(X)
         END IF
         RETURN
         END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PRMOUT(MUNIT)
C outputs parameter summary to the specified unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E5'
C Passed variables:               
      INTEGER MUNIT            !unit number for output
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (MUNIT .EQ. OUNIT) CALL CLEAR
C 
      WRITE (MUNIT,2)
      WRITE (MUNIT,4)
      WRITE (MUNIT,6) TARGET
      WRITE (MUNIT,8) EBEAM
      WRITE (MUNIT,10) RMAX
      WRITE (MUNIT,12) NBASIS
      WRITE (MUNIT,14) QMAX 
      WRITE (MUNIT,16) NBASIS*PI/RMAX
      WRITE (MUNIT,2)
C
2     FORMAT (' ')
4     FORMAT (' Output from example 5: Nuclear Charge Densities')
6     FORMAT (' For the target ',A10)
8     FORMAT (' the data are at a beam energy of ',F7.3,' MeV')
10    FORMAT (' the radial cutoff for the charge density=',F7.3,' fm')
12    FORMAT (' the number of sines used=', I2)
14    FORMAT (' the maximum experimental momentum transfer='
     +        ,F7.3,' fm**-1')
16    FORMAT (' the maximum numerical momentum transfer=',F7.3,'fm**-1')
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TXTOUT(MUNIT,ITER,CHISQ,NDOF,NSINE,QBASIS,CZERO,A,Z)
C outputs the charge density parameters and goodness of fit to MUNIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E5'
C Input variables:
      INTEGER MUNIT                        !unit to write to
      REAL CZERO(CMAX)                     !Fourier coefficients
      INTEGER NSINE                        !number of coefficients
      INTEGER ITER                         !number of iterations
      REAL CHISQ                           !goodness of fit
      INTEGER NDOF                         !degrees of freedom in fit
      DOUBLE PRECISION A(CMAX+1,CMAX+1)    !matrix we're inverting
      REAL QBASIS                          !max num momentum transfer
      REAL Z                               !total charge
C Local variables:
      INTEGER N                            !index of CZERO
      REAL SIGMA                           !error of CZERO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE (MUNIT,*) ' '
      WRITE (MUNIT,*) ' '
      WRITE (MUNIT,10) ITER
      WRITE (MUNIT,15) CHISQ,NDOF
      WRITE (MUNIT,20) NSINE,Z
      WRITE (MUNIT,25) QBASIS
      WRITE (MUNIT,*) ' '
      WRITE (MUNIT,30)
C
      DO 50 N=1,NSINE
         SIGMA=0.
         IF (A(N,N) .GT. 0) SIGMA=SQRT(A(N,N))
         WRITE (MUNIT,40) N,CZERO(N),SIGMA
50    CONTINUE
C
      IF (MUNIT .EQ. OUNIT) CALL PAUSE('to continue...',1)
C
10    FORMAT (' Iteration ', I4)
15    FORMAT (' Chi**2 = ',1PE15.8,' for ',I4,' degrees of freedom')
20    FORMAT (' Number of sines = ',I2, '  total charge =',F10.4)
25    FORMAT (' Maximum momentum transfer for this basis = ',1PE15.8,
     +        ' fm**-1')
30    FORMAT (' The expansion coefficients of the charge density '
     +        'and their errors are:')
40    FORMAT (' C(',I2,') = ',1PE15.8,' +- ',1PE15.8)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE GRFOUT(DEVICE,RHO,SIGT,NSINE,DRHO)
C outputs differential cross section vs. momentum transfer
C and nuclear charge density (with error bars) vs. radius
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables                                                   
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.E5'
      INCLUDE 'GRFDAT.ALL'
C Input variables:
      INTEGER DEVICE         !which device is being used?
      REAL SIGT(DATMAX)      !total cross section
      REAL RHO(NGRF)         !density for graphing
      REAL DRHO(NGRF)        !error in density 
      INTEGER NSINE          !number of sines   
C Local variables
      INTEGER IDATA,IR       !indexes data, density
      INTEGER EXPMAX,EXPMIN  !min and max exp for diff cross section
      CHARACTER*9 CEBEAM     !EBEAM as a character string 
      INTEGER LENGTH         !length of character strings
      CHARACTER*9 CSINE      !number of sines as a string
      INTEGER SINLEN
      REAL X(2),Y(2)         !arrays for error bars
      INTEGER SCREEN                  !send to terminal
      INTEGER PAPER                   !make a hardcopy
      INTEGER FILE                    !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     calculate parameters for graphing the cross sections
      IF (DEVICE .NE. FILE) THEN
          NPLOT=1                        !how many plots?
          IPLOT=1
C 
          YMAX=0.                        !find limits on data points
          YMIN=SIGT(1)
          DO 20 IDATA=1,NPTS
             IF (SIGT(IDATA) .GT. YMAX) YMAX=SIGT(IDATA)
             IF (SIGE(IDATA) .GT. YMAX) YMAX=SIGE(IDATA)
             IF (SIGT(IDATA) .LT. YMIN) YMIN=SIGT(IDATA)
             IF (SIGE(IDATA) .LT. YMIN) YMIN=SIGE(IDATA)
20        CONTINUE                                             
C         find integer limits on exponent
          EXPMAX=INT(LOG10(YMAX))
          IF (YMAX .GT. 1.) EXPMAX =EXPMAX+1
          EXPMIN=INT(LOG10(YMIN))
          IF (YMIN .LT. 1.) EXPMIN=EXPMIN-1
          YMAX=10.**EXPMAX
          YMIN=10.**EXPMIN
C 
          XMIN=QEFF(1)                 !more limits
          XMAX=QEFF(NPTS)
          Y0VAL=XMIN
          X0VAL=YMIN
C          
          NPOINT=NPTS
C 
          ILINE=1                      !line and symbol styles
          ISYM=1
          IFREQ=1
          NXTICK=5
          NYTICK=EXPMAX-EXPMIN
          IF (NYTICK .GT. 8) THEN      !keep number of ticks small
             IF (MOD(NYTICK,2) .EQ. 0) THEN
                NYTICK=NYTICK/2
             ELSE            
                NYTICK=8
             END IF
          END IF
C
          CALL CONVRT(EBEAM,CEBEAM,LENGTH)            !titles and labels
          INFO = 'Calculated(X) and Experimental(0)'
          LABEL(1)= 'Qeffective (fm**-1)'
          LABEL(2)= 'Differential Cross Section (mBarnes/sr)'
          TITLE = ' Electron Scattering on '//TARGET
     +            //' at '//CEBEAM(1:LENGTH)//' MeV'
C
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
          CALL LGLNAX                          !draw axes
      END IF
C                                                      
C     output cross sections vs. angle
      IF (DEVICE .EQ. FILE) THEN
          WRITE (GUNIT,200)
          WRITE (GUNIT,210)
          WRITE (GUNIT,220) (QEFF(IDATA),SIGT(IDATA),IDATA=1,NPTS)
      ELSE                                     
         CALL XYPLOT(QEFF,SIGE)               !plot experimental data
         NPOINT=2
         IFREQ=0
         DO 80 IDATA=1,NPTS                    !with error bars
            X(1)=QEFF(IDATA)
            X(2)=QEFF(IDATA)
            Y(1)=SIGE(IDATA)+DSIGE(IDATA)
            Y(2)=SIGE(IDATA)-DSIGE(IDATA)
            CALL XYPLOT(X,Y)
80       CONTINUE
         NPOINT=NPTS
         IFREQ=1
         ISYM=4
         CALL XYPLOT (QEFF,SIGT)              !plot calculated sigma
         CALL GPAGE(DEVICE)
      END IF
C
C     calculate parameters for charge density
      IF (DEVICE .NE. FILE) THEN
          YMAX=0.                        !find limits on data points
          YMIN=0.
          DO 120 IR=1,NGRF
             IF (RHO(IR) .GT. YMAX) YMAX=RHO(IR)
120       CONTINUE
          XMIN=0.
          XMAX=RMAX
          Y0VAL=XMIN
          X0VAL=YMIN
C          
          NPOINT=NGRF
C 
          ILINE=1                      !line and symbol styles
          ISYM=4
          IFREQ=1
          NXTICK=5
          NYTICK=5
          IFREQ=0
C
          CALL ICNVRT(NSINE,CSINE,SINLEN)   !titles and labels
          INFO=' '
          LABEL(1)= 'Radius (fermis)'
          LABEL(2)= 'Nuclear charge density (fermi**-3)'
          TITLE = TARGET//' using '//CSINE(1:SINLEN)//' sine functions'
C
          CALL LNLNAX                          !draw axes
      END IF
C                                                      
C     output charge density and its error
      IF (DEVICE .EQ. FILE) THEN
         WRITE (GUNIT,*) '  '
         WRITE (GUNIT,300)
         WRITE (GUNIT,310)
         WRITE (GUNIT,320)(RGRF(IR),RHO(IR),DRHO(IR),IR=1,NGRF)
      ELSE
         CALL XYPLOT (RGRF,RHO)               !charge density
         NPOINT=2
         IFREQ=0
         DO 180 IR=1,NGRF                      !with error bars
            X(1)=RGRF(IR)
            X(2)=X(1) 
            Y(1)=RHO(IR)-DRHO(IR)
            Y(2)=Y(1)+2*DRHO(IR)
            CALL XYPLOT(X,Y)
180      CONTINUE
      END IF
C 
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)  !end graphics package
      IF (DEVICE .EQ. SCREEN) CALL TMODE        !switch to text mode
C
100   FORMAT (/,' Patience, please; output going to a file.')
200   FORMAT (10X,'Qeff (1/fermi)',9X,'Sigma (mBarnes/sr)')
210   FORMAT (10X,'--------------',9X,'------------------')
220   FORMAT (2(10X,1PE15.8))
300   FORMAT (12X,'R (fermi)',8X,'Charge density (fm**-3)',10X,'Error')
310   FORMAT (12X,'---------',8X,'-----------------------',10X,'-----')
320   FORMAT (3(9X,1PE15.8))
C
      RETURN
      END
 
