CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM PROJ5
C     PROJECT 5: Solution of a schematic shell model
C  COMPUTATIONAL PHYSICS (FORTRAN VERSION)
C  by Steven E. Koonin and Dawn C. Meredith
C  Copyright 1989, Addison-Wesley Publishing Company Inc.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL INIT           !display header screen, setup parameters
5     CONTINUE            !main loop/ execute once for each set of param
        CALL PARAM        !get input from screen
        CALL ARCHON       !find eigenvalues and eigenvectors
      GOTO 5
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ARCHON
C calculates the eigenvalues and (optionally) the eigenvectors of a
C schematic shell model for several values of the coupling strength CHI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P5'
C Local variables
      REAL CEVEN(MAXDIM,MAXDIM),CODD(MAXDIM,MAXDIM)  !eigenvectors
      REAL EVEVEN(MXNCHI,MAXDIM),EVODD(MXNCHI,MAXDIM)!eigenvalues
      INTEGER JCHI                !CHI index
      REAL V                      !coupling strength
      INTEGER NLINES              !number of lines printed out
      INTEGER SCREEN              !send to terminal
      INTEGER PAPER               !make a hardcopy
      INTEGER FILE                !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     output summary of parameters
      IF (TTERM) CALL PRMOUT(OUNIT,NLINES)
      IF (TFILE) CALL PRMOUT(TUNIT,NLINES)
      IF (GFILE) CALL PRMOUT(GUNIT,NLINES)
C
      DO 100 JCHI=1,NCHI     !for each value of CHI find the eigenvalues
         V=CHI(JCHI)/(NPART-1)    !and eigenvectors; output text results
         CALL HAMLTN(JCHI,V,CEVEN,CODD,EVEVEN,EVODD)
         CALL TXTOUT(JCHI,V,CEVEN,CODD,EVEVEN,EVODD,NLINES)
100   CONTINUE
      IF (TTERM) CALL PAUSE('to continue...',1)
      IF (TTERM) CALL CLEAR
C
C     graphics output; plot is different if there is only one CHI value
      IF (NCHI .NE. 1) THEN
        IF (GTERM) CALL GRFOUT(SCREEN,EVEVEN,EVODD)
        IF (GFILE) CALL GRFOUT(FILE,EVEVEN,EVODD)
        IF (GHRDCP) CALL GRFOUT(PAPER,EVEVEN,EVODD)
      ELSE
        IF (GTERM) CALL ONEOUT(SCREEN,EVEVEN,EVODD)
        IF (GFILE) CALL ONEOUT(FILE,EVEVEN,EVODD)
        IF (GHRDCP) CALL ONEOUT(PAPER,EVEVEN,EVODD)
      END IF
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE HAMLTN(JCHI,V,CEVEN,CODD,EVEVEN,EVODD)
C subroutine to set up and diagonalize the Hamiltonian for one
C value of the coupling, V
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P5'
C Input Variables
      INTEGER JCHI                    !index of the coupling
      REAL V                          !coupling strength
C Output variables:
      REAL CEVEN(MAXDIM,MAXDIM),CODD(MAXDIM,MAXDIM)  !eigenvectors
      REAL EVEVEN(MXNCHI,MAXDIM),EVODD(MXNCHI,MAXDIM)!eigenvalues
C Local Variables
      INTEGER IDIM,ICOMP              !index on matrix elements
      INTEGER NDIM                    !size of matrix
      REAL DIAG(MAXDIM),LDIAG(MAXDIM) !diagonal,off diagonal matrix elem
      INTEGER M                       !eigenvalue of Jz
      REAL TEMP                       !temp storage for LDIAG
      INTEGER NFIND                    !number of eigenvalues to find
      REAL EVAL(MAXDIM),TMPVEC(MAXDIM) !temp storage for values and vect
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     find evectors and eigenvalues for even-m states
      NDIM=JJ+1                            !number of even-m states
      DO 100 IDIM=1,NDIM                   !set up the matrix elements
         M=-JJ+2*(IDIM-1)
         DIAG(IDIM)=M                      !diagonal matrix element
         IF (IDIM.NE.1) THEN               !lower off-diag matrix elem
            TEMP=JJ1-M*(M-1)
            TEMP=TEMP*(JJ1-(M-1)*(M-2))
            LDIAG(IDIM)=-V/2*SQRT(ABS(TEMP))
         ENDIF
100   CONTINUE
C
      IF (MOD(JJ,2) .EQ. 0)THEN            !find only neg eigenvalues
         NFIND=1+JJ/2
      ELSE
         NFIND=(JJ+1)/2
      ENDIF
      CALL DGNLZ(DIAG,LDIAG,NDIM,NFIND,EVAL)  !find the eigenvalues
C
      DO 200 IDIM=1,NFIND
         EVEVEN(JCHI,IDIM)=EVAL(IDIM)      !save eigenvalues
         IF (VECTRS) THEN                  !sometimes find eigenvectors
            CALL INVITR(DIAG,LDIAG,EVAL(IDIM),NDIM,TMPVEC)
            DO 300 ICOMP=1,NDIM
              CEVEN(ICOMP,IDIM)=TMPVEC(ICOMP)
300         CONTINUE
         ENDIF
200   CONTINUE
C
C     find evectors and eigenvalues for odd-m states
      NDIM=JJ                              !number of odd-m states
      DO 400 IDIM=1,NDIM                   !set up the matrix elements
         M=-JJ-1+2*IDIM
         DIAG(IDIM)=M                      !diagonal matrix element
         IF (IDIM .NE. 1) THEN             !lower off-diag matrix elem
            TEMP=JJ1-M*(M-1)
            TEMP=TEMP*(JJ1-(M-1)*(M-2))
            LDIAG(IDIM)=-V/2*SQRT(ABS(TEMP))
         ENDIF
400   CONTINUE
C
      IF (MOD(JJ,2) .EQ. 0) THEN           !find only neg eigenvalues
         NFIND=JJ/2
      ELSE
         NFIND=(JJ+1)/2
      ENDIF
      CALL DGNLZ(DIAG,LDIAG,NDIM,NFIND,EVAL)  !find the eigenvalues
C
      DO 500 IDIM=1,NFIND
         EVODD(JCHI,IDIM)=EVAL(IDIM)       !save eigenvalues
         IF (VECTRS) THEN                  !sometimes find eigenvectors
            CALL INVITR(DIAG,LDIAG,EVAL(IDIM),NDIM,TMPVEC)
            DO 600 ICOMP=1,NDIM
              CODD(ICOMP,IDIM)=TMPVEC(ICOMP)
600         CONTINUE
         ENDIF
500   CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DGNLZ(DIAG,LDIAG,NDIM,NFIND,EVAL)  
C finds NFIND eigenvalues (EVAL) of a NDIM x NDIM symmetric tridiagonal 
C martix with DIAG and LDIAG matrix elements
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global parameters:
      INCLUDE 'PARAM.P5'
C Input Variables
      REAL DIAG(MAXDIM),LDIAG(MAXDIM)   !diag,off-diag matrix elem
      INTEGER NDIM,NFIND    !matrix dim, number of values to find
C Output Variables:
      REAL EVAL(MAXDIM)     !eigenvalues
C Local Variables
      REAL UBOUND,LBOUND    !bounds on eigenvalues
      REAL RAD,GER          !temp storage for computing bounds
      REAL SPACNG           !estimate of eigenvalue spacing
      REAL DLAM             !step in eigenvalue search
      INTEGER IDIM          !index
      REAL LAMBDA           !current guess for the eigenvalue
      INTEGER COUNT,L       !how many egnvls are .lt. LAMBDA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
C     find Gerschgorin bounds on eigenvalues (Eq. 5.11)
      LBOUND=DIAG(1)-ABS(LDIAG(2))
      UBOUND=DIAG(1)+ABS(LDIAG(2))
      DO 100 IDIM=2,NDIM-1
         RAD=ABS(LDIAG(IDIM+1))+ABS(LDIAG(IDIM))
         GER=DIAG(IDIM)-RAD
         LBOUND=AMIN1(GER,LBOUND)
         GER=DIAG(IDIM)+RAD
         UBOUND=AMAX1(GER,UBOUND)
100   CONTINUE
      GER=DIAG(NDIM)-ABS(LDIAG(NDIM))
      LBOUND=AMIN1(LBOUND,GER)
      GER=DIAG(NDIM)+ABS(LDIAG(NDIM))
      UBOUND=AMAX1(UBOUND,GER)
C
      LAMBDA=LBOUND                 !guess for first eigenvalue
      SPACNG=(UBOUND-LBOUND)/NDIM   !guess for eigenvalue spacing
C
      DO 200 L=1,NFIND              !loop to find NFIND eigenvalues
         DLAM=SPACNG                !initial guess for step
C
201      CONTINUE                   !coarse search to find upper
            LAMBDA=LAMBDA+DLAM      !bound for this eigenvalue
            CALL PLYNML(DIAG,LDIAG,LAMBDA,NDIM,COUNT)
         IF (COUNT .LT. L)  GOTO 201
C
         LAMBDA=LAMBDA-DLAM         !restart search
202      CONTINUE                   !search on a much finer scale
            CALL PLYNML(DIAG,LDIAG,LAMBDA,NDIM,COUNT)
            DLAM=DLAM/2             !next step is half as big
            IF (COUNT .LE. L-1) THEN
               LAMBDA=LAMBDA+DLAM   !not there yet, step forward
            ELSE
               LAMBDA=LAMBDA-DLAM   !went too far, step backward
            ENDIF
         IF (DLAM .GT. LTOLE)  GOTO 202
         EVAL(L)=LAMBDA
200   CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PLYNML(DIAG,LDIAG,LAMBDA,NDIM,COUNT)
C counts the eigenvalues less than LAMBDA
C by evaluating the terms in the sequence of polynomials used to eval
C the determinant and looking for changes in sign of those terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P5'
C Input Variables
      REAL  DIAG(MAXDIM),LDIAG(MAXDIM)   !diag,off-diag matrix elem
      INTEGER NDIM,NFIND        !matrix dim, number of values to find
      REAL LAMBDA               !guess for the eigenvalue
C Output Variables
      INTEGER COUNT             !number sign changes in the sequence
C Local Variables
      REAL TEMP,TEMP1,DET       !terms in the sequence
      INTEGER IDIM              !index                               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      TEMP1=1                   !first term in the sequence
      DET=DIAG(1)-LAMBDA        !second term
      IF (DET .LT. 0) THEN      !initialize COUNT
         COUNT=1
      ELSE 
         COUNT=0
      ENDIF
C
      DO 100 IDIM=2,NDIM         !recursion relation for determinant
         TEMP=(DIAG(IDIM)-LAMBDA)*DET-LDIAG(IDIM)**2*TEMP1
         TEMP1=DET               !roll values
         DET=TEMP
         IF (DET*TEMP1 .LT. 0.) COUNT=COUNT+1    !sign change?
         IF (TEMP1 .EQ. 0.)     COUNT=COUNT+1    !count zeros once
         IF (ABS(DET) .GE. 100000) THEN          !keep things small
            DET=DET/100000
            TEMP1=TEMP1/100000
         ENDIF
100   CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INVITR(DIAG,LDIAG,EVAL,NDIM,TMPVEC)
C finds the eigenvector TMPVEC corresponding to the eigenvalue EVAL
C of a symmetric, tri-diagonal NDIM x NDIM matrix with matrix elements
C DIAG and LDIAG 
C eigenvector is found by two inverse iterations of (E-H)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'PARAM.P5'
C Input Variables
      INTEGER NDIM                    !size of matrix
      REAL DIAG(MAXDIM),LDIAG(MAXDIM) !diagonal,off diagonal matrix elem
      REAL EVAL                       !eigenvalue 
C Output Variables
      REAL TMPVEC(MAXDIM)             !eigenvector
C Local Variables
      DOUBLE PRECISION A(MAXDIM,MAXDIM)!(E-H) and its inverse
      DOUBLE PRECISION NEWVEC(MAXDIM)  !temp storage for the evector
      INTEGER I,J                      !indices for matrix elements
      REAL NORM,SUM                    !variables for normalizing
      LOGICAL SINGLR                   !is the matrix singular?
      REAL TMPEPS                      !temp factor to avoid singularity
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      TMPEPS=EPS                       !initial value is input value
C     find the inverse of (E-H)
10    DO 100 I=1,NDIM
         DO 200 J=1,NDIM                  
200      A(I,J)=0.0                    !set A=(E-H)
         A(I,I)=(1.D+00+TMPEPS)*DBLE(EVAL)-DBLE(DIAG(I)) !diagonal elem
         IF (I .LT. NDIM) A(I,I+1)=-DBLE(LDIAG(I+1))  !off diagonal elem
         IF (I .GT. 1)    A(I,I-1)=-DBLE(LDIAG(I))
         TMPVEC(I)=1.                  !arbitrary initial vector
100   CONTINUE
      CALL MATINV(A,NDIM,SINGLR)       !find A inverse
      IF (SINGLR) THEN                 !if A is singular
          TMPEPS=TMPEPS*5.D+00         !try again with bigger EPS
          GOTO 10
      END IF                               
C
      DO 300 I=1,NDIM                  !first inverse iteration
         SUM=0.0
         DO 400 J=1,NDIM
400      SUM=SUM+A(I,J)*TMPVEC(J)
         NEWVEC(I)=SUM
300   CONTINUE
C
      NORM=0.0
      DO 700 I=1,NDIM                  !second inverse iteration
         SUM=0.0                       !and normalization
         DO 800 J=1,NDIM
800      SUM=SUM+A(I,J)*NEWVEC(J)
         TMPVEC(I)=SUM
         NORM=NORM+SUM*SUM
700   CONTINUE
      NORM=1/SQRT(NORM)
      DO 900 I=1,NDIM
        TMPVEC(I)=TMPVEC(I)*NORM
900   CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE MATINV(A,NDIM,SINGLR)
C inverts the matrix A of dimension NDIM using Gauss-Jordan elimination
C A is replaced by its inverse
C SINGLR is true if A is singular
C on exit, A is replaced by its inverse
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P5'
C Input/output variables:
      DOUBLE PRECISION A(MAXDIM,MAXDIM) !matrix for inverting
      INTEGER NDIM                 !dimension of matrix          
      LOGICAL SINGLR               !is the matrix singular?
C Local variables:                 
      DOUBLE PRECISION P(MAXDIM),Q(MAXDIM) !temporary storage
      LOGICAL ZEROED(MAXDIM)       !keeps track of who's zeroed
      INTEGER I,J,K                !indices
      INTEGER KBIG                 !index of row being zeroed
      DOUBLE PRECISION PIVOT       !largest diag element in A
      REAL BIG                     !largest diag elem of non zeroed rows
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SINGLR=.FALSE.
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

C        if all diagonals are zero, then the matrix is singular
         IF (KBIG .EQ. 0) THEN
            WRITE (OUNIT,*) ' Matrix is singular'
            SINGLR=.TRUE.
            RETURN
         END IF
C        matrix is ill conditioned if the size of elements varies greatly
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
      SUBROUTINE TXTOUT(JCHI,V,CEVEN,CODD,EVEVEN,EVODD,NLINES)
C calculates <Jz> and outputs results
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
      INCLUDE    'IO.ALL'      
      INCLUDE    'PARAM.P5'
C Input Variables
      REAL       CEVEN(MAXDIM,MAXDIM),CODD(MAXDIM,MAXDIM)  !eigenvectors
      REAL       EVEVEN(MXNCHI,MAXDIM),EVODD(MXNCHI,MAXDIM)!eigenvalues
      INTEGER    JCHI                                      !CHI index
      REAL       V                   !coupling strength
      INTEGER    NLINES              !number of lines printed out
C Local Variables
      INTEGER    I,K,IDIM       !indices
      INTEGER    M              !eigenvalue of Jz
      REAL       E,JZ           !eigenvalue and Jz
      INTEGER    DELINE         !number of lines this time round
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DELINE=3+JJ+1             !number of lines written per call
C     clear the screen at a convenient place
      IF ((NLINES+DELINE  .GT. TRMLIN-6) .AND. (TTERM)) THEN
          CALL PAUSE('to continue...',1)
          CALL CLEAR
          NLINES=0
      END IF
C
      IF (TTERM) THEN            !write out parameters
         WRITE (OUNIT,10) CHI(JCHI),V
         IF (VECTRS) THEN
             WRITE (OUNIT,20)
             WRITE (OUNIT,30)
         ELSE
             WRITE (OUNIT,50)
             WRITE (OUNIT,60)
         END IF
         NLINES=NLINES+3
      END IF
      IF (TFILE) THEN
         WRITE (TUNIT,10) CHI(JCHI),V
         IF (VECTRS) THEN
             WRITE (TUNIT,20)
             WRITE (TUNIT,30)
         ELSE
             WRITE (TUNIT,50)
             WRITE (TUNIT,60)
         END IF
      END IF
C
      DO 100 IDIM=1,JJ+1        !loop over all eigenvalues
         IF (MOD(IDIM,2) .EQ. 0) THEN     !putting them in order
            K=IDIM/2
            E=EVODD(JCHI,K)
            IF (VECTRS) THEN              !find the expectation value
               JZ=0                       ! of Jz if have eigenvectors
               DO 300 I=1,JJ
                  M=-JJ-1+2*I
                  JZ=JZ+M*CODD(I,K)**2    !contribution to <Jz>
300            CONTINUE
            ENDIF
         ELSE                             !even m states
            K=(IDIM+1)/2
            E=EVEVEN(JCHI,K)
            IF (VECTRS) THEN              !their contribution to <Jz>
               JZ=0
               DO 200 I=1,JJ+1
                  M=-JJ+2*(I-1)
                  JZ=JZ+M*CEVEN(I,K)**2
200            CONTINUE
            ENDIF
         ENDIF
C        write it all out
         IF (VECTRS .AND. TTERM) WRITE (OUNIT,40) E,JZ
         IF (VECTRS .AND. TFILE) WRITE (TUNIT,40) E,JZ
         IF (.NOT. VECTRS .AND. TTERM) WRITE (OUNIT,70) E
         IF (.NOT. VECTRS .AND. TFILE) WRITE (TUNIT,70) E
         IF (TTERM) NLINES=NLINES+1
         IF ((TTERM) .AND. (NLINES .EQ. TRMLIN-6)) THEN
            CALL PAUSE('to continue...',1)
            CALL CLEAR
            NLINES=0
         END IF
100   CONTINUE
C
10    FORMAT (23X,'Chi=',1PE11.4,5X,'V=',1PE11.4)
20    FORMAT (20X,'Eigenvalue',25X,'<Jz>')
30    FORMAT (20X,'----------',25X,'----')
40    FORMAT (17X,1PE15.8,17X,1PE15.8)                            
50    FORMAT (35X,'Eigenvalue')
60    FORMAT (35X,'----------')
70    FORMAT (32X,1PE15.8)
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
      INCLUDE 'PARAM.P5'
C Local parameters:
      CHARACTER*80 DESCRP          !program description
      DIMENSION DESCRP(20)
      INTEGER NHEAD,NTEXT,NGRAPH   !number of lines for each description
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL SETUP                !get environment parameters
C
C     display header screen     
      DESCRP(1)= 'PROJECT 5'
      DESCRP(2)= 'Solution of a Schematic Shell Model'
      NHEAD=2
C
C     text output description
      DESCRP(3)= 'Negative eigenvalues and <Jz>'
      NTEXT=1
C
C     graphics output description
      DESCRP(4)= 'Eigenvalues vs. coupling strength'
      NGRAPH=1
C 
      CALL HEADER(DESCRP,NHEAD,NTEXT,NGRAPH)
C 
      CALL MENU                      !setup constant part of menu
C                
      MTYPE(13)=NUM
      MPRMPT(13)='Enter the number of particles (must be even)'
      MTAG(13)='Number of particles'
      MLOLIM(13)=2
      MHILIM(13)=MAXN
      MINTS(13)=14
C
      MTYPE(14)=FLOAT
      MPRMPT(14)='Enter the intial value of coupling constant (Chi)'
      MTAG(14)= 'Initial coupling constant'
      MLOLIM(14)=0.
      MHILIM(14)=1000.
      MREALS(14)=1.0
C
      MTYPE(15)=FLOAT
      MPRMPT(15)='Enter increment in coupling constant (Chi)'
      MTAG(15)='Increment in coupling constant'
      MLOLIM(15)=0.
      MHILIM(15)=1000.
      MREALS(15)=.5
C
      MTYPE(16)=NUM
      MPRMPT(16)='Enter the number of values for coupling constant'
      MTAG(16)='Number of values of coupling constant'
      MLOLIM(16)=1
      MHILIM(16)=MXNCHI
      MINTS(16)=10
C
      MTYPE(17)=BOOLEN
      MPRMPT(17)='Do you want the eigenvectors <Jz> calculated?'
      MTAG(17)='Calculate eigenvectors'
      MINTS(17)=1
C
      MTYPE(18)=SKIP
      MREALS(18)=35
C
      MTYPE(38)=FLOAT
      MPRMPT(38)='Enter the tolerance for eigenvalue search'
      MTAG(38)='Tolerance for eigenvalue search'
      MLOLIM(38)=1.E-07
      MHILIM(38)=1.
      MREALS(38)=1.E-05
C
      MTYPE(39)=FLOAT
      MPRMPT(39)='Enter small factor to keep (E-H) nonsingular' 
      MTAG(39)='Factor to keep (E-H) nonsingular'
      MLOLIM(39)=1.E-07
      MHILIM(39)=1.
      MREALS(39)=1.E-05
C
      MTYPE(40)=SKIP
      MREALS(40)=60
C
      MSTRNG(MINTS(75))= 'proj5.txt'
C                      
      MTYPE(76)=SKIP
      MREALS(76)=80.
C 
      MSTRNG(MINTS(86))= 'proj5.grf'
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
      INCLUDE 'PARAM.P5'
C Local variables:
      INTEGER JCHI               !CHI index
C     map between menu indices and input parameters
      INTEGER INPART,ICHI,IDELCH,INCHI,IVEC,ILTOLE,IEPS,INGRF
      PARAMETER (INPART = 13)
      PARAMETER (ICHI   = 14)
      PARAMETER (IDELCH = 15)
      PARAMETER (INCHI  = 16)
      PARAMETER (IVEC   = 17)
      PARAMETER (ILTOLE = 38)
      PARAMETER (IEPS   = 39)
      PARAMETER (INGRF  = 87 )
C Functions:
      LOGICAL LOGCVT             !converts 1 and 0 to true and false 
      INTEGER GETINT             !get integer input from screen
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
      NPART=MINTS(INPART)
      CHI(1)=MREALS(ICHI)
      NCHI=MINTS(INCHI)
      DELCHI=MREALS(IDELCH)
      VECTRS=LOGCVT(MINTS(IVEC))
      LTOLE=MREALS(ILTOLE)
      EPS=MREALS(IEPS)
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
C     make sure that NPART is even (to avoid lots of bookkeeping)
10    IF (MOD(NPART,2) .NE. 0) THEN
         NPART=GETINT(NPART-1,2,MAXN,
     +            ' Number of particles must be even, try again')
         MINTS(INPART)=NPART
      GOTO 10
      END IF
      CALL CLEAR
C
C     calculated parameters
      JJ=NPART/2
      JJ1=JJ*(JJ+1)
      DO 100 JCHI=2,NCHI
         CHI(JCHI)=CHI(1)+(JCHI-1)*DELCHI
100   CONTINUE
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE PRMOUT(MUNIT,NLINES)
C outputs parameter summary to the specified unit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global variables:
       INCLUDE 'IO.ALL'
       INCLUDE 'PARAM.P5'
C Input variables:               
       INTEGER MUNIT            !unit number for output
       INTEGER NLINES           !number of lines printed out (I/O)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (MUNIT .EQ. OUNIT) CALL CLEAR
C 
       WRITE (MUNIT,2)
       WRITE (MUNIT,4)
       WRITE (MUNIT,6)NPART
       WRITE (MUNIT,10)JJ
C
       NLINES=4
C
2      FORMAT(1X,/)
4      FORMAT (' Output from project 5:',
     +     '  Solution of a schematic shell model')
6      FORMAT (1X,'Number of Particles = ',1I3)
10     FORMAT (1X,'Quasi Spin = ',1I3)
C
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE GRFOUT(DEVICE,EVEVEN,EVODD)
C graphs eigenvalues vs. CHI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global parameters:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P5'
      INCLUDE 'GRFDAT.ALL'
C Input variables
      INTEGER DEVICE                 !which device are we calling
      REAL EVEVEN(MXNCHI,MAXDIM),EVODD(MXNCHI,MAXDIM)!eigenvalues
C Local parameters:                         
      INTEGER JCHI,IDIM              !CHI and DIM indices
      CHARACTER*9 CNPART             !NPART as a character string
      INTEGER LNPART                 !length of that string
      REAL EVAL(MAXDIM)              !one eigenvalue at several CHI
      INTEGER SCREEN                 !send to terminal
      INTEGER PAPER                  !make a hardcopy
      INTEGER FILE                   !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     calculate parameters for graphing
      IF (DEVICE .NE. FILE) THEN
          NPLOT=1                       !how many plots
          IPLOT=1
C 
          YMAX=0.                       !limits
          YMIN=EVEVEN(NCHI,1)
          DO 10 JCHI=1,NCHI-1
             IF (EVEVEN(JCHI,1) .LT. YMIN) YMIN=EVEVEN(JCHI,1)
10        CONTINUE
          YMIN=YMIN+.01*(YMIN)
          XMIN=CHI(1)                   !leave a little extra room
          XMAX=CHI(NCHI)
          Y0VAL=XMIN
          X0VAL=YMIN
C 
          NPOINT=NCHI
C 
          ILINE=1                       !line and symbol styles
          ISYM=1
          IFREQ=1
          NXTICK=5
          NYTICK=5
C 
          CALL ICNVRT(NPART,CNPART,LNPART)     !titles
          TITLE='Schematic Shell Model with '
     +           //CNPART(1:LNPART)//' particles'
          LABEL(1)='CHI'
          LABEL(2)='eigenvalues'
C 
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
          CALL LNLNAX                          !draw axes
      END IF                                   
C 
C     output results
      DO 60 IDIM=1,JJ+1        !loop over all eigenvalues
        DO 50 JCHI=1,NCHI      !  and values of CHI
         IF (MOD(IDIM,2) .EQ. 0) EVAL(JCHI)=EVODD(JCHI,IDIM/2)
         IF (MOD(IDIM,2) .EQ. 1) EVAL(JCHI)=EVEVEN(JCHI,(IDIM+1)/2)
50      CONTINUE     
        IF (DEVICE .EQ. FILE)THEN
          WRITE (GUNIT,70)
          WRITE (GUNIT,80) (CHI(JCHI),EVAL(JCHI),JCHI=1,NCHI)
        ELSE
          CALL XYPLOT(CHI,EVAL)
        END IF
60    CONTINUE
C 
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)   !end graphics package
      IF (DEVICE .EQ. SCREEN) CALL TMODE         !switch to text mode
C 
70    FORMAT (5X,'    CHI       ',5X,'  Eigenvalue   ')
80    FORMAT (2(5X,1PE15.8))
100   FORMAT (/,' Patience, please; output going to a file.')
C 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE ONEOUT(DEVICE,EVEVEN,EVODD)
C graphs eigenvalues for one value of CHI 
C allows user to visually inspect eigenvalue spacings
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Global parameters:
      INCLUDE 'IO.ALL'
      INCLUDE 'PARAM.P5'
      INCLUDE 'GRFDAT.ALL'
C Input variables
      INTEGER DEVICE                 !which device are we calling
      REAL EVEVEN(MXNCHI,MAXDIM),EVODD(MXNCHI,MAXDIM)!eigenvalues
C Local parameters:                         
      INTEGER JCHI,IDIM              !CHI and DIM indices
      CHARACTER*9 CNPART             !NPART as a character string
      INTEGER LNPART                 !length of that string
      CHARACTER*9 CCHI               !CHI as a character string
      INTEGER LCHI                   !length of that string
      REAL EVAL(2),X(2)              !eigenvalue, dummy variable X
      INTEGER SCREEN                 !send to terminal
      INTEGER PAPER                  !make a hardcopy
      INTEGER FILE                   !send to a file
      DATA SCREEN,PAPER,FILE/1,2,3/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     messages for the impatient
      IF (DEVICE .NE. SCREEN) WRITE (OUNIT,100) 
C
C     calculate parameters for graphing
      IF (DEVICE .NE. FILE) THEN
          NPLOT=1                       !how many plots
          IPLOT=1
C 
          YMAX=0.                       !limits
          YMIN=(1.01)*EVEVEN(1,1)   
          XMIN=0.                       !x-axis is 'dummy' axis
          XMAX=10.
          Y0VAL=XMIN
          X0VAL=YMIN 
          X(1)=XMIN
          X(2)=1.
C 
          NPOINT=2
C 
          ILINE=1                       !line and symbol styles
          ISYM=1
          IFREQ=0
          NXTICK=5
          NYTICK=5
C 
          CALL ICNVRT(NPART,CNPART,LNPART)  !titles
          CALL CONVRT(CHI,CCHI,LCHI)
          TITLE='Schematic Shell Model with '
     +           //CNPART(1:LNPART)//' particles and CHI='//CCHI(1:LCHI)
          LABEL(1)=' '
          LABEL(2)='eigenvalues'
C 
          CALL GTDEV(DEVICE)                   !device nomination
          IF (DEVICE .EQ. SCREEN) CALL GMODE   !change to graphics mode
          CALL LNLNAX                          !draw axes
      END IF                                   
C 
C     output results
      DO 60 IDIM=1,JJ+1        !loop over all eigenvalues
        DO 50 JCHI=1,2
         IF (MOD(IDIM,2) .EQ. 0) EVAL(JCHI)=EVODD(1,IDIM/2)
         IF (MOD(IDIM,2) .EQ. 1) EVAL(JCHI)=EVEVEN(1,(IDIM+1)/2)
50      CONTINUE     
        IF (DEVICE .EQ. FILE) THEN
            WRITE (GUNIT,70)
            WRITE (GUNIT,80)  (X(JCHI),EVAL(JCHI),JCHI=1,2)
        ELSE
          CALL XYPLOT(X,EVAL)
        END IF
60    CONTINUE
C 
C     end graphing session
      IF (DEVICE .NE. FILE) CALL GPAGE(DEVICE)   !end graphic package
      IF (DEVICE .EQ. SCREEN) CALL TMODE         !switch to text mode
C 
70    FORMAT (5X,'    Dummy     ',5X,'  Eigenvalue   ')
80    FORMAT (2(5X,1PE15.8))
100   FORMAT (/,' Patience, please; output going to a file.')
C 
      RETURN
      END
                                          
