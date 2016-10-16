C subroutine to return 2 normally distributed numbers, GAUSS1 and GAUSS2
       SUBROUTINE DIST(GAUSS1,GAUSS2)
       INTEGER SEED
       DATA SEED/39249187/
       DATA PI/3.1415926/
  
       TWOU=-2.*LOG(1.-RAN(SEED))
       RADIUS=SQRT(TWOU)
       THETA=2*PI*RAN(SEED)
       GAUSS1=RADIUS*COS(THETA)
       GAUSS2=RADIUS*SIN(THETA)

       RETURN
       END
