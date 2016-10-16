       DIMENSION PHI(0:200)

       EXACT(R)=1.-(R+2)*EXP(-R)/2
       SOURCE(R)=-R*EXP(-R)/2
       H=.1
       NSTEP=20./H
       CONST=H**2/12

       SM=0.
       SZ=SOURCE(H)
       PHI(0)=0
       PHI(1)=EXACT(H)
*       PHI(1)=.95*PHI(1)

       DO 10 IR=1,NSTEP-1
          R=(IR+1)*H
          SP=SOURCE(R)

          PHI(IR+1)=2*PHI(IR)-PHI(IR-1)+CONST*(SP+10.*SZ+SM)
          SM=SZ
          SZ=SP
10     CONTINUE

       SLOPE=(PHI(NSTEP)-PHI(NSTEP-10))/(10*H)
       DO 20 IR=1,NSTEP
          R=IR*H
          PHI(IR)=PHI(IR)-SLOPE*R
          DIFF=EXACT(R)-PHI(IR)
          PRINT *, R,EXACT(R),PHI(IR),DIFF
20     CONTINUE

       STOP
       END
