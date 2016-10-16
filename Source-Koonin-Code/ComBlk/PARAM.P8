C param.p8
C
      REAL S                 !inter-proton separation
      REAL S2                !S/2
      REAL BETA              !variational parameter
      INTEGER METHOD,CALC    !which method?  which quant to calculate?
      INTEGER NENSEM         !ensemble size
      REAL DT                !PIMC step size
      REAL DELTA             !step size in configuration space
      INTEGER NCORR          !max correlation length
      INTEGER NGROUP         !initial group size
      DOUBLE PRECISION DSEED !random number seed  
      INTEGER NTHERM         !number of thermalization sweeps
      INTEGER NFREQ          !freq of sweeps to avoid correlations
      INTEGER NSMPL          !size of groups
      LOGICAL TERSE          !terse output?
C
      REAL A,ALPHA           !constants in PHI
      REAL HBM               !hbar**2 divided by electron mass
      REAL E2                !electron charge squared 
      REAL ABOHR             !Bohr radius (Angstroms)
      REAL HBMDT             !hbar**2*dt/m
      REAL SQHBDT            !sqrt(hbar**2*dt/m)
C
      INTEGER MAXENS         !maximum ensemble size
      INTEGER MAXCRR         !max number of groups for correlation
      INTEGER NCOORD         !number of coordinates
      PARAMETER (MAXENS=100,MAXCRR=500,NCOORD=6)
C
      COMMON / PPARAM / S,BETA
      COMMON / NPARAM / DSEED,NTHERM,NFREQ,NSMPL,METHOD,CALC,NENSEM,
     +                  DT,DELTA,NCORR,NGROUP,TERSE
      COMMON / PCALC / A,S2,HBMDT,SQHBDT
      COMMON / CONST / E2,ABOHR,ALPHA,HBM
