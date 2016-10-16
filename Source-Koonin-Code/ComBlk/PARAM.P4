C param.p4
C
      INTEGER MAXN,MAXE              !limits on integration steps
      INTEGER MAXL                   !limits on partial waves
      INTEGER MAXANG                 !limits on number of angles
C
      PARAMETER (MAXN=1000)
      PARAMETER (MAXE=100)
      PARAMETER (MAXL=100)
      PARAMETER (MAXANG=100)
C
      INTEGER Z                      !charge of nucleus
      REAL E                         !energy of particle
      REAL VZERO                     !depth of potential
      INTEGER POT                    !which potential
      INTEGER LENZ,SQUARE,GAUSS      !types of potentials
      REAL V(MAXN+MAXE)              !potential function
      LOGICAL ATTRCT                 !is V attractive?
C     
      PARAMETER (LENZ=1)
      PARAMETER (SQUARE=2)
      PARAMETER (GAUSS=3)
C
      INTEGER NPTS,NXTRA             !number of integration points
      INTEGER LSTART,LSTOP           !range of partial waves
      INTEGER NANG                   !number of angles
C
      REAL E2                        !square of charge on electron
      REAL HBARM                     !Planck's constant**2/mass of elec
                                     !in eV * Angstroms**2
      REAL RMAX                      !maximum radius of potential
      REAL PI                        !constants
C
      REAL K,K2                      !wave number
      REAL Z6                        !sixth root of Z
      REAL DR                        !step size
      INTEGER LMAX                   !estimate of largest partial wave
      REAL RXTRA                     !value for r2-r1
      INTEGER LTABLE                 !how many Leg. pol. are in table
      REAL THETA(0:MAXANG),CTHETA(0:MAXANG) !theta and cos(theta)
      REAL DEGREE(0:MAXANG)          !theta in degrees
      REAL PL(0:MAXL,0:MAXANG)       !Legendre polynomials 
      REAL JL(2,0:MAXL),NL(2,0:MAXL) !spherical Bessel func at r1,r2
      REAL R(MAXN+MAXE)              !array of radii
C
      COMMON/PPARAM/Z,E,VZERO,POT
      COMMON/NPARAM/NPTS,NXTRA,LSTART,LSTOP,NANG
      COMMON/CONSTS/E2,HBARM,RMAX,PI                                        
      COMMON/PCALC/Z6,K,K2,DR,LMAX,RXTRA,V,LTABLE,PL,JL,NL,THETA,
     +        CTHETA,R,DEGREE,ATTRCT
