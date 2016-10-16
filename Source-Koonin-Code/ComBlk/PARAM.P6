C param.p6                  
C
      INTEGER NXMAX,NYMAX             !maximum lattice size
      PARAMETER (NXMAX=100,NYMAX=100)         
      INTEGER MAXLEV                  !maximum number of contour levels
      PARAMETER (MAXLEV=20)
      INTEGER VRTCTY,STREAM           !display flags
C
      LOGICAL FIRST                   !first time through menu?
      INTEGER NX,NY                   !lattice size
      REAL VOMEGA                     !vorticity relaxation parameter
      REAL SOMEGA                     !stream relaxation parameter
      REAL MVOMEG                     !one minus vort omega
      REAL MSOMEG                     !one minus stream omega
      INTEGER HFWID                   !half width of plate
      INTEGER LENGTH                  !length of plate
      INTEGER FRONT                   !front edge of plate
      INTEGER BACK                    !back edge of plate
      REAL REYNLD                     !lattice Reynolds number
      REAL REYND4                     !Reynolds number / 4
C
      INTEGER PTYPE                   !choice for initial field
      CHARACTER*12 PFILE              !file to input init data 
      LOGICAL SUBLAT                  !is there a sublattice?
      INTEGER NXLL,NYLL,NXUR,NYUR     !sublattice parameters
C
      INTEGER NFREQ                   !frequency of display
      INTEGER NLEV                    !number of contour levels
C
      INTEGER BNDCND(NXMAX,NYMAX)     !flag to indicate bound cond
      !0=not a boundary; 1=boundary; 2=plate interior
C
      LOGICAL XSKIP,YSKIP             !skip spaces or lines in display
      INTEGER XCNTR,YCNTR             !how to center display
C
      COMMON / FLAG   / FIRST
      COMMON / PPARAM / NX,NY,REYNLD
      COMMON / NPARAM / VOMEGA,SOMEGA,SUBLAT,NXLL,NYLL,NXUR,NYUR,PTYPE
      COMMON / BCPRM / BNDCND,HFWID,LENGTH,FRONT
      COMMON / CPARAM / BACK,MVOMEG,MSOMEG,REYND4
      COMMON / GPARAM / NLEV,NFREQ,XSKIP,YSKIP,XCNTR,YCNTR              
      COMMON / ASCII / PFILE
