c     COMMON blocks for ZELIG version 2.3, modified for runon
c
c     Arrays are orderedso the innermost loop
c     index is the leftmost element of the array.
c
c     Parameters to define I/O logical units:
      integer RUNCON, INSITE, INSPP, PCHFILE, PRTFILE, TRAFILE,
     2  PROFILE, LOGFILE, INCFILE
      PARAMETER (RUNCON=2,INSITE=3,INSPP=4,PCHFILE=7,PRTFILE=8,
     2  TRAFILE=9,PROFILE=10,LOGFILE=11,INCFILE=12)
c     Parameters to dimension ZELIG version 2.3:
c     MR=max rows, MC=max columns in grid;
c     MT=max trees per plot;
c     MS=max number of species;
c     MH=max height of foliage profiles;
c     MST=max number of soil types on grid;
c     ML=max number of soil layers.
      PARAMETER (MR=10,MC=10,MT=200,MS=20,MH=85,MST=9,ML=10)
      character*4 msp
c     Variables in COMMON:
      COMMON/RUN/ mode, indata, nrows, ncols, nspp, nyrs, kyr, iprt,
     2     ilog, ipch, ilai, itrx, nmax, nsoils, mxs, idum(3)
c   msrt and ssmf added for water logging in this common block
      COMMON/SPP/ msp(MS), amax(MS), dmax(MS), hmax(MS), b2(MS), b3(MS),
     2     g(MS), lf(MS), ddmin(MS), ddmax(MS), ddf(MS),
     3     light(MS), mdrt(MS), msrt(MS), nutri(MS), seed(MS),
     4     nsprt(MS),sdmax(MS), smf(2,MS), ssmf(2,MS), nstmp(2,MS)
      COMMON/TREES/ isp(MT,MR,MC), dbh(MT,MR,MC), ibc(MT,MR,MC),
     2       nogro(MT,MR,MC)
      COMMON/PLOTS/ ind(MR,MC), msol(MR,MC), maxht(MR,MC), biom(MR,MC),
     2       sba(MS,MR,MC), sw(ML,MR,MC), 
     3       elai(12,MR,MC), ala(MH,MR,MC), al(0:MH,MR,MC)
      COMMON/SITE/ xtree, area, xk, theta, phib, phid, frunon,
     2       ct, tx, bgs, egs, tgs, xt(12), vt(12), xr(12), vr(12),
     3       b(12),ia(12),sun(12), t(12), r(12), runon(12), sff(3),
     4       nsl(MST), sf(MST), dl(ML,MST), fc(ML,MST), wp(ML,MST),
     5       wk(ML,MST), frft(ML,MST), frfs(ML,MST)

      
c     Definitions of variables in COMMON:
c
c     /CONTROL/ run configuration and control ...
c     MODE:  =0 if independent plots; =1 if interactive shading;
c     INDATA:  =0 to start from bare ground; 1, to read external initials;
c       =2 to write tree dump in final year (to be used as initials);
c     NROWS, NCOLS:  number of plots (grid rows, columns) to simulate;  
c     NSPP:  number of species in species driver file;
c     NYRS:  number of years to simulate (KYR, the counter);
c     IPRT:  ILOG, IPCH, ILAI, ITRX, interval for print, diagnostics log,
c       plot punch, leaf-area profile, and tracer file;
c     NMAX:  maximum number of trees per plot (set by AREA);
c     NSOILS:  number of soil types;
c     NSL:  number of soil layers per soil type (ML dimensions);
c     IDUM:  seed for random number generators; 
c     MXS:  height (m) of cross-section through leaf area profile (MODE=1)
c
c     /SPECIES/  dimensioned by parameter MS ...
c     MSP, AMAX, DMAX:  species mnemonic, max age (yr), and max dbh (cm);
c     HMAX, B2, B3:  max height (m) and coefficients of height allometry;
c     G, LF:  growth rate, lifeform;
c     DDMIN, DDMAX:  min, max degree-days;
c     LIGHT, MDRT, NUTRI:  shade-, drought-, nutrient response class;
c     SEED, NSPRT, SDMAX:  seedling establishment rate, sprouts, and
c       max sproutable diameter;
c     The following are derived and renewed for each plot each year:
c     DDF:  degree-day factor;
c     NSTMP:  number of new stumps from mortality (1=healthy, 2=stressed).
c
c     /TREES/  dimensioned by MT trees for plot (MR,MC) ...
c     ISP:  species label (A4);
c     DBH, IBC:  diameter (cm) and height (m) to base-of-crown;
c     DN:  demand for N for tree (sum of folaige, wood, and root N);
c     NOGRO:  number of no-growth years.
c
c     /PLOTS/  all dimensioned (at least) by plot (MR,MC) ... 
c     IND:  number of individual trees; 
c     MSOL:  mapped soil type (integer label);
c     MAXHT:  maximum canopy height;
c     BIOM:  above-ground woody biomass;
c     SBA:  basal area per species;
c     ALA:  actual leaf area per height per plot;
c     ELAI:  effective leaf area index per month, for interception;
c     AL:  available light at height m; 
c     SFF:  soil fertility factors, per response class;
c     SW:  soil water per plot and per soil layer;
c     SMF:  soil-moisture factor (1=topsoil, 2=average of all soil layers).
c   
c     /SITE/  site variables; some by soil layer, some by month ...
c     XTREE:  mean canopy area (sq.m) of typical full-sized tree;
c     AREA:  plot area (sq.m); SF:  soil fertility (Mg/ha/yr);
c     THETA, PHIB, PHID:  sun angle (rad), fraction beam and diffuse sun;
c     CT, TX:  coefficients for Priestley-Taylor PET estimate;
c     BGS, EGS:  day growing season begins and ends;
c     XT, VT; XR, VR:  mean temperature (C) and std. dev.,
c       mean total precip (cm) and std. dev.; 
c     T, R, RUNON:  annual temperature, rain, and runon (each year);
c     SUN:  solar radiation (W/sq.cm);
c     NSL:  number of soil layers per soil type;
c     SF:  soil fertility per soil type (Mg woody production /ha/year);
c     DL, FC, WP:  depth, field capacity, and wilting point (cm),
c       per soil layer (dimensioned by layer per soil type).
