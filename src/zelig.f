C     ZELIG version 2.3 by D.L. Urban, April 1993
C     Modified by M.F. Acevedo to include runon in order to account
C     for extra supply of water in bottomland conditions. Also modify
C     to have limited I/O and diganostics. No write statements in the
C     fortran. All writing to file is done by calling C functions. 
      
C     ZELIG is a gridded version of JABOWA/FORET.
C     In default mode, shading effects are implemented by integrating
C     available light over each tree's canopy (in 1-m increments).
C     In interactive-shading mode, the light regime is partitioned 
C     into diffuse and direct-beam radiation, shading is based
C     on sun angle and tree height, and plots shade each other.
C     Below-ground constraints include soil fertility (as maximum
C     annual woody production) and soil moisture (as drought-
C     days per growing season, for a multi-layer soil). 
C
C     HEIGHT equation LA function ALEAF are used in BOOKS and GROW
C     Optimal increment function DINCO is used in GROW; 
C     Function FERTF is sidestepped if a species is an N-fixer (NUTRI=0)
C     REGEN uses top-soil dry-days to compute SMF
C     Species are filtered through 3 years of weather before they
C     are planted as saplings.  Stump-sprouting is included by default.
C     WEATHR simulates weather
C     SOLWAT uses a Priestly-Taylor estimate of PET and an n-layer soil
C     (each with D, WP, FC); and uses canopy cover to partition PET into
c     transpiration and surface evaporation.
      
      subroutine ZELIG
      include 'zcommon.h'
      character*72 title
      character*20 driver1, driver2, driver3
c     I/O files:
      driver1 = 'z_control_inp.txt'
      driver2 = 'z_site_inp.txt'     
      driver3 = 'z_species_inp.txt'
c     Unit 2=RUNCON is run-control parameters:
      open(RUNCON,file=driver1,status='old')
c     Unit 3=INSITE is site parameters:
      open(INSITE,file=driver2,status='old')
c     Unit 4=INSPP is species driver-data:
      open(INSPP,file=driver3,status='old')
c     Read run title and control parameters
      read(RUNCON,2000) title, mode, nrows, ncols, nyrs, iprt, itrx
 2000 format(a72,10(/i5))

      call INITL
      do 1 kyr=1,nyrs
        call WEATHR
        do 2 kc=1,ncols
        do 2 kr=1,nrows
          call SOLWAT(kr,kc)
          call MORTAL(kr,kc)
          call GROW(kr,kc)
          call REGEN(kr,kc)
          call BOOKS(kr,kc)
    2   continue
        call GRID
    1   continue
        close(runcon)
        close(insite)
        close(inspp)
        end

c     INITL reads input data and initializes the model run.
      Subroutine INITL
      include 'zcommon.h'
      character*24 locale, species, cname
      character*20 soil(MST)
      data drzs /20.0/

c     Read site data: runon accounts for bottomland position
      read(INSITE,3001) locale, plat, plon, elev, theta, phib, phid,
     2  xk, xtree, xcht, nsoils, frunon
      nslm=0
      do 101 ksol=1,nsoils
        read(INSITE,3002) nsl(ksol), sf(ksol), soil(ksol)
        nslk=nsl(ksol)
        if (nslk.gt.nslm) nslm=nslk
        read(INSITE,3003) (dl(l,ksol), fc(l,ksol), wp(l,ksol),l=1,nslk)
  101   continue
 3001   format(a24/2f6.1,f7.1/3f6.2,/f6.3,/2f4.0,/i5,4f6.2)
 3002   format(i5,f6.2,2x,a20)
 3003   format(3f6.2)

c     Adjust for independent vs interactive-plot mode ...
c     Compute zone-of-influence from sun angle and max canopy height:
c       Shadow length for canopy tree:
        sl=xcht/TAN(theta)
c     No. of cells to aggregate:
        ts=SQRT(xtree)
        nac=INT(sl/ts+0.5)+1
      if (mode.eq.0) then
c     Independent plots:
        area=xtree*FLOAT(nac)
c       User override to hard-wire plot area (negate XTREE in driver):
          if (xtree.lt.0.0) area=ABS(xtree)
        xarea=FLOAT(nrows*ncols)*area/10000.0
        nmax=INT(area)
        else
c     Interactive-shading mode:
c     MXS= ht (m) of diagonal x.s. thru a plot at this sun angle:
        mxs=ts*TAN(theta)+0.5
        area=xtree
        nmax=INT(xtree)
        xarea=FLOAT(nrows*ncols)*area/10000.0
        aa=xtree*FLOAT(nac)
        endif

        do 102 ksol=1,nsoils
        nslk=nsl(ksol)
c     Correct SF units to plot size (Mg/ha->kg/plot)
        sf(ksol)=sf(ksol)*area/10.0
  102   continue
 8005   format(10x,i3,2x,f6.2,2x,2f6.2)

c     Read temperature and precip data (means and s.d.)
      read(INSITE,3004) (xt(mo),mo=1,12), (vt(mo),mo=1,12)
      read(INSITE,3004) (xr(mo),mo=1,12), (vr(mo),mo=1,12)
 3004 format(12f5.1)
      read(INSITE,3005) (sun(mo),mo=1,12), tl, tu
 3005 format(12f6.1,/2f6.2)

c     Parameters of gamma distribution for precipitation and initialize temp
        do 1111 mo=1,12
          b(mo)=xr(mo)/vr(mo)**2
          a=xr(mo)**2/vr(mo)**2
          ia(mo)=int(a+0.5)
          if (ia(mo).eq.0) ia(mo)=1
          t(mo)=xt(mo)
 1111     continue

c     Tally soil depths and Wk's:
        do 100 ksol=1,nsoils
          dtot=0.0
          do 110 l=1,nsl(ksol)
            dtot=dtot+dl(l,ksol)
            wk(l,ksol)=0.75*fc(l,ksol)
  110       continue
c     ... Set fine-root fraction per layer (over all layers):
          d1=0.0
          do 120 l=1,nsl(ksol)
            d2=d1+dl(l,ksol)
            frft(l,ksol)=(2.0*d2/dtot)*(1.0-d2/(2.0*dtot)) -
     &        (2.0*d1/dtot)*(1.0-d1/(2.0*dtot)) 
            if (frft(l,ksol).lt.0.0) frft(l,ksol)=0.0
            d1=d2
  120       continue
c     ... and for fine-root fractions for seedlings (to DRZS only):
          d1=0.0
          do 130 l=1,nsl(ksol)
            d2=d1+dl(l,ksol)
            if (d2.le.drzs) then
              frfs(l,ksol)=(2.0*d2/drzs)*(1.0-d2/(2.0*drzs)) -
     &          (2.0*d1/drzs)*(1.0-d1/(2.0*drzs)) 
              if (frfs(l,ksol).lt.0.0) frfs(l,ksol)=0.0
              else
              frfs(l,ksol)=0.0
              endif
            d1=d2
  130       continue
  100   continue

c     Read soils map for the grid:
      do 103 kr=1,nrows
        read(INSITE,3006) (msol(kr,kc),kc=1,ncols)
 3006   format(50i1)
  103   continue
c     Summarize grid size for run ...
      nplots=nrows*ncols

c     Compute constants for PET estimates:
      e1=33.8639*((0.00738*tl+0.8072)**8.0 -
     &  0.000019*ABS(1.8*tl+48.0)+0.001316)
      e2=33.8639*((0.00738*tu+0.8072)**8.0 -
     &  0.000019*ABS(1.8*tu+48.0)+0.001316)
      ct=1.0/(38.0-2.0*elev/305.0+380.0/(e2-e1))
      tx=-2.5-0.14*(e2-e1)-elev/550.0
c     IDUM is the seed for the random number generator; seed (1) goes to
c     WEATHR, (2) to REGEN, and (3) to MORTAL:
c     WEATHER uses RAN1, GAUSS1, and GAMMA; REGEN uses RAN2 and GAUSS2;
c     and MORTAL uses RAN3.
      idum(1)=-1
      idum(2)=-1
      idum(3)=-1
c     Start simulation with a spin of the random-number generator:
c     Enter seed for random weather start-up [0 to skip]:  '
      istart = 1
      if (istart.gt.0) then
        do 105 n=1,istart
          y=RAN1(idum(1))
  105     continue
        endif

c     Initialize plot arrays ...
      do 10 kc=1,ncols
      do 11 kr=1,nrows
        ind(kr,kc)=0
        biom(kr,kc)=0.0
        maxht(kr,kc)=2
c     ... Tree arrays:
        do 111 ki=1,MT
          isp(ki,kr,kc)=0
          dbh(ki,kr,kc)=0.0
          ibc(ki,kr,kc)=0
          nogro(ki,kr,kc)=0
  111     continue
c     ... Light profile:
        al(0,kr,kc)=1.0
        do 112 m=1,MH
          al(m,kr,kc)=1.0
          ala(m,kr,kc)=0.0
  112     continue
c     ... Effective leaf area for interception:
        do 113 mo=1,12
          elai(mo,kr,kc)=0.0
  113     continue
c     ... Soil water per layer:
        ksol=msol(kr,kc)
        do 114 l=1,nsl(ksol)
          sw(l,kr,kc)=fc(l,ksol)
  114     continue
c     ... Species basal area:
        do 115 ks=1,MS
          sba(ks,kr,kc)=0.0
  115     continue
c     ... and soil fertility factors:
   11   continue
   10   continue
        do 116 l=1,3
          sff(l)=1.0
  116     continue

c     Read and compute species parameters ...
c     Read number of species ...
      read(INSPP,4001) nspp
 4001 format(i4)
      kseed=0
c     LSpp is local species included (non-0 seed)
      lspp=0
c     Seed(k) is set from assumption of 1 stem/sq.m, and a 10-yr
c       stocking time; rank seeding rates are then adjusted to this.
      smax=area/10.0
      stot=0.0
      do 12 k=1,nspp
        read(INSPP,4002) msp(k), species, cname, amax(k), dmax(k), 
     2    hmax(k), b2(k), b3(k), g(k), lf(k), ddmin(k), ddmax(k),  
     3    light(k), mdrt(k), nutri(k), seed(k), nsprt(k), sdmax(k)
 4002   format(a4,2x,2a24,/f5.0,f4.0,f7.3,f8.4,f7.4,f5.0,i3,1x,2f5.0,
     2    1x,3i2,f3.0,i2,f4.0)
        if (k.eq.1.and.seed(k).lt.0.0) kseed=1
        if (kseed.eq.1.and.seed(k).ne.0.0) seed(k)=1.0
        if (seed(k).gt.0.0) then
          lspp=lspp+1
          stot=stot+seed(k)
          endif
c     Initialize stumps, soil moisture, and degree-day arrays:
        do 121 l=1,2
          nstmp(l,k)=0
          smf(l,k)=1.0
  121     continue
        ddf(k)=1.0
   12   continue
c     Relativize seed(k) and adjust to plot size:
      do 13 k=1,nspp
        if (seed(k).eq.0.0) go to 13
        seed(k)=smax*seed(k)/stot
   13   continue
      kyr=0
      call GRID
      return
      end

c     BOOKS keeps track of biomass, leaf area, and available light
      Subroutine BOOKS(kr,kc)
      include 'zcommon.h'
c     CP, compensation point for 5 tolerance classes:
      real cp(5)
c     Month that growing season begins, ends (for ELAI, used in SOLWAT);
c       MCH, max canopy height (used if MODE=0 to compute light profile).
      integer mgsb, mgse, mch
      data cp/0.15,0.12,0.09,0.06,0.03/

c     Initialize, and set months of growing season:
        mgsb=INT(bgs/30.0)+1
        mgse=INT(egs/30.0)+1
          if (mgse.gt.12) mgse=12  
        mch=2
        bmi=0.0
        biom(kr,kc)=0.0
        ksol=msol(kr,kc)
c     ... and initialize actual LAI:
      do 301 m=1,mh
        ala(m,kr,kc)=0.0
  301   continue
c     ... and effective LAI per month, for interception:
      do 302 mo=1,12
        elai(mo,kr,kc)=0.0
  302   continue
c     ... and basal area per species:
      do 303 ks=1,nspp
        sba(ks,kr,kc)=0.0
  303   continue
      if (ind(kr,kc).eq.0) go to 33

c     Individual tree loop ...
      do 32 ki=1,ind(kr,kc)
        ks=isp(ki,kr,kc)
        kf=lf(ks)
        d=dbh(ki,kr,kc)
c     Compute tree height and hang leaves:
        ht=HEIGHT(d,hmax(ks),b2(ks),b3(ks)) 
        iht=INT(ht+0.5)
        if (iht.gt.mh) iht=mh
        if (iht.gt.mch) mch=iht
        mbc=ibc(ki,kr,kc)
        hc=FLOAT(mbc)
c     If a tree's lost all its leaves to shading, skip it:
        if (hc.ge.ht) go to 32
c     Hang leaves in canopy:
c       fd=unit foliage density; f, as lai: 
        cl=FLOAT(iht-mbc+1)
        fd=ALEAF(d,ht,hc,kf)/cl
        f=fd/area
        cpl=cp(light(ks))
        do 321 m=iht,mbc,-1
c     Stop when out of light or below old canopy:
          if (al(m,kr,kc).lt.cpl) go to 322
          ala(m,kr,kc)=ala(m,kr,kc)+f
  321     continue
  322   nbc=m+1
        tla=fd*FLOAT(iht-m)
c     Reassign ibc:
        ibc(ki,kr,kc)=nbc
c     Effective LAI for interception is 10% of TLA for deciduous forms
c       in the offseason (to account for some branch/stem surface):
        do 323 mo=1,12
          sai=tla/area
          if ((mo.lt.mgsb.or.mo.gt.mgse).and.kf.gt.5) sai=sai*0.10
          elai(mo,kr,kc)=elai(mo,kr,kc)+sai
  323     continue
c     Increment plot woody biomass (kg):
        wb=WOOD(d,kf)
        biom(kr,kc)=biom(kr,kc)+wb
c     Potential biomass increment from current and potential dbh: 
        di=DINCO(g(ks),tla,d,ht,dmax(ks),hmax(ks),b2(ks),b3(ks))
        dx=d+di
        wbi=WOOD(dx,kf)-wb
        bmi=bmi+wbi
c     Increment basal area:
        sba(ks,kr,kc)=sba(ks,kr,kc)+3.14159*(d/200.0)**2
   32   continue

c     Compute available light for independent-plot mode ...
   33   al(mch,kr,kc)=1.0
        clai=ala(mch,kr,kc)
        maxht(kr,kc)=mch
        do 331 m=(mch-1),1,-1
          al(m,kr,kc)=EXP(-xk*clai)
          clai=clai+ala(m,kr,kc)
  331     continue
        al(0,kr,kc)=EXP(-xk*clai)

c     Compute soil fertility multipliers:
        if (bmi.gt.sf(ksol)) then
          rf=sf(ksol)/bmi
          else
          rf=1.0
          endif
        do 332 l=1,3
          sff(l)=FERTF(l,rf)
  332     continue
      return
      end

c     WEATHR generates temperature, DegD, and precipitation
c     This is weather version 2.2
      Subroutine WEATHR
      include 'zcommon.h'
      real ddbase, rtot, fs, stot
      real days(12)
      external GAMMA
      data days/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./
c     DDbase, base (C) for degree-day heat sum;
      data ddbase /5.56/
c     Ptot, total annual precipitation; Rtot, total annual rain
c     Degd, degree-days;  bgs (egs), beginning (end) of growing season;
c       (igs, a flag to catch the season's first warming);
      degd=0.0
      ptot=0.0
      rtot=0.0
      stot=0.0
      bgs=1.0
      egs=365.0
      igs=0

c     Generate weather:
      do 12 mo=1,12
c     Temperatures are normally distributed:
        t(mo)=xt(mo)+GAUSS1(idum(1))*vt(mo)
c     Ppt is from a gamma distribution:
        z=GAMMA(ia(mo),idum(1))
        r(mo)=z/b(mo)
        if (r(mo).lt.0.0) r(mo)=0.0
        ptot=ptot+r(mo)
c     fs artifically constant to emulate flooding by runon
        fs=frunon
        runon(mo)=r(mo)*fs
        rtot=rtot+r(mo)
        stot=stot+runon(mo)

c     Find bgs, egs, and total growing season length tgs
        lm=mo-1
        if (lm.eq.0) lm=12
        if (t(lm).lt.ddbase.and.t(mo).ge.ddbase.and.igs.eq.0) then
          if (mo.eq.1) then
            bgs=1.0
          else
            bgs=TLINE(lm,t(lm),mo,t(mo),ddbase)
            endif
          igs=1
          endif
        if (mo.gt.6.and.t(lm).ge.ddbase.and.t(mo).lt.ddbase)
     2    egs=TLINE(lm,t(lm),mo,t(mo),ddbase)
c     Compute growing degree-days:
        if (t(mo).gt.ddbase) degd=degd+(t(mo)-ddbase)*days(mo)
   12   continue
      tgs=egs-bgs+1
        if (tgs.le.0.0) tgs=tgs+365.0

c     Compute species degree-day factors:
      do 13 ks=1,nspp
        ddf(ks)=DEGDF(ddmin(ks),ddmax(ks),degd)
   13   continue

      return
      end

c     Leaf area (sq.m), on dbh (cm), keyed by tree form, and
c     adjusted for height-to-crown, via sapwood area.
      Real Function ALEAF(d,ht,hc,lf)
      real d, rc, ht, hc, pi
      integer lf
      real slr(9), t1(9), t2(9), t3(9), s1(9), s2(9), b(9)
c     Leaf form labels: 1=ABsp=THpl, 2=PIpo, 3=PSme, 4=TShe, 5=PIen, 
c       6=ACma, 7=ALru, 8=QUal (=ACma), 9=average (=ACma):
c     SLR: sapwood:leaf area ratio (Waring and Schlesinger 1985)
        data slr /.56, .18, .54, .46, .45, .21, .10, .40, .33/
c     T's are taper coefficients (Kozak et al. 1969 and S. Garman)
      data t1 /0.96904,0.90895,0.87150,1.13877,1.14542,
     &  0.95997,0.97576,0.95997,0.95997/
      data t2 /-1.64733,-1.46990,-1.50058,-1.66358,-2.10213,
     &  -1.46336,-1.22922,-1.46336,-1.46336/
      data t3 /0.67829,0.56095,0.62908,0.52481,0.95671, 
     &  0.50339,0.25347,0.50339,0.50339/
c     S's are sapwood coefficents (regressed from various data)
c       (all hardwoods use pipo values) 
      data s1 /4.67,18.19,5.43,16.3,5.72,15.1,15.1,15.1,15.1/
      data s2 /.0341,.0108,.0460,.0178,.0393,.0314,.0314,.0314,.0314/
c     B's are bark coefficients:  b inside bark=b*dbh outside bark ...
c       (all hardwoods are faked at .95)
      data b /.95, .95, .89, .96, .96, .95, .95, .95, .95/
      data pi /3.14159/
c     LF=0:  use original FORET allometry (or internal default) ...
c       ALEAF is total leaf area multiplied by crown ratio:
      if (lf.eq.0) then
        ALEAF=0.160694*d**2.129
        ALEAF=ALEAF*(ht-hc+1.0)/ht
        return
        endif
c     All explicit tree forms ...
c     Compute D at base of crown, using generic taper equations;  
c       then use slr ratios to convert this area (sq cm)
c       to leaf area (sq m) ...
c     Estimate d at crown base=dc, radius=dc/2
c       (if no pruning, =dbh-bark thickness)
      if (hc.gt.1.0) then
        dc=d*SQRT(t1(lf)+t2(lf)*(hc/ht)+t3(lf)*(hc**2/ht**2))
        else
        dc=b(lf)*d
        endif
        rc=dc/2.0
c     Compute sapwood width sw from d:
        sw=s1(lf)*(1.0-EXP(-s2(lf)*d))
c     Check to make sure sw isn't > rc:
        if (sw.gt.rc) sw=rc
c     Compute sapwood x.s.area by differencing:
        sa=pi*rc**2-pi*(rc-sw)**2
        ALEAF=slr(lf)*sa
      return
      end

c     Available light factor, by tolerance class and light:
      Real Function ALF(kt,al)
      real c1(5), c2(5), c3(5)
c     Parameters revised January 1992;
c     Tolerance classes reversed March 1993 
c     (1=very intolerant; 5=very tolerant).
      data c1 /1.57851,1.25977,1.12598,1.04689,1.02046 /
      data c2 /1.18855,1.78588,2.43920,3.29031,4.16533 /
      data c3 / 0.15,0.12,0.09,0.06,0.03/
      ALF=c1(kt)*(1.0-EXP(-c2(kt)*(al-c3(kt))))
        if (ALF.lt.0.0) ALF=0.0
      return
      end
c     AMORT is for ambient (natural) mortality
c
      Logical Function AMORT(agemax,idum)
      integer idum
      real agemax
      if (RAN3(idum).lt.(4.605/agemax)) then
        AMORT=.true.
        else
        AMORT=.false.
        endif
      return
      end

c     Degree-day Factor, on min, max, and actual degree-days
      Real Function DEGDF(ddmin,ddmax,degd)
      real ddmin, ddmax, degd
      DEGDF=4.0*(degd-ddmin)*(ddmax-degd)/(ddmax-ddmin)**2
      if (DEGDF.lt.0.0) DEGDF=0.0
      return
      end
c     Diameter growth as a function of leaf area, etc.
      Real Function DINCO(g,tl,d,h,dmax,hmax,h2,h3)
      real g, tl, d, h, dmax, hmax, h2, h3
      real hc, h1, c
c     Convert heights from m to cm; c is a shorthand constant.
        hc=100.0*h
        h1=100.0*hmax
      c=EXP(h2*d)
      DINCO=g*tl*(1.0-d*hc/(dmax*h1)) /
     2  (d*h1*(-1.0*d*h3*h2*c*(1.0-c)**(h3-1.0)+2.0*(1.0-c)**h3))
      return
      end

c     DLA returns a diagonal leaf-area profile through ALA,
c       with vertical step-size lxs (m),
c       in look-direction ldir (1-4=nesw; 0=vertical)
c       (called only if model in is interactive-shading mode)
c
      Real Function DLA(ir,ic,iht,lxs,ldir)
      include 'zcommon.h'
      real xla
      integer ir, ic, iht, jr, jc, lxs, jxs, ldir
        xla=0.0
        jr=ir
        jc=ic
        jxs=lxs
c     Vertical profile:  set top-height to top-of-canopy ...
      if (ldir.eq.0) jxs=MH+1
        m1=iht
        m2=m1+jxs
        do 101 mm=m1,m2-1
          if (mm.gt.MH) go to 105
          xla=xla+ala(mm,jr,jc)
  101     continue
c     Diagonal profiles:  rove until XLA pops out of top-of-canopy,
c       wrapping grid at edges ...
  102 if (ldir.eq.1) then
        jr=jr-1
        if (jr.lt.1) jr=jr+nrows
        endif
      if (ldir.eq.2) then
        jc=jc+1
        if (jc.gt.ncols) jc=jc-ncols
        endif
      if (ldir.eq.3) then
        jr=jr+1
        if (jr.gt.nrows) jr=jr-nrows
        endif
      if (ldir.eq.4) then
        jc=jc-1
        if (jc.lt.1) jc=jc+ncols
        endif
      m1=m2
      m2=m1+jxs
      do 103 mm=m1,m2-1
        if (mm.gt.mh) go to 105
        xla=xla+ala(mm,jr,jc)
  103   continue
      go to 102
  105 DLA=xla
      return
      end
c     Drought Factor, given max drought tolerance and dry-days:
c
      Real Function DRTF(mdt,ddays)
      real dt, ddays, drt
      integer mdt
      dt=FLOAT(mdt)/10.0
      drt=MIN(dt,ddays)
      DRTF=SQRT((dt-drt)/dt)
      return
      end
c     Fertility factor, by nutrient response class (1=intol, 3=tol)
c     and soil fertility (relative, on [0,1]).
c
      Real Function FERTF(nrc,sf)
      integer nrc
      real b1(3), b2(3), b3(3)
      real sf
      data b1 /1.03748,1.00892,1.01712/
      data b2 /-4.02952,-5.38804,-4.12162/
      data b3 /0.17588,0.12242,0.00898/
      if (nrc.eq.0) then
        FERTF=1.0
        else
        FERTF=b1(nrc)*(1.0-exp(b2(nrc)*(sf-b3(nrc))))
        if (FERTF.lt.0.0) FERTF=0.0
        endif
      return
      end
c     GAMMA(ia,idum) returns a random deviate, z, from a 1-parameter,
c       standardized gamma distribution with integer shape parameter ia. 
c     This routine is from Press et al. (1987), 
c       as implemented by Gordy Bonan in his boreal forest model.
c
      Real Function GAMMA(ia,idum) 
      if (ia.lt.6) then
    5    x=1. 
         do 10 j=1,ia 
         x=x*RAN1(idum)
   10    continue 
         if (x.lt.(1.e-6)) go to 5
         x=-ALOG(x)
      else
   20    v1=2.*RAN1(idum)-1. 
         v2=2.*RAN1(idum)-1.
         if ((v1**2.+v2**2.).gt.(1.)) go to 20
         y=v2/v1
         am=ia-1
         s=SQRT(2.*am+1.)
         x=s*y+am
         if (x.le.(0.)) go to 20
c        domain of exp function is (-675.81, 741.66)
         coeff=am*ALOG(x/am)-s*y
         coeff=AMAX1(coeff,-675.81)
         coeff=AMIN1(coeff,741.66)
         e=(1.+y**2.)*EXP(coeff)
         if (RAN1(idum).gt.e) go to 20 
      endif
      GAMMA=x
      return
      end

c     GAUSS1(idum) returns a normal random number, x,s=(0,1).
c     This function is dedicated to WEATHR (with RAN1).
      Real Function GAUSS1(idum)
      data iset/0/
      if (iset.eq.0) then
   21   v1=2.0*RAN1(idum)-1
        v2=2.0*RAN1(idum)-1
        r=v1**2+v2**2
        if (r.ge.1.) go to 21
        fac=SQRT(-2.*ALOG(r)/r)
        gset=v1*fac
        GAUSS1=v2*fac
        iset=1
      else
        GAUSS1=gset
        iset=0
        endif
      return
      end

c     GAUSS2(idum) returns a normal random number, x,s=(0,1).
c     GAUSS2 is dedicated to REGEN (with RAN2).
      Real Function GAUSS2(idum)
      data iset/0/
      if (iset.eq.0) then
   21   v1=2.0*RAN2(idum)-1
        v2=2.0*RAN2(idum)-1
        r=v1**2+v2**2
        if (r.ge.1.) go to 21
        fac=SQRT(-2.*ALOG(r)/r)
        gset=v1*fac
        GAUSS2=v2*fac
        iset=1
      else
        GAUSS2=gset
        iset=0
        endif
      return
      end

c     GRID does all the between-plot interactions and aggregation
      Subroutine GRID
      include 'zcommon.h'
c     XBA, total basal area per species over all plots
      real xba(MS)

c     Initialize totals for tracer ...
c     B is biomass; D, density; BA, basal area, L, LAI; H, canopy height:
        sumb=0.0
        ssqb=0.0
        sumd=0.0
        suml=0.0
        sumh=0.0
        ss=FLOAT(nrows*ncols)
        acf=10000.0/area
      do 301 ks=1,nspp
        xba(ks)=0.0
  301   continue
c     Set constants of interactive shading:
        pd=phid/5.0
        km=10
      do 30 kc=1,ncols
      do 31 kr=1,nrows

c     Do diagonal leaf area profiles if in interactive mode ...
c       DLA returns a diagonal LAI profile in steps thru array ALA;
c       directions 1-4 are N,E,S, and W ...
c     Direct-beam insolation accrues in steps of MXS southward;
c       diffuse light is tallied across a sky arc in 9 samples ...
c       ... as 1 DLA in each direction, with height 10 m ...
c       ... + the vertical profile (direction 0):
        if (mode.eq.1) then
          do 311 m=0,(maxht(kr,kc)-1)
            xla=DLA(kr,kc,m+1,mxs,3)
            alb=phib*EXP(-xk*xla)
            ald=0.0
            do 312 ldir=0,4
              xla=DLA(kr,kc,m+1,km,ldir)
              ald=ald+pd*EXP(-xk*xla)
  312         continue
c     Total light is direct-beam + diffuse:
            al(m,kr,kc)=alb+ald
  311       continue
          endif

c     Tallies for tracer ...
c     ... Density:
        sumd=sumd+FLOAT(ind(kr,kc))
c     ... Biomass (Mg):
        bm=biom(kr,kc)/1000.0
        sumb=sumb+bm
        ssqb=ssqb+bm**2
c     ... LAI:
        suml=suml+DLA(kr,kc,1,MH,0)
c     ... Canopy height:
        sumh=sumh+FLOAT(maxht(kr,kc))
c     ... Basal area (sq.m) per species:
        do 313 ks=1,nspp
          xba(ks)=xba(ks)+sba(ks,kr,kc)
  313     continue
   31 continue
   30 continue

c     Print aggregate stand summary:
        if (MOD(kyr,iprt).eq.0) call PRINTING

c     Tracers, averaged and converted to per-ha:
      if (MOD(kyr,itrx).eq.0) then
        xb=(sumb/ss)*acf
        sdb=SQRT((ssqb-sumb**2/ss)/(ss-1.0))*acf
        xd=(sumd/ss)*acf
        xl=suml/ss
        xh=sumh/ss
        xbatot=0.0
        do 302 ks=1,nspp
          xba(ks)=(xba(ks)/ss)*acf
          xbatot=xbatot+xba(ks)
  302     continue
c      call C function to perform IO         
        call ztracer(nspp, kyr, xd, xb, sdb, xbatot, xl, xh, xba) 
        endif

      return
      end

c     GROW increments tree diameters
      Subroutine GROW(kr,kc)
      include 'zcommon.h'

c     Individual tree loop:
      do 62 ki=1,ind(kr,kc)
        ks=isp(ki,kr,kc)
        d=dbh(ki,kr,kc)
        mbc=ibc(ki,kr,kc)
c     Compute tree height:
        ht=HEIGHT(d,hmax(ks),b2(ks),b3(ks)) 
        iht=INT(ht+0.5)
        if (iht.gt.mh) iht=mh
c     Available light factor:
        algf=0.0
        if (FLOAT(mbc).gt.ht) then
          tla=0.0
          sfgf=0.0
          smgf=0.0
          gf=0.0
          dinc=0.0
          nogro(ki,kr,kc)=nogro(ki,kr,kc)+1
          go to 62
          endif
  620   do 621 m=mbc,iht
          algf=algf+ALF(light(ks),al(m,kr,kc))
  621     continue
        cl=FLOAT(iht-mbc+1)
        algf=algf/cl
c     Tree's leaf area:
        tla=ALEAF(d,ht,FLOAT(mbc),lf(ks))
c     Assign SFF (N-fixers are immune to soil fertility):
        if (nutri(ks).eq.0) then
          sfgf=1.0
          else
          sfgf=sff(nutri(ks))
          endif
c     Assign soil moisture growth factor:
        smgf=smf(2,ks)
c     Optimal growth increment as f(d,g, ..., and leaf area):
        dm=DINCO(g(ks),tla,d,ht,dmax(ks),hmax(ks),b2(ks),b3(ks))
c     Growth factor is light*min(water,nutrients)*temperature:
        gf=algf*min(smgf,sfgf)*ddf(ks)
c     Realized increment:
        dinc=dm*gf
        dbh(ki,kr,kc)=dbh(ki,kr,kc)+dinc
c     NoGro, if growth<10% of opt, or dinc<.1 mm
        if (gf.lt.0.10.or.dinc.lt.0.01) then
          nogro(ki,kr,kc)=nogro(ki,kr,kc)+1
          else
          nogro(ki,kr,kc)=0
          endif
   62   continue
      return
      end

c     Tree height (m), allometric on dbh (cm), per species
      Real Function HEIGHT(d,h1,h2,h3)
      real d, h1, h2, h3
      HEIGHT=h1*(1.0-EXP(h2*d))**h3
        if (HEIGHT.lt.1.37) HEIGHT=1.37
      return
      end

c     MORTAL kills trees.
      Subroutine MORTAL(kr,kc)
      include 'zcommon.h'
      integer live(15), ndead(2,15)
      logical AMORT, SMORT
c     IDead and ndead tally dead trees for diagnostics
      idead=0
      do 411 l=1,15
        live(l)=0
        ndead(1,l)=0
        ndead(2,l)=0
  411   continue
c     Initialize stump arrays:
      do 412 ks=1,nspp
        nstmp(1,ks)=0
        nstmp(2,ks)=0
  412   continue
c     Individual tree loop:
        if (ind(kr,kc).eq.0) return
      do 42 ki=1,ind(kr,kc)
        ks=isp(ki,kr,kc)
        d=dbh(ki,kr,kc)
c     Tally stems in dbh classes (for diagnostics)
        l=INT(d/10.0)+1
          if (l.gt.15) l=15
        live(l)=live(l)+1
c     Mortality can come from stress, or ambient ...
c     Unhealthy tree dies:
        if (SMORT(nogro(ki,kr,kc),idum(3))) then
c     Label dbh and tally 1 dead from stress
          dbh(ki,kr,kc)=-2.0
          ndead(2,l)=ndead(2,l)+1
          idead=idead+1
c     Check if sproutable:
          if (d.ge.5.0.and.d.le.sdmax(ks)) then
            nstmp(2,ks)=nstmp(2,ks)+1
            endif
          go to 42
          endif
c     Healthy tree dies:
        if (AMORT(amax(ks),idum(3))) then
c     Label dbh and tally 1 dead of natural causes
          dbh(ki,kr,kc)=-1.0
          ndead(1,l)=ndead(1,l)+1
          idead=idead+1
c     Check if sproutable:
          if (d.ge.5.0.and.d.le.sdmax(ks)) then
            nstmp(1,ks)=nstmp(1,ks)+1
            endif
          endif
   42   continue

c     Reshuffle tree arrays if there was any mortality:
        if (idead.eq.0) return
      in=0
      do 43 io=1,ind(kr,kc)
        if (dbh(io,kr,kc).eq.0.0) go to 44
        if (dbh(io,kr,kc).lt.0.0) go to 43
        in=in+1
        dbh(in,kr,kc)=dbh(io,kr,kc)
        nogro(in,kr,kc)=nogro(io,kr,kc)
        isp(in,kr,kc)=isp(io,kr,kc)
        ibc(in,kr,kc)=ibc(io,kr,kc)
   43   continue
   44 ind(kr,kc)=in
c     Void remainder of tree arrays
      inxs=in+1
      do 45 ix=inxs,nmax
        dbh(ix,kr,kc)=0.0
        nogro(ix,kr,kc)=0
        isp(ix,kr,kc)=0
        ibc(ix,kr,kc)=0
   45   continue
      return
      end

c     PRINTING writes aggregate stand statistics for all plots
      Subroutine PRINTING
      include 'zcommon.h'
      real ba(MS), sd(MS), freq(MS), dc(10), sdc(10,MS), tsdc(10)
      real rsd(MS), rba(MS), riv(MS)
      integer inc(MS,MR,MC)
      character*5 msp1(MS)

      plots=FLOAT(nrows*ncols)
      xarea=plots*area/10000.0
c     Variable definitions ...
c       ITD is (integer) total density; ITB, density of trees>10 cm;
c       ISD, SBA, RIV are species density, basal area, and importance;
c       IDC is stem tally by 10-cm dbh class (ISDC, by spp).
c       T* is total (density, BA, biomass, ...), which are converted to
c       mean values X*.
c     Initialize to zero:
        td=0.0
        tb=0.0
        tba=0.0
        tbiom=0.0
        tla=0.0
        tht=0.0
        tdbh=0.0
        ssqd=0.0
        sdd=0.0
        xdbh=0.0
        do 701 ks=1,nspp
          ba(ks)=0.0
          sd(ks)=0.0
          freq(ks)=0.0
          do 702 l=1,10
            sdc(l,ks)=0.0
            dc(l)=0.0
  702     continue
          do 703 kr=1,nrows
          do 703 kc=1,ncols
            inc(ks,kr,kc)=0
  703     continue
  701     continue

      do 70 kc=1,ncols
      do 71 kr=1,nrows
      if (ind(kr,kc).eq.0) go to 71
c     Individual tree loop ...
        do 72 ki=1,ind(kr,kc)
          ks=isp(ki,kr,kc)
          d=dbh(ki,kr,kc)
c     ... Increment diameter classes:
          l=(d/10.0)+1
          if (l.gt.10) l=10
          sdc(l,ks)=sdc(l,ks)+1.0
          dc(l)=dc(l)+1.0
c     ... and stem density:
          sd(ks)=sd(ks)+1.0
          td=td+1.0
          if (d.ge.10.0) tb=tb+1.0
c     ... and dbh tallies:
          tdbh=tdbh+d
          ssqd=ssqd+d**2
c     Set inc=1 if species is included on this plot:
          inc(ks,kr,kc)=1
   72     continue
c     Increment biomass and LAI (using a vertical LAI profile):
        tbiom=tbiom+biom(kr,kc)
        tla=tla+DLA(kr,kc,1,MH,0)
c     ... and basal area:
        do 711 ks=1,nspp
          ba(ks)=ba(ks)+sba(ks,kr,kc)
          tba=tba+sba(ks,kr,kc)
  711     continue
c     ... and canopy height:
        tht=tht+FLOAT(maxht(kr,kc))
   71   continue
   70   continue

      sumriv=0.0
c     Compute species frequency:
      do 741 ks=1,nspp
        do 742 kc=1,ncols
        do 742 kr=1,nrows
          if (inc(ks,kr,kc).ne.0) freq(ks)=freq(ks)+1.0
  742   continue
        freq(ks)=freq(ks)/plots
  741   continue

        do 74 ks=1,nspp
          if (td.gt.0.0) then
            rsd(ks)=100.0*sd(ks)/td
            else
            rsd(ks)=0.0
          endif
          if (tba.gt.0.0) then
            rba(ks)=100.0*ba(ks)/tba
            else
            rba(ks)=0.0
          endif        
          riv(ks)=(rsd(ks)+rba(ks))/2.0
          sumriv=sumriv+riv(ks)
   74     continue

c     Convert stand aggregates to per-ha ...
        tbiom=tbiom/1000.0
        tbmha=tbiom/xarea
        tdha=td/xarea
        tdb=tb/xarea
        tbaha=tba/xarea
        xlai=tla/plots
        xht=tht/plots
        if (td.gt.0.0) then 
          xdbh=tdbh/td
          sdd=SQRT((ssqd-tdbh**2/td)/(td-1.0))
          endif
c       add comma delimiter to facilitate parsing in C
        do 77 ks=1,nspp
         msp1(ks)= msp(ks)//','
  77    continue 
c      call C function to write in order to avoid IO in fortran
      call zprint(nspp,msp1,kyr,xarea,sd,sdc,dc,rsd,ba,rba,
     2            riv,freq, tdha, tdb, tbaha, xdbh, sdd, tbmha, 
     3            xlai, xht)
      
      return
      end


c 7.1 Uniform Deviates 273 Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF c SCIENTIFIC COMPUTING (ISBN c 0-521-43064-X)
c Copyright (C) 1986-1992 by Cambridge University Press. Programs Copyright (C) 1986-1992 by c c Numerical Recipes Software.
c Long period (> 2×1018) random number generator of L’Ecuyer with Bays-Durham shuffle
c and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
c of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
c alter idum between successive deviates in a sequence. RNMX should approximate the largest
c floating value that is less than 1.

c     RAN1 is dedicated to weather.


        REAL FUNCTION RAN1(idum)
        INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *  IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     *  IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER idum1,j,k,iv1(NTAB),iy1
c        SAVE iv1,iy1,idum1
        DATA idum1/123456789/, iv1/NTAB*0/, iy1/0/

c 	Initialize
        if (idum.le.0) then 
         idum=max(-idum,1)
c 	 prevent idum=0 
         idum1=idum
         do 11 j=NTAB+8,1,-1
c 	  Load the shuffle table (after 8 warm-ups).
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
           if (j.le.NTAB) iv1(j)=idum
  11     continue
        iy1=iv1(1)
        endif

c	Start here when not initializing.
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
c	Compute idum=mod(IA1*idum,IM1) without overflowsby Schrage’s method. 
        if (idum.lt.0) idum=idum+IM1 
        k=idum1/IQ2
        idum1=IA2*(idum1-k*IQ2)-k*IR2 
c	Compute idum1=mod(IA2*idum1,IM2) likewise.
        if (idum1.lt.0) idum1=idum1+IM2
        j=1+iy1/NDIV
c	Will be in the range 1:NTAB.
        iy1=iv1(j)-idum1
c	Here idum is shuffled, idum and idum1 are combined to generate output.
        iv1(j)=idum
        if(iy1.lt.1)iy1=iy1+IMM1
        RAN1=min(AM*iy1,RNMX)
c	Because users don’t expect endpoint values.
        return
        END

c 7.1 Uniform Deviates 273 Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF c SCIENTIFIC COMPUTING (ISBN c 0-521-43064-X)
c Copyright (C) 1986-1992 by Cambridge University Press. Programs Copyright (C) 1986-1992 by c c Numerical Recipes Software.
c Long period (> 2×1018) random number generator of L’Ecuyer with Bays-Durham shuffle
c and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
c of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
c alter idum between successive deviates in a sequence. RNMX should approximate the largest
c floating value that is less than 1.

c     RAN2 is dedicated to REGEN.


        REAL FUNCTION RAN2(idum)
        INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *  IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     *  IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER idum2,j,k,iv2(NTAB),iy2
c        SAVE iv2,iy2,idum2
        DATA idum2/123456789/, iv2/NTAB*0/, iy2/0/

c 	Initialize
        if (idum.le.0) then 
         idum=max(-idum,1)
c 	 prevent idum=0 
         idum2=idum
         do 11 j=NTAB+8,1,-1
c 	  Load the shuffle table (after 8 warm-ups).
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
           if (j.le.NTAB) iv2(j)=idum
  11      continue
        iy2=iv2(1)
        endif

c	Start here when not initializing.
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
c	Compute idum=mod(IA1*idum,IM1) without overflowsby Schrage’s method. 
        if (idum.lt.0) idum=idum+IM1 
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2 
c	Compute idum2=mod(IA2*idum2,IM2) likewise.
        if (idum2.lt.0) idum2=idum2+IM2
        j=1+iy2/NDIV
c	Will be in the range 1:NTAB.
        iy2=iv2(j)-idum2
c	Here idum is shuffled, idum and idum2 are combined to generate output.
        iv2(j)=idum
        if(iy2.lt.1)iy2=iy2+IMM1
        RAN2=min(AM*iy2,RNMX)
c	Because users don’t expect endpoint values.
        return
        END

c 7.1 Uniform Deviates 273 Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF c SCIENTIFIC COMPUTING (ISBN c 0-521-43064-X)
c Copyright (C) 1986-1992 by Cambridge University Press. Programs Copyright (C) 1986-1992 by c c Numerical Recipes Software.
c Long period (> 2×1018) random number generator of L’Ecuyer with Bays-Durham shuffle
c and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
c of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
c alter idum between successive deviates in a sequence. RNMX should approximate the largest
c floating value that is less than 1.

c     RAN3 is dedicated to MORTAL.


        REAL FUNCTION RAN3(idum)
        INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *  IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     *  IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER idum3,j,k,iv3(NTAB),iy3
        SAVE iv3,iy3,idum3
        DATA idum3/123456789/, iv3/NTAB*0/, iy3/0/

c 	Initialize
        if (idum.le.0) then 
         idum=max(-idum,1)
c 	 prevent idum=0 
         idum3=idum
         do 11 j=NTAB+8,1,-1
c 	  Load the shuffle table (after 8 warm-ups).
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
           if (j.le.NTAB) iv3(j)=idum
  11    continue
         iy3=iv3(1)
        endif

c	Start here when not initializing.
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
c	Compute idum=mod(IA1*idum,IM1) without overflowsby Schrage’s method. 
        if (idum.lt.0) idum=idum+IM1 
        k=idum3/IQ2
        idum3=IA2*(idum3-k*IQ2)-k*IR2 
c	Compute idum3=mod(IA2*idum3,IM2) likewise.
        if (idum3.lt.0) idum3=idum3+IM2
        j=1+iy3/NDIV
c	Will be in the range 1:NTAB.
        iy3=iv3(j)-idum3
c	Here idum is shuffled, idum and idum3 are combined to generate output.
        iv3(j)=idum
        if(iy3.lt.1)iy3=iy3+IMM1
        RAN3=min(AM*iy3,RNMX)
c	Because users don’t expect endpoint values.
        return
        END

c     REGEN sprouts stumps and plants saplings.
      Subroutine REGEN(kr,kc)
      include 'zcommon.h'
      parameter (MY=7)
c     Num, number to plant; Cohorts are seedlings; Probs, the
c       relative chance per species; XS and XSS, expected saplings
c       and stump sprouts.
      integer num(MS)
      real rf(MS), probs(MS), xss(MS) 
      real cohorts(MY,MS,MR,MC)
c  *  LYrs, lag years before seedlings are available for planting,
c       is now an array, keyed by shade tolerance class,
c       assuming more tolerant species start more slowly.
      integer lyrs(5)
      character*50 lfmt(5)
      integer ns, nss, nstot, nsstot
      data lyrs / 2,3,4,5,6 /

c     Initialize cohorts in year 1:
      if (kyr.eq.1) then
        do 510 ks=1,nspp
        do 511 l=1,MY
          cohorts(l,ks,kr,kc)=0.0
  511     continue
  510     continue
        endif
c     Zero tallies for saplings:
c     XStot is expected total saplings; xsstot, expected stump-sprouts;
c       NS is number of saplings to plant; nss, number of sprouts;
c       NPoss the number possible, nposs2, possible after sprouts;
c       NPlant, the total number planted.
c     BAtot is total basal area, which is accounted when available space
c       is computed:  a significant % of ground area may be BA.
        batot=0.0
        xstot=0.0
        xsstot=0.0
        ns=0
        nss=0
        nstot=0
        nsstot=0
        nposs=0
        nposs2=0
        nplant=0
      do 52 ks=1,nspp
c     Total basal area:
        batot=batot+sba(ks,kr,kc)
        lt=light(ks)
c     Available light regen factor, based on ground-level light ...
        alrf=ALF(lt,al(0,kr,kc))
c     Soil moisture factor is based on topsoil dry-days;
c     Below-ground factor is MIN(water, fertility):
        if (nutri(ks).ne.0) then
          sfrf=sff(nutri(ks))
          else
          sfrf=1.0
          endif
        bgf=MIN(sfrf,smf(1,ks))
c     Regen factor is light*BGF*temperature:
        rf(ks)=alrf*bgf*ddf(ks)
c     Saplings expected from seeds, filtered for lyrs ...
c     Work thru cohorts, voiding each if rf=0 this year:
        ly=lyrs(lt)
        do 521 l=(ly-1),1,-1
          if (rf(ks).eq.0.0) then
            cohorts(l+1,ks,kr,kc)=0.0
            else
            cohorts(l+1,ks,kr,kc)=cohorts(l,ks,kr,kc)
            endif
  521     continue
c     This year's seedling establishment:
        cohorts(1,ks,kr,kc)=seed(ks)*rf(ks)
c     Total saplings surviving from seed:
        xstot=xstot+cohorts(ly,ks,kr,kc)
c     Expected number of stump sprouts ... 
c     Healthy stumps sprout to their full potential;
c       stress-killed stumps are subject to the regen multiplier: 
        xss(ks)=FLOAT(nstmp(1,ks)*nsprt(ks)) +
     2    FLOAT(nstmp(2,ks)*nsprt(ks))*rf(ks) 
        xsstot=xsstot+xss(ks) 
        num(ks)=INT(xss(ks)+RAN2(idum(2)))
        nsstot=nsstot+num(ks) 
   52   continue

c     Stump sprouts first:  
c     Number of sprouts is MIN(available,possible);
c       this is 1 per sq.m, minus trees present and sq.m of basal area:
        mba=INT(batot)
        nposs=nmax-ind(kr,kc)-mba
        nss=MIN(nsstot,nposs) 
c     No room at all:
        if (nposs.eq.0) go to 55
c     If NSStot is OK, keep them and go to seedlings: 
        if (nsstot.le.nposs) go to 54 
c     Else, NSStot gt NPoss, so remove some ...
   53 nssxs=nsstot-nposs
c     Compute equal-probability array (only for species with sprouts) 
      spp=0.0
      do 530 ks=1,nspp
        if (xss(ks).gt.0.0) spp=spp+1.0
  530   continue
      do 531 ks=1,nspp
        if (xss(ks).gt.0.0) then
          probs(ks)=1.0/spp
          else
          probs(ks)=0.0
          endif
        if (ks.ge.2) probs(ks)=probs(ks)+probs(ks-1)
  531   continue
c     Remove sprouts probabilistically: 
      nssr=0
  532 y=RAN2(idum(2))
      ks=1
      if (y.le.probs(ks)) go to 535 
  533 do 534 ks=2,nspp
        if (y.gt.probs(ks-1).and.y.le.probs(ks)) go to 535
  534   continue
  535 if (num(ks).gt.0) then
        num(ks)=num(ks)-1 
        nssr=nssr+1 
        endif 
      if (nssr.lt.nssxs) go to 532
c     ... and now there's no room left for seedlings,
c       so go plant sprouts:
        go to 55  

c     Space left after sprouting, so, add seedlings:  
   54 nposs2=nposs-nss
      nstot=INT(xstot+RAN2(idum(2)))
c     Number of seedlings is MIN(available,possible)
        ns=MIN(nstot,nposs2)
c     If no seedlings, plant any sprouts:
        if (ns.eq.0) go to 55 
c     Else, probabilistically generate seedlings
  541 do 542 ks=1,nspp
        probs(ks)=cohorts(lyrs(light(ks)),ks,kr,kc)/xstot
        if (ks.ge.2) probs(ks)=probs(ks)+probs(ks-1)
  542   continue
      ins=0 
  543 y=RAN2(idum(2))
      ks=1
      if (y.le.probs(ks)) go to 546 
  544 do 545 ks=2,nspp
        if (y.gt.probs(ks-1).and.y.le.probs(ks)) go to 546
  545   continue
  546 num(ks)=num(ks)+1 
      ins=ins+1 
      if (ins.lt.ns) go to 543

   55 in=ind(kr,kc) 
      nplant=ns+nss

c     Plant trees and assign initial diameters
      do 56 ks=1,nspp
        if (seed(ks).eq.0.0) go to 56
        if (nplant.eq.0.or.num(ks).eq.0) go to 56
c     Plant trees and assign initial values:
        do 561 ki=1,num(ks)
          in=in+1
          dbh(in,kr,kc)=2.50+0.5*ABS(GAUSS2(idum(2)))
          isp(in,kr,kc)=ks
          ibc(in,kr,kc)=1
          nogro(in,kr,kc)=0
  561     continue
   56   continue
      ind(kr,kc)=in 
      return
      end

c     SMORT is for stress mortality
      Logical Function SMORT(nogro,idum)
      integer nogro, idum
      if (nogro.gt.1.and.RAN3(idum).lt.(0.369)) then
        SMORT=.true.
        else
        SMORT=.false.
        endif
      return
      end

c     SOLWAT simulates the soil water balance.
c     This is SOLWAT version 2.3. Modified by Acevedo to include runon.
      Subroutine SOLWAT(kr,kc)
      include 'zcommon.h'
c     Water is current content (cm); OWater, old water (last month's);
c     Pond is current ponding (cm water);
c     Wk's are critical water contents, =75% of FC (from Sellers);
c     Dtot is total depth of soil profile; 
c     FRFT=fine-root fraction for each layer, for trees 
c       (roots decrease linearly with depth);
c     FRFS, the same, to depth of rooting zone for seedlings DRZS;
c     Effective LAI for interception, ELAI, is computed in BOOKS.
c     the water balance uses NTS timesteps
c     within each month, to better approximate water fluxes;
c     NTS is set to 16, and used as a real number TS.
      integer nts
      real sm, smr, wsc, xc, xi, dtot, xddt, drzs, xdds, ddm, ts
      real water(ML), owater(ML), ddays(ML)
      real days(12)
      real pond(MR,MC)
c     introduced for the water ponding
      real deficit
      data days/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./
c     XC, interception loss rate (.%/LAI) (from Aston 1979);
      data ddbase /5.56/
      data wsc /0.08/
      data xc /0.05/
      data nts /16/
      data ts /16.0/

c     Initialize ponding at 0.0:
      pond(kr,kc)=0.0
c     soil type
      ksol=msol(kr,kc)
      nl=nsl(ksol)
c     Soil water balance ...
c       Use local variables for water balance:
      do 210 l=1,nl
        water(l)=sw(l,kr,kc)
        ddays(l)=0.0
  210   continue
c     Cday, calendar day; PETA (AETA), annual PET (AET);
c       Runoff, excess percolation (cumro, cumulative for year).
        cday=15.0
        peta=0.0
        aeta=0.0
        runoff=0.0
        cumro=0.0
        toti=0.0
c     Monthly loop:
      do 21 mo=1,12
        do 220 l=1,nl
        owater(l)=water(l)
  220   continue
      aet=0.0
c     Adjust for within-month timestep:
        rain=r(mo)/ts
        srunon=runon(mo)/ts
      do 22 kt=1,nts
c     pond dynamics:
      pond(kr,kc)=pond(kr,kc)+srunon

c     Priestley-Taylor PET ...
c     PET is decomposed into fractions driven by transmitted (TR) vs 
c       absorbed (AR) radiation; the former ET comes from the surface soil;
c       the latter comes from all layers based on their root density.
c     Total PET:
        if (t(mo).le.0.0) then
          p=0.0
          else
          hvap=597.391-0.568*t(mo)
          h=ct*(t(mo)-tx)*sun(mo)
          p=h/hvap*days(mo)/ts
          if (p.lt.0.0) p=0.0
          endif
        pet=p
        peta=peta+p

c     Interception loss, and adjust PET:
        xi=elai(mo,kr,kc)*xc*(rain)
        if (xi.gt.pet) xi=pet
        pet=pet-xi
        toti=toti+xi

c     Partition PET according to radiation transmitted to ground:
        tr=al(0,kr,kc)
        ar=1.0-tr
        pet1=tr*pet
        pet2=ar*pet

c     Water Input (wi) =  Rain - Interception + pond ...
c       Water is intercepted first; then fills layer 1, 
c         then percolates to layer 2 ...
         pond(kr,kc) = pond(kr,kc) + rain - xi
         if(pond(kr,kc).le.0.0) pond(kr,kc) = 0.0

c     AET and soil water drawdown ...
c       Demand varies for each layer:  
c       For first layer, demand is PET1 + its share of PET2;
c         for each lower layer, demand is its share of PET2.
c       Water percolates only if PET is met for the layer.
c       Drawdown is linear.
c     ... First layer:
c     check water first layer for infilt
          deficit = 0.0
          do 23 l=1,nl
           deficit = deficit + (fc(l,ksol) - water (l))
   23     continue
          if(deficit.lt.0.0) then
           wi=0.0
          else
           wi = MIN(pond(kr,kc),deficit)
          endif
           pond(kr,kc) = pond(kr,kc) - wi
          water(1)=water(1)+wi
          aet1=pet1+pet2*frft(1,ksol)
          if (water(1).lt.wk(1,ksol)) aet1=aet1*water(1)/wk(1,ksol)
          water(1)=water(1)-aet1
          aet=aet+aet1
          if (water(1).lt.0.0) water(1)=0.0
          perc=0.0
          if (water(1).gt.fc(1,ksol)) then
            perc=water(1)-fc(1,ksol)
            water(1)=fc(1,ksol)
            endif
c     ... Lower layers:  
        do 222 l=2,nl
          water(l)=water(l)+perc
          aetx=pet2*frft(l,ksol)
          if (water(l).lt.wk(l,ksol)) aetx=aetx*water(l)/wk(l,ksol)
          water(l)=water(l)-aetx
          aet=aet+aetx
          if (water(l).lt.0.0) water(l)=0.0
          perc=0.0
          if (water(l).gt.fc(l,ksol)) then
            perc=water(l)-fc(l,ksol)
            water(l)=fc(l,ksol)
            endif

  222     continue
   22   continue
        aeta=aeta+aet+xi
c     ... Adjust runoff:
        runoff=perc
        cumro=cumro+runoff

c     Interpolate drought-days between months as necessary
        ocday=cday
        cday=cday+days(mo)
c     If in growing season, adjust dry-days:
        if (cday.gt.bgs.and.ocday.lt.egs) then
c     Do each layer ...
        do 224 l=1,nl
c     Else, interpolate ...
          w=water(l)
          ow=owater(l)
          wpl=wp(l,ksol)
c     Sufficient soil water both months:
          if (ow.gt.wpl.and.w.gt.wpl) ddm=0.0
c     Dry both months:
          if (ow.le.wpl.and.w.le.wpl) then
            ddm=days(mo)
            if (ocday.lt.bgs.and.cday.gt.bgs) ddm=cday-bgs
            if (ocday.lt.egs.and.cday.gt.egs) ddm=egs-ocday
            endif
c     Drying period:
          if (ow.gt.wpl.and.w.le.wpl) then
            ddm=days(mo)*(wpl-w)/(ow-w)
            if (ocday.lt.bgs.and.cday.gt.bgs) ddm=MIN(ddm,cday-bgs)
            if (ocday.lt.egs.and.cday.gt.egs) ddm=egs-cday+ddm
            if (ddm.lt.0.0) ddm=0.0
            endif
c     Wetting period:
          if (ow.le.wpl.and.w.gt.wpl) then
            ddm=days(mo)*(wpl-ow)/(w-ow)
            if (ocday.lt.bgs.and.cday.gt.bgs) ddm=ocday+ddm-bgs
            if (ddm.lt.0.0) ddm=0.0
            if (ocday.lt.egs.and.cday.gt.egs) ddm=MIN(ddm,egs-ocday)
            endif
c     ... and tally up dry-days this month:
          ddays(l)=ddays(l)+ddm
  224     continue
          endif
   21 continue

c     Relativize dry-days to growing-season length, and average 
c       (XDDt is integrated by FRFT over the whole profile; 
c       the topsoil XDDs is integrated over FRFS):
        xdds=0.0
        xddt=0.0
      do 225 l=1,nl
        ddays(l)=ddays(l)/tgs
        xddt=xddt+ddays(l)*frft(l,ksol)
        xdds=xdds+ddays(l)*frfs(l,ksol)
c     ... and reassign local water status to COMMON:
        sw(l,kr,kc)=water(l)
  225   continue
c     Compute species soil moisture factors:
      do 226 ksp=1,nspp
        smf(1,ksp)=DRTF(mdrt(ksp),xdds)
        smf(2,ksp)=DRTF(mdrt(ksp),xddt)
  226   continue

      return
      end

c     TLINE does linear interpolations on temperatures,
c       to find bgs and egs, given monthly temps and the gdd base
      Real Function TLINE(m1,t1,m2,t2,ddb)
      real d1, t1, d2, t2, ddb
      integer m1, m2
      real cd(12)
      data cd/15.0,46.0,74.0,105.0,135.0,166.0,196.0,227.0,258.0,
     2  288.0,319.0,349.0/
        d1=cd(m1)
        d2=cd(m2)
        if (m1.eq.12) d1=d1-365.0
c     Regression slope is a; intercept is b:
        a=(t2-t1)/(d2-d1)
        b=t2-a*d2
c     Solve for calendar day at which t=ddbase:
      d=(ddb-b)/a
      TLINE=d
      return
      end

      Real Function WOOD(dbh,kf)
      real dbh
      integer kf
      real c1(9), c2(9)
      data c1 / -3.4053, -3.4032, -1.7383, -2.5095, -1.8817, 
     &  -2.7894, -3.0703, -3.3763, -2.7894 /
      data c2 /2.6933, 2.6031, 2.3320, 2.5638, 2.3293,
     &  2.6273, 2.6870, 2.7242, 2.6273 /
c     Form labels: 1=Abies, 2=Pinus, 3=PSme, 4=Tsuga, 5=Picea; 
c       6=Acer, 7=Alnus, 8=Quercus, 9=average or default.
c     ABspp, PIspp, PSme, TShe, PIsi, 
c       ACma, ALru, CAch, and ACma (all from BioLib).
c     WOOD returns BOLE+BARK+BRANCH biomass, estimated by regressing
c       this sum compiled from separate allometries.
      WOOD=EXP(c1(kf)+c2(kf)*LOG(dbh))
      return
      end
