      subroutine cerio

c     global prefixed parameters number of stages, delay order
      parameter(mt=4,mo=3000)

c     global, simulation control parameters
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv

c     global, model variables
      common/recru1/ram(mo,mt),rm(mo,mt),raf(mo,mt),rf(mo,mt)
      common/recru2/recrum(mt),recruf(mt)
      common/growth/g(mo,mt),ga(mo,mt)
      common/compu1/tau(mt),eff(mt),deltau(mt),s(mt),delmax(mt)
      common/compu2/growth(mt),consu(mt),grorat(mt),delold(mt)
      common/compu3/delay(mt),starv(mt),temp,phper,slight
      common/compu4/gsligh,gxkd,gext,sratio,slat,clat
      common/compu5/tharv(mt),ampmu(mt),amu(mt),sm(mo,mt),sa(mo,mt)
      common/morta2/smu(mt)
      common/compu6/rematf(mt),rematm(mt),resexf(mt),resexm(mt)
      common/state1/pnum(mt),pnumm(mt),pnumf(mt),totnum,eggs,dorm
      common/state2/food1,food2,foodt,taccu,paccu
      common/state/yd(10)
      common/index/x_max(10),x_min(10),x_avg(10),x_rms(10),x_end(10)

c     multiple runs, variables for control of sensitivity analysis
      type sens
         character*17 name
         real max,min,step
         integer irun,k_flag
      end type
      type(sens) sen

c     model parameters by members of array
      type par_number
      sequence        
          real para(11*mt+54)
      end type

c     model parameters by name
      type par_name
      sequence        
          real    chlong(mt),taumin(mt),taumax(mt),scalsx
          real    cmax(mt),pkhalf(mt),tcmax,tcmin,pmax,pmin
          real    stcri(mt),stsen(mt)
          real    fertmx,broodx,fercri,fersen,bdcri,bdsen,hmh,fmale
          real    swsex,swsen
          real    pmu(mt),pmegg,ampmax(mt),thresh(mt),restre(mt)
          real    tini,trest,thold,phini,phold
c     model parameters: food and environmental
          real    fsmax1,fsmax2,fmax1,fmax2,tf1,tf2,plgs,dlay12
          real    qual1,qual2, fgmax1,fgmax2
          real    tmed,tamp,slmed,slamp,tdel,sldel,plat,cphper,tox
          real    skcch,pgccel,slopt,skext,eupho
          real    tmax1,tmax2,tmin1,tmin2
c     initial conditions
          real    food1i,food2i,pinocul
      end type

      type(par_number) eqn
      type(par_name) equ
      equivalence(eqn,equ)

c     variables for control of multiple runs
      type simul
          integer nruns,ndata,nstate
      end type
      type(simul) sim

      character*17 sub_title
     
c     hard coded number of states
      mst=mt-1
      sim%nstate=mst+7

c     files
      open(4,file='cerio_inp.txt')

c     default value for sensitivity flag
      sen%k_flag=1

c     reading controls and parameters
      call readin(sen%name,sen%max,sen%min,sen%step,sub_title)
      call readpar(equ,sen%name)

c     calculating number of runs and data
      call num_round(sen%max,sen%min,sen%step,sim%nruns,sim%ndata)

c     initialization in multiple runs: first param value, run#, modify params
      eqn%para(sen%k_flag)=sen%min

      sen%irun=1

 10   call calcpar(equ)

c     initialization of one run
      t=time0
      ict=0
      call initial(equ)
      call assign_state
      call calc_index(sim%nstate)
  1   call cerioprint(sen%irun, sen%name, eqn%para(sen%k_flag),t,time0,
     &                sim%nruns, sim%nstate, sim%ndata, sub_title,yd)


c     calculations
  2   call calcu(equ)
      call update
      call assign_state
      call calc_index(sim%nstate)

c     iterate within a run
      t=t+dt
      ict=ict+1
      if(ict.lt.mdt) goto 2
      ict=0
      abt=abs(t-tfin)
      if (t.lt.tfin.or.abt.lt.0.001) goto 1

c     end of a run
      call cerioindex(sen%irun, sen%name, eqn%para(sen%k_flag),t,time0,
     &                sim%nruns, sim%nstate, sim%ndata, sub_title,yd,
     &                x_max,x_min,x_avg,x_rms,x_end)

c     iterate parameter value for next run
      eqn%para(sen%k_flag)=eqn%para(sen%k_flag)+sen%step
      sen%irun=sen%irun+1
      if(sen%irun.le.sim%nruns) goto 10

c     end of last run, close files and end program
      close(4)
      end

      subroutine readin(senname,senmax,senmin,senstep,sub_title)
c     reads the simulation control parameters
      parameter(mt=4,mo=3000)

c     simulation control
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv

      type sens 
         character*17 name
         real max,min,step
         integer irun,k_flag
      end type
      type(sens) sen

      character*17 label
      character*17 time_unit,sub_title
      character*17 senname
      real f_mdt, senmax,senmin,senstep

      read(4,3) label, senname
      read(4,1) label, senmax
      read(4,1) label, senmin
      read(4,1) label, senstep
      read(4,3) label, sub_title
      read(4,1) label, time0
      read(4,1) label, tfin
      read(4,1) label, dt
      read(4,1) label, f_mdt
      mdt = aint(f_mdt)
      read(4,3) label, time_unit
      read(4,2) label,lab
      read(4,2) label,lfv

 1    format(a17, f17.0)
 2    format(a17, i17)
 3    format(a17, a17)

      sen%name = senname
      sen%max = senmax
      sen%min = senmin
      sen%step = senstep

      end

      subroutine num_round(senmax,senmin,senstep,simnruns,simndata)
c     calculates number of runs, data points

c     simulation control
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv

      type sens
         character*17 name
         real max,min,step
         integer irun,k_flag
      end type
      type(sens) sen

      type simul
          integer nruns,ndata,nstate
      end type
      type(simul) sim

      real senmax, senmin,senstep
      integer simnruns, simndata

      if (senstep.ne.0.00) then
          runs=(senmax-senmin)/senstep
      else
          runs=0.00
      endif

      if((runs-aint(runs)).ge.0.5) then
          simnruns=anint(runs)+1
        else
          simnruns=aint(runs)+1
      endif

      data=(tfin-time0)/(dt*float(mdt))
      if((data-aint(data)).ge.0.5) then
          simndata=anint(data)+1
        else
          simndata=aint(data)+1
      endif
      sim%nruns = simnruns
      sim%ndata = simndata  

      end

      subroutine readpar(equ,senname)
c     reads the parameters and check sensitivity flag, serviced by f_sen
      parameter(mt=4,mo=3000)

c     simulation control
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv

c     shared with f_sen
      common/counter/i_cnt_par

      type sens
         character*17 name
         real max,min,step
         integer irun,k_flag
      end type
      type(sens) sen

c     model parameters by name       
      type par_name
      sequence        
          real    chlong(mt),taumin(mt),taumax(mt),scalsx
          real    cmax(mt),pkhalf(mt),tcmax,tcmin,pmax,pmin
          real    stcri(mt),stsen(mt)
          real    fertmx,broodx,fercri,fersen,bdcri,bdsen,hmh,fmale
          real    swsex,swsen
          real    pmu(mt),pmegg,ampmax(mt),thresh(mt),restre(mt)
          real    tini,trest,thold,phini,phold
c     model parameters: food and environmental
          real    fsmax1,fsmax2,fmax1,fmax2,tf1,tf2,plgs,dlay12
          real    qual1,qual2, fgmax1,fgmax2
          real    tmed,tamp,slmed,slamp,tdel,sldel,plat,cphper,tox
          real    skcch,pgccel,slopt,skext,eupho
          real    tmax1,tmax2,tmin1,tmin2
c     initial conditions
          real    food1i,food2i,pinocul
      end type
      type(par_name) equ

c     local
      character*17 label
      character*17 senname

      i_cnt_par = 1

      do 10 j=1,mt
      call f_sen(j,equ%chlong(j), senname)
 10   continue
      read(4,2) label, unit_1

      do 20 j=1,mt
      call f_sen(j,equ%taumin(j), senname)
 20   continue
      do 30 j=1,mt
      call f_sen(j,equ%taumax(j), senname)
 30   continue
      read(4,2) label, unit_2

      call f_sen(1,equ%scalsx, senname)
      read(4,2) label, unit_3

      do 40 j=1,mt
      call f_sen(j,equ%cmax(j), senname)
 40   continue
      read(4,2) label, unit_4

      do 50 j=1,mt
      call f_sen(j,equ%pkhalf(j), senname)
 50   continue
      read(4,2) label, unit_5

      call f_sen(1,equ%tcmax, senname)
      call f_sen(1,equ%tcmin, senname)
      read(4,2) label, unit_6

      call f_sen(1,equ%pmax, senname)
      call f_sen(1,equ%pmin, senname)
      read(4,2) label, unit_7

      do 60 j=1,mt
      call f_sen(j,equ%stcri(j), senname)
 60   continue
      do 70 j=1,mt
      call f_sen(j,equ%stsen(j), senname)
 70   continue
      read(4,2) label, unit_8

      call f_sen(1,equ%fertmx, senname)
      call f_sen(1,equ%broodx, senname)
      read(4,2) label, unit_9

      call f_sen(1,equ%fercri, senname)
      call f_sen(1,equ%fersen, senname)
      call f_sen(1,equ%bdcri, senname)
      call f_sen(1,equ%bdsen, senname)
      call f_sen(1,equ%hmh, senname)
      call f_sen(1,equ%fmale, senname)
      call f_sen(1,equ%swsex, senname)
      call f_sen(1,equ%swsen, senname)
      read(4,2) label, unit_10

      do 80 j=1,mt
      call f_sen(j,equ%pmu(j), senname)
 80   continue
      call f_sen(1,equ%pmegg, senname)
      read(4,2) label, unit_11

      do 90 j=1,mt
      call f_sen(j,equ%ampmax(j), senname)
 90   continue
      read(4,2) label, unit_12

      do 100 j=1,mt
      call f_sen(j,equ%thresh(j), senname)
 100  continue
      do 110 j=1,mt
      call f_sen(j,equ%restre(j), senname)
 110  continue
      read(4,2) label, unit_13

      call f_sen(1,equ%tini, senname)
      call f_sen(1,equ%trest, senname)
      call f_sen(1,equ%thold, senname)
      read(4,2) label, unit_14

      call f_sen(1,equ%phini, senname)
      call f_sen(1,equ%phold, senname)
      read(4,2) label, unit_15

      if(lab.eq.1) then
          call f_sen(1,equ%fsmax1, senname)
          call f_sen(1,equ%fsmax2, senname)
          read(4,2) label, unit_16
          call f_sen(1,equ%fmax1, senname)
          call f_sen(1,equ%fmax2, senname)
          read(4,2) label, unit_17
              if(lfv.eq.1) then
                  call f_sen(1,equ%tf1, senname)
                  call f_sen(1,equ%tf2, senname)
                  read(4,2) label, unit_18
                  call f_sen(1,equ%plgs, senname)
                  call f_sen(1,equ%dlay12, senname)
                  read(4,2) label, unit_19
              else
                  i_cnt_par=i_cnt_par+4
              endif
       else
          i_cnt_par=i_cnt_par+4
      endif

      call f_sen(1,equ%qual1, senname)
      call f_sen(1,equ%qual2, senname)
      read(4,2) label, unit_20

      if(lab.eq.0) then
          call f_sen(1,equ%fgmax1, senname)
          call f_sen(1,equ%fgmax2, senname)
          read(4,2) label, unit_21
      else
          i_cnt_par=i_cnt_par+2
      endif

      call f_sen(1,equ%tmed, senname)

      if(lab.eq.1) then
          read(4,2) label, unit_211
          i_cnt_par=i_cnt_par+5
      else
          call f_sen(1,equ%tamp, senname)
          read(4,2) label, unit_22
          call f_sen(1,equ%slmed, senname)
          call f_sen(1,equ%slamp, senname)
          read(4,2) label, unit_23
          call f_sen(1,equ%tdel, senname)
          call f_sen(1,equ%sldel, senname)
          read(4,2) label, unit_24
      endif

      if(lab.eq.0) then
          call f_sen(1,equ%plat, senname)
          read(4,2) label, unit_25
          i_cnt_par=i_cnt_par+1
      else
          call f_sen(1,equ%cphper, senname)
          read(4,2) label, unit_26
          i_cnt_par=i_cnt_par+1
      endif

      call f_sen(1,equ%tox, senname)
      read(4,2) label, unit_27

      if(lab.eq.0) then
          call f_sen(1,equ%skcch, senname)
          read(4,2) label, unit_28
          call f_sen(1,equ%pgccel, senname)
          read(4,2) label, unit_29
          call f_sen(1,equ%slopt, senname)
          read(4,2) label, unit_30
          call f_sen(1,equ%skext, senname)
          read(4,2) label, unit_31
          call f_sen(1,equ%eupho, senname)
          read(4,2) label, unit_32
          call f_sen(1,equ%tmax1, senname)
          call f_sen(1,equ%tmin1, senname)
          call f_sen(1,equ%tmax2, senname)
          call f_sen(1,equ%tmin2, senname)
          read(4,2) label, unit_33
      else
         i_cnt_par=i_cnt_par+9
      endif

      call f_sen(1,equ%food1i, senname)
      call f_sen(1,equ%food2i, senname)
      read(4,2) label, unit_34

      call f_sen(1,equ%pinocul, senname)
      read(4,2) label, unit_35

 2    format (a17,a17)
      end

      subroutine f_sen(j, par_value, senname)
c     services readpar

c     simulation control
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv

c     shared with readpar
      common/counter/i_cnt_par

      type sens
         character*17 name
         real max,min,step
         integer irun,k_flag
      end type
      type(sens) sen
c     local
      character*17 label
      character*17 senname

      if(j.le.mst) then
          read(4,1) label, par_value
 1        format (a17,f17.0)
          if(label.eq.senname) sen%k_flag=i_cnt_par
      endif

      i_cnt_par = i_cnt_par + 1
      return
      end

      subroutine calcpar(equ)
c     calculates some parameters from the values read from input file
      parameter(mt=4,mo=3000)

c     simulation control
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv

c     model variables
      common/compu1/tau(mt),eff(mt),deltau(mt),s(mt),delmax(mt)
      common/compu4/gsligh,gxkd,gext,sratio,slat,clat
      common/compu5/tharv(mt),ampmu(mt),amu(mt),sm(mo,mt),sa(mo,mt)
      common/compu6/rematf(mt),rematm(mt),resexf(mt),resexm(mt)

c     model parameters by name       
      type par_name
      sequence        
          real    chlong(mt),taumin(mt),taumax(mt),scalsx
          real    cmax(mt),pkhalf(mt),tcmax,tcmin,pmax,pmin
          real    stcri(mt),stsen(mt)
          real    fertmx,broodx,fercri,fersen,bdcri,bdsen,hmh,fmale
          real    swsex,swsen
          real    pmu(mt),pmegg,ampmax(mt),thresh(mt),restre(mt)
          real    tini,trest,thold,phini,phold
c     model parameters: food and environmental
          real    fsmax1,fsmax2,fmax1,fmax2,tf1,tf2,plgs,dlay12
          real    qual1,qual2, fgmax1,fgmax2
          real    tmed,tamp,slmed,slamp,tdel,sldel,plat,cphper,tox
          real    skcch,pgccel,slopt,skext,eupho
          real    tmax1,tmax2,tmin1,tmin2
c     initial conditions
          real    food1i,food2i,pinocul
      end type
      type(par_name) equ

      do 1 j=2,mst
          eff(j)=equ%chlong(j)/(equ%taumin(j)*equ%cmax(j))
          equ%thresh(j)=equ%thresh(j)*eff(j)*equ%cmax(j)
          equ%restre(j)=equ%restre(j)*eff(j)*equ%cmax(j)
          tharv(j)=equ%swsex*eff(j)*equ%cmax(j)
  1   continue
      equ%fercri=equ%fercri*equ%cmax(mst)*equ%swsex
      equ%bdcri=equ%bdcri*equ%cmax(mst)*equ%swsex
      equ%fersen=equ%fersen*equ%cmax(mst)
      equ%bdsen=equ%bdsen*equ%cmax(mst)
      slat=sin((3.14/180.0)*equ%plat)
      clat=cos((3.14/180.0)*equ%plat)

      end

      subroutine  initial(equ)
      parameter(mt=4,mo=3000)
c     simulation control
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv
c     model variables
      common/recru1/ram(mo,mt),rm(mo,mt),raf(mo,mt),rf(mo,mt)
      common/growth/g(mo,mt),ga(mo,mt)
      common/compu1/tau(mt),eff(mt),deltau(mt),s(mt),delmax(mt)
      common/compu2/growth(mt),consu(mt),grorat(mt),delold(mt)
      common/compu3/delay(mt),starv(mt),temp,phper,slight
      common/compu4/gsligh,gxkd,gext,sratio,slat,clat
      common/compu5/tharv(mt),ampmu(mt),amu(mt),sm(mo,mt),sa(mo,mt)
      common/state1/pnum(mt),pnumm(mt),pnumf(mt),totnum,eggs,dorm
      common/state2/food1,food2,foodt,taccu,paccu

c     model parameters by name       
      type par_name
      sequence        
          real    chlong(mt),taumin(mt),taumax(mt),scalsx
          real    cmax(mt),pkhalf(mt),tcmax,tcmin,pmax,pmin
          real    stcri(mt),stsen(mt)
          real    fertmx,broodx,fercri,fersen,bdcri,bdsen,hmh,fmale
          real    swsex,swsen
          real    pmu(mt),pmegg,ampmax(mt),thresh(mt),restre(mt)
          real    tini,trest,thold,phini,phold
c     model parameters: food and environmental
          real    fsmax1,fsmax2,fmax1,fmax2,tf1,tf2,plgs,dlay12
          real    qual1,qual2, fgmax1,fgmax2
          real    tmed,tamp,slmed,slamp,tdel,sldel,plat,cphper,tox
          real    skcch,pgccel,slopt,skext,eupho
          real    tmax1,tmax2,tmin1,tmin2
c     initial conditions
          real    food1i,food2i,pinocul
      end type
      type(par_name) equ

      do 1 j=1,mst
          pnumm(j)=0.0
          pnumf(j)=0.0
          pnum(j)=pnumm(j)+pnumf(j)
          tau(j)=equ%taumin(j)
          deltau(j)=tau(j)/dt
          delmax(j)=equ%taumax(j)/dt
          idmax=int(delmax(j))
          do 2 i=1,idmax
              g(i,j)=equ%chlong(j)/tau(j)
              rf(i,j)=0.0
              rm(i,j)=0.0
              sm(i,j)=equ%pmu(j)
 2        continue
          starv(j)=0.0
          grorat(j)=1.0
          s(j)=exp(-equ%pmu(j)*tau(j))
 1    continue
      if (lab.eq.0) then
          soldec=0.4093*sin((6.28/365.0)*(t-82.2))
          soldlv=(-slat*sin(soldec)-0.1047)/(clat*cos(soldec))
          phper=7.639*acos(soldlv)
          temp=equ%tmed+equ%tamp*sin((6.28/365.0)*(t-equ%tdel))
          slight=equ%slmed+equ%slamp*sin((6.28/365.0)*(t-equ%sldel))
      else
          phper=equ%cphper
          temp=equ%tmed
          slight=0.0
      endif
      taccu=0.0
      paccu=0.0
      food1=equ%food1i
      food2=equ%food2i
      foodt=food1+food2
      dorm=equ%pinocul
      eggs=0.0

      delay(1)=deltau(1)
      grorat(1)=1.00
      do 3 j=2,mst
          delay(j)=delay(j-1)+deltau(j)
          delold(j)=deltau(j)-(1.0-grorat(j))
          if(delold(j).gt.deltau(j)) then
              idta=int(delold(j))
          else
              idta=int(deltau(j))
          endif
      gtest=g(idta,j)
      if(gtest.eq.0.0) then
          grorat(j)=1.00
      else
          grorat(j)=g(1,j)/gtest
      endif
      if(pnum(j).eq.0.0) then
          grorat(j)=1.00
      endif
      if(grorat(j).gt.2.0) then
          grorat(j)=2.00
      endif
  3   continue
      end

      subroutine assign_state
      parameter(mt=4,mo=3000)
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv
      common/state1/pnum(mt),pnumm(mt),pnumf(mt),totnum,eggs,dorm
      common/state2/food1,food2,foodt,taccu,paccu
      common/state/yd(10)
      common/compu3/delay(mt),starv(mt),temp,phper,slight

      dimension propm(mt),propf(mt),prop(mt)

      totnum=0.0
      do 1 j=1,mst
      totnum=totnum+pnum(j)+eggs
  1   continue

c     ------remove the branching to calculate the proportions --
c     ------assign proportions to state as an option-------------
      goto 3
      do 2 j=1,mst
      if(totnum.eq.0.0) then
          propf(j)=0.0
          propm(j)=0.0
      else
          propf(j)=pnumf(j)/totnum
          propm(j)=pnumm(j)/totnum
      endif
      prop(j)=propm(j)+propf(j)
  2   continue
      if(totnum.eq.0.0) then
          propeg=0.0
      else
          propeg=eggs/totnum
      endif
c ----end of calculation of proportions---------------------

 3    yd(1)=eggs
      yd(2)=pnumf(2)
      yd(3)=pnumf(3)
      yd(4)=pnumm(2)+pnumm(3)
      yd(5)=food1
      yd(6)=food2
      yd(7)=totnum
      yd(8)=temp
      yd(9)=phper
      yd(10)=slight
      end

      subroutine calcu(equ)
      parameter(mt=4,mo=3000)
c     simulation control
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv
c     model variables
      common/recru1/ram(mo,mt),rm(mo,mt),raf(mo,mt),rf(mo,mt)
      common/recru2/recrum(mt),recruf(mt)
      common/growth/g(mo,mt),ga(mo,mt)
      common/compu1/tau(mt),eff(mt),deltau(mt),s(mt),delmax(mt)
      common/compu2/growth(mt),consu(mt),grorat(mt),delold(mt)
      common/compu3/delay(mt),starv(mt),temp,phper,slight
      common/compu4/gsligh,gxkd,gext,sratio,slat,clat
      common/compu5/tharv(mt),ampmu(mt),amu(mt),sm(mo,mt),sa(mo,mt)
      common/morta2/smu(mt)
      common/compu6/rematf(mt),rematm(mt),resexf(mt),resexm(mt)
      common/state1/pnum(mt),pnumm(mt),pnumf(mt),totnum,eggs,dorm
      common/state2/food1,food2,foodt,taccu,paccu

c     model parameters by name       
      type par_name
      sequence        
          real    chlong(mt),taumin(mt),taumax(mt),scalsx
          real    cmax(mt),pkhalf(mt),tcmax,tcmin,pmax,pmin
          real    stcri(mt),stsen(mt)
          real    fertmx,broodx,fercri,fersen,bdcri,bdsen,hmh,fmale
          real    swsex,swsen
          real    pmu(mt),pmegg,ampmax(mt),thresh(mt),restre(mt)
          real    tini,trest,thold,phini,phold
c     model parameters: food and environmental
          real    fsmax1,fsmax2,fmax1,fmax2,tf1,tf2,plgs,dlay12
          real    qual1,qual2, fgmax1,fgmax2
          real    tmed,tamp,slmed,slamp,tdel,sldel,plat,cphper,tox
          real    skcch,pgccel,slopt,skext,eupho
          real    tmax1,tmax2,tmin1,tmin2
c     initial conditions
          real    food1i,food2i,pinocul
      end type
      type(par_name) equ

c     temperature, light and photoperiod calculations
      if (lab.eq.0) then
          soldec=0.4093*sin((6.28/365.0)*(t-82.2))
          soldlv=(-slat*sin(soldec)-0.1047)/(clat*cos(soldec))
          phper=7.639*acos(soldlv)
          temp=equ%tmed+equ%tamp*sin((6.28/365.0)*(t-equ%tdel))
          slight=equ%slmed+equ%slamp*sin((6.28/365.0)*(t-equ%sldel))
      else
          phper=equ%cphper
          temp=equ%tmed
          slight=0.0
      endif

c     algae growth multipliers for field conditions
      if (lab.eq.0) then
          fgrtm1=4.0*(temp-equ%tmin1)*(equ%tmax1-temp)
          fgrtm1=fgrtm1/(equ%tmax1-equ%tmin1)**2
          fgrtm2=4.0*(temp-equ%tmin2)*(equ%tmax2-temp)
          fgrtm2=fgrtm2/(equ%tmax2-equ%tmin2)**2
      endif

c     next season for field conditions
      if(lab.eq.0) then
          time01=abs(time0+365.0-t)
          if (time01.le.0.1) then
              food1=equ%food1i
              food2=equ%food2i
              foodt=food1+food2
          endif
      endif

c     ------ calculation of food supply rate

c     for field conditions
      if(lab.eq.0) then
          foodc=foodt*equ%pgccel
          gext=0.0088*foodc/equ%skcch
          gext=gext + 0.053*(foodc/equ%skcch)**(2.0/3.0)+equ%skext
          sratio=slight/equ%slopt
          gexkd=gext*equ%eupho
          gsligh = (exp(-sratio * exp(-gexkd))-exp(-sratio))
          gsligh=gsligh  * (exp(1.0)/gexkd)
          fsupl1=equ%fgmax1*fgrtm1*gsligh*food1
          fsupl2=equ%fgmax2*fgrtm2*gsligh*food2

c     for laboratory conditions
      else

c     for variable staggered food supply
      if(lfv.eq.1) then
             fsupl1=equ%fsmax1*sin(3.14*t/(equ%tf1*equ%plgs))
             fsupl2=sin(3.14*(t-equ%dlay12)/(equ%tf2*equ%plgs))
             fsupl2=equ%fsmax2*fsupl2

c     for constant supply
      else
              fsupl1=equ%fsmax1
              fsupl2=equ%fsmax2
          endif
      endif

c     make sure to exclude negative supply
      if(fsupl1.lt.0.0) then
          fsupl1=0.0
      endif
      if(fsupl2.lt.0.0) then
          fsupl2=0.0
      endif
c     -------- end of calculation of food supply

c     calculating total food and ratio of food items
      foodt=food1+food2
      if(foodt.eq.0.0) then
          ratio=0.0
      else
          ratio=food1/foodt
      endif

c     calculating growth limiting factors (multipliers)

      do 1 j=2,mst
       flim=foodt/(equ%pkhalf(j)+foodt)
       tlim=4.0*(temp-equ%tcmin)*(equ%tcmax-temp)
       tlim=tlim/(equ%tcmax-equ%tcmin)**2
       plim=4.0*(phper-equ%pmin)*(equ%pmax-phper)
       plim=plim/(equ%pmax-equ%pmin)**2
       slim= 1.0/(1.0 + exp((equ%tox-equ%stcri(j))/equ%stsen(j)))
       consu(j)=equ%cmax(j)*flim*tlim*plim*slim
       growth(j)=consu(j)*eff(j)*(ratio*equ%qual1+(1.0-ratio)*equ%qual2)
 1    continue

c     proceed to calculate the switching function
      quot=consu(mst)/equ%cmax(mst)
      asex=1.0/(1.0 + exp((equ%swsex-quot)/equ%swsen))

c     now mortalities will be calculated taking into account potential
c     increase due to starvation

c     for eggs there is no potential increase
      amu(1)=equ%pmu(1)
c     but for other classes there is an adjustment to be made
      do 2 j=2,mst
          ampmu(j)=(1.0+exp((equ%thresh(j)-starv(j))/equ%restre(j)))
          ampmu(j)=equ%ampmax(j)/ampmu(j)
          amu(j)=equ%pmu(j)*(1.0+ampmu(j))
 2    continue

c     ----- recruitments

      if(taccu.ge.equ%thold.and.paccu.ge.equ%phold) then
          hatch=1.0
      else
          hatch=0.0
      endif
      ferti=(1.0+exp((-consu(mst)+equ%fercri)/equ%fersen))
      ferti=equ%fertmx/ferti
      brood=equ%broodx/(1.0+exp((-consu(mst)+equ%bdcri)/equ%bdsen))
      if(pnumf(mst).ne.0.0) then
          sexrat=(pnumm(mst)+pnumm(2))/pnum(mst)
      else
          sexrat=0.0
      endif
      emate= sexrat/(equ%hmh+sexrat)
      do 3 j=1,mst
          deltau(j)=tau(j)/dt
 3    continue
      delay(1)=deltau(1)
      grorat(1)=1.00
      do 4 j=2,mst
          delay(j)=delay(j-1)+deltau(j)
          delold(j)=deltau(j)-(1.0-grorat(j))
          if(delold(j).gt.deltau(j)) then
              idta=int(delold(j))
          else
              idta=int(deltau(j))
          endif
          gtest=g(idta,j)
          if(gtest.eq.0.0) then
              grorat(j)=1.00
          else
              grorat(j)=g(1,j)/gtest
          endif
          if(pnum(j).eq.0.0) then
              grorat(j)=1.00
          endif
          if(grorat(j).gt.2.0) then
              grorat(j)=2.00
          endif
  4   continue
      recruf(1)=dorm*hatch
      resexf(2)=asex*brood*pnumf(mst)
      resexm(2)=(1.0-asex)*brood*pnumf(mst)*equ%fmale
      idtu1=int(deltau(1))
      rematf(2)=grorat(1)*rf(idtu1,1)*s(1)
      recruf(2)=resexf(2)+rematf(2)
      rematm(2)=grorat(1)*rm(idtu1,1)*s(1)
      recrum(2)=resexm(2)+rematm(2)
      do 9 j=1,mst
      idtuj=int(deltau(j))
      smu(j)=sm(idtuj,j)
 9    continue
      do 5 j=mst,mt
          if(delold(j-1).gt.deltau(j-1)) then
              idtuj=int(delold(j-1))
          else
              idtuj=int(deltau(j-1))
          endif
          rematf(j)=grorat(j-1)*rf(idtuj,j-1)*s(j-1)
          recruf(j)=rematf(j)
          rematm(j)=grorat(j-1)*rm(idtuj,j-1)*s(j-1)
          recrum(j)=rematm(j)
  5   continue

c     integration
      if(temp.le.equ%trest) then
          dorm=dorm+eggs
          eggs=0.00
      endif

      if(temp.ge.equ%tini) then
          taccu=taccu + dt*(temp-equ%tini)
      else
          taccu=0.0
      endif

      if(phper.ge.equ%phini) then
          paccu=paccu + dt*(phper-equ%phini)
      endif

      eggs=eggs + dt*(1.0-asex)*ferti*pnumf(mst)*emate
      eggs=eggs - dt*equ%pmegg*eggs
      if(dorm.gt.0.0) then
          dorm= dorm -hatch*dorm - dt*equ%pmegg*dorm
      endif
      if(dorm.lt.0.0) then
          dorm=0.0
      endif

c     ---- food dynamics, calculating balance

c     consumption
      totacon=0.0
      do 6 j=2,mst
          totacon=totacon+pnumf(j)*consu(j)
          totacon=totacon+pnumm(j)*consu(j)*equ%scalsx
  6   continue

c     supply
      if(lab.eq.1) then
          if(food1.ge.equ%fmax1) then
              fsupl1=0.0
          endif
          if(food2.ge.equ%fmax2) then
              fsupl2=0.0
          endif
      endif

c     balance, supply - consumption
      food1=food1+dt*(fsupl1-totacon*ratio)
      food2=food2+dt*(fsupl2-totacon*(1.0-ratio))
      if(food1.lt.0.0) then
          food1=0.00
      endif
      if(food2.lt.0.0) then
          food2=0.00
      endif
c     ------- end of food dynamics

c     ----- integration to update variables that accum info

c     accum starvation for effect on mortality thru ampmu
      do 7 j=2,mst
          starv(j)=starv(j)+dt*(tharv(j)-growth(j))
      if(starv(j).lt.0.0) then
          starv(j)=0.0
      endif
  7    continue

      do 8 j=1,mst
          s(j)=s(j)+dt*(s(j)*(grorat(j)*smu(j)-amu(j)))
          tau(j)=tau(j)+dt*(1.0-grorat(j))
          if (tau(j).lt.equ%taumin(j)) then
              tau(j)= equ%taumin(j)
          endif
          pnumf(j)=pnumf(j)+dt*(recruf(j)-rematf(j+1)-amu(j)*pnumf(j))
          if(pnumf(j).lt.0.0) then
              pnumf(j)=0.0
          endif
          pnumm(j)=pnumm(j)+dt*(recrum(j)-rematm(j+1)-amu(j)*pnumm(j))
          if(pnumm(j).lt.0.0) then
              pnumm(j)=0.0
          endif
          pnum(j)=pnumf(j)+pnumm(j)
  8   continue
      end

      subroutine update
      parameter(mt=4,mo=3000)
c     simulation control
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv
c     model variables
      common/recru1/ram(mo,mt),rm(mo,mt),raf(mo,mt),rf(mo,mt)
      common/recru2/recrum(mt),recruf(mt)
      common/growth/g(mo,mt),ga(mo,mt)
      common/compu1/tau(mt),eff(mt),deltau(mt),s(mt),delmax(mt)
      common/compu2/growth(mt),consu(mt),grorat(mt),delold(mt)
      common/compu5/tharv(mt),ampmu(mt),amu(mt),sm(mo,mt),sa(mo,mt)
      do 1 j=1,mst
          idmx=int(delmax(j))
          do 2 m=2,idmx
              ga(m,j)=g(m-1,j)
  2       continue
          g(1,j)=growth(j)
          do 3 m=2,idmx
              g(m,j)=ga(m,j)
  3       continue
          do 4 m=2,idmx
              raf(m,j)=rf(m-1,j)
              ram(m,j)=rm(m-1,j)
  4       continue
          rf(1,j)=recruf(j)
          rm(1,j)=recrum(j)
          do 5 m=2,idmx
              rf(m,j)=raf(m,j)
              rm(m,j)=ram(m,j)
  5       continue
          do 6 m=2,idmx
              sa(m,j)=sm(m-1,j)
  6       continue
          sm(1,j)=amu(j)
          do 7 m=2,idmx
              sm(m,j)=sa(m,j)
  7       continue
  1       continue
      end


      subroutine calc_index(simnstate)
      parameter(mt=4,mo=3000)
      common/control/dt,mdt,tfin,t,ict,time0,mst,lab,lfv
      common/state1/pnum(mt),pnumm(mt),pnumf(mt),totnum,eggs,dorm
      common/state2/food1,food2,foodt,taccu,paccu
      common/state/yd(10)
      common/index/x_max(10),x_min(10),x_avg(10),x_rms(10),x_end(10)

      integer simnstate

      if(t.eq.time0) then
          do 1 i=1,simnstate
              x_max(i)=yd(i)
              x_min(i)=yd(i)
              x_avg(i)=yd(i)
              x_rms(i)=0.00
              x_end(i)=yd(i)
  1       continue
      else
          do 2 i=1,simnstate
              if(yd(i).ge.x_max(i)) x_max(i)=yd(i)
              if(yd(i).le.x_min(i)) x_min(i)=yd(i)
              x_avg(i)=(x_avg(i)*(t-time0)/dt+yd(i))*(dt/(dt+t-time0))
              rms = (yd(i) - x_avg(i))**2
              x_rms(i)=(x_rms(i)*(t-time0)/dt+rms)*(dt/(dt+t-time0))
              x_end(i)=yd(i)
  2      continue
      endif
      end
