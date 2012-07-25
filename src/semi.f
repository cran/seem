      subroutine semi
C     '******************************
C     ' Miguel F. Acevedo
C     ' Program to simulate semi-Markov process 
C     ' Uses linear chain trick to emulate distributed delays
C     ' and fixed delays concatenated with them
C     ' System scaled in time to decades (ten years). Iterations dt in years
C     '*********************************
C     '
c     ----------parameter declarations--------------------      
      PARAMETER(MT=20,MTF=500,MO=10,ML=300,MTL=100)
C     mt is max number of states, mtf is max time of simulation,
c     mo is max order of distributed lags, 
C     ml is max order of fixed lags(dim of largest fixed delay array)
C     mtl is max number of fixed lags, (max number of transitions 
c         which have a fixed latency)
      COMMON/CONTRO/MST,DT,MDT,TFIN,T,ICT
c     simulation control parameters,      
      COMMON/STATES/XT(MT,MT,MO),XF(MT),X(MT),XL(MTL,ML),EXIT(MT,MT)
c     state parameters
      COMMON/PARAME/P(MT,MT), D(MT,MT), K(MT,MT), AV(MT,MT),
     &    IQ(MT,MT), QD(MT,MT),MARK(MT,MT)
c     model parameters: prob, rates, order dist lags, etc
c     more info provided in readin routine
c     used for calculations
c     graphics control parameters
c     variables for managing exiting fractions from each state
c     ----------------main--------------------
      CHARACTER*20 filein
c     open the input file
      OPEN(4,FILE="semi_inp.txt")
C     read from filein: the number states, time step, time to print, final time 
      CALL READIND
c     read parameters and initial condition
      CALL READPARAM
C     the following subroutine selects only the non zero fixed lags  
c     to save time during execution
      CALL SELECT
C     sets the initial values for time and states
      CALL INITIALIZE
c     simulation run
      CALL SIMULA
c     finish
      CLOSE(4)
      END
c     *******************end of main*******************
c     -------------------------------------------------------------      
c     The subroutine readIND reads in number of states
c     time step, time to print/plot and simulation time
      SUBROUTINE READIND
      COMMON/CONTRO/MST,DT,MDT,TFIN,T,ICT
      read(4,2050)MST
 2050 FORMAT(I7)     
      read(4,2060)DT,MDT,TFIN
 2060 FORMAT(F7.2,I7,F7.2)
      END
C     ----------------------------------------------
C     The subroutine readPARAM reads in the initial condition and the
c     parameter values
      SUBROUTINE READPARAM
      PARAMETER(MT=20,MTF=500,MO=10,ML=300,MTL=100)
      COMMON/CONTRO/MST,DT,MDT,TFIN,T,ICT
      COMMON/PARAME/P(MT,MT), D(MT,MT), K(MT,MT), AV(MT,MT),
     &     IQ(MT,MT), QD(MT,MT),MARK(MT,MT)
      COMMON/STATES/ XT(MT,MT,MO),XF(MT),X(MT),XL(MTL,ML),EXIT(MT,MT)
      read(4,1000) (XF(I),I=1,MST)
c     the xf vector contains the initial condition for the simulation
 1000 FORMAT(16F7.3)
 1001 FORMAT(16F7.2)
      DO 10 J=1,MST
      read(4,1000) (P(J,I),I=1,MST)
  10  CONTINUE 
C     'transition probabilities
      DO 20 J=1,MST
      read(4,1002) (K(J,I),I=1,MST)
 1002 FORMAT(16I7)      
  20  CONTINUE
C     'order of distributed time lags i.e. number of first order stages         C     one value for every pair of states
      DO 30 J=1,MST      
      read(4,1001) (AV(J,I),I=1,MST)
  30  CONTINUE 
C     'mean values for distributed lags                 
c     with the order and the mean, the rates are computed for every pair of 
C     states, in the following nested loop by the ratio order/mean 
      DO 50 J=1,MST                                 
      DO 55 I=1,MST                               
       IF (AV(J,I).NE.0) THEN
         D(J,I)=K(J,I)/AV(J,I)
        ELSE
         D(J,I)=0.0
       ENDIF
   55     CONTINUE                                     
   50   CONTINUE                                       
C     'rates
c     the values for the fixed latency for every pair of states;   
C     it is converted to an integer so that an array can be managed later
      DO 70 J=1,MST                     
      read(4,1001) (QD(J,I),I=1,MST)
      DO 75 I=1,MST                   
      FDT=1/DT  
      IDT=INT(FDT)
      DDT=FDT-IDT
      IF (DDT.GT.0.5) IDT=IDT+1
      IQ(J,I)=IDT*QD(J,I)
  75    CONTINUE                      
  70  CONTINUE                        
C     'fixed lags
      END
C       -------------------------------------------------------
c       the following routine scans all the transition pairs, checks whether
C       there is fixed delay for that transition, to increase a counter
C       kcuenta of the number of non-zero fixed delays and to assign a 
C       number 'mark(i,j)' to every pair of states i,j 
C
      SUBROUTINE SELECT
      PARAMETER(MT=20,MTF=500,MO=10,ML=300,MTL=100)
      COMMON/CONTRO/MST,DT,MDT,TFIN,T,ICT
      COMMON/PARAME/P(MT,MT), D(MT,MT), K(MT,MT), AV(MT,MT),
     &     IQ(MT,MT), QD(MT,MT),MARK(MT,MT)
  77      KCUENTA=1
      DO 73 J=1,MST
       DO 73 I=1,MST
        IF(IQ(J,I).NE.0) THEN
         MARK(I,J)=KCUENTA
          KCUENTA=KCUENTA+1
        ENDIF
 73      CONTINUE
      END
C     ---------------------------------------------------
C     the following routine initializes the simulation      
c     to t=0, x=init cond, etc 
c      
      SUBROUTINE INITIALIZE
      PARAMETER(MT=20,MTF=500,MO=10,ML=300,MTL=100)   
      COMMON/CONTRO/MST,DT,MDT,TFIN,T,ICT
      COMMON/PARAME/P(MT,MT), D(MT,MT), K(MT,MT), AV(MT,MT),
     &     IQ(MT,MT), QD(MT,MT),MARK(MT,MT)
      COMMON/STATES/ XT(MT,MT,MO),XF(MT),X(MT),XL(MTL,ML),EXIT(MT,MT)
c      initial time for simulation is t=0
      T=0.0
      ICT=0     
C         ict is another counter of iterations before printing or plot      
      DO 89 I=1,MST
      X(I)=XF(I)
 89   CONTINUE
C     this last loop initilaize the state variables to the initial condit
      DO 700 J=1,MST
      DO 700 I=1,MST
c       all of the exiting fractions will be made zero at time0 
      EXIT(J,I)=0.0 
      DO 710 L=1,MO
        XT(J,I,L)=0.0
c         this statement makes equal to zero all the stages of distributed 
C         delay up to a maximun value mo
 710     CONTINUE
       DO 700 L=1,ML
         XL(MARK(J,I),L)=0.0
C          the fixed lag stages are also made zero up to ml; 
c          only the ones which have a nonzero latency
 700  CONTINUE
      END
c     ---------------------------------------------
c     This is the simulation subroutine      
c     it contains the iteration loop      
      SUBROUTINE SIMULA           
      PARAMETER(MT=20,MTF=500,MO=10,ML=300,MTL=100)
      COMMON/CONTRO/MST,DT,MDT,TFIN,T,ICT
      COMMON/STATES/ XT(MT,MT,MO),XF(MT),X(MT),XL(MTL,ML),EXIT(MT,MT)
      COMMON/PARAME/P(MT,MT), D(MT,MT), K(MT,MT), AV(MT,MT),
     &       IQ(MT,MT), QD(MT,MT),MARK(MT,MT)

c     ******************starts simulation loop************
c      print routine is a function in C to avoid fortran write statements
 499    call semiprint(MST,DT,TFIN,XF,P,T,X) 
C               ' calculation loop
c               *******inner iteration loop********* 
 500            CALL calcul
c               take another time step dt
                T=T+DT 
                ICT=ICT+1
c               check if it is time to save variables 
c               for printing or plotting       
                IF (ICT.LT.MDT) GOTO 500
c               *********ends inner iteration loop*****
       ICT=0
c      if not yet the repeat iteration loop, but if yes then reset counter
C     ---------iterates  if t<=tfin----------------
C     'aggregation loop
       CALL AGREGA
      ABT=ABS(T-TFIN)
       IF (T.LT.TFIN.OR.ABT.LT.0.001) GOTO 499
c     *********************end of simulation loop***************
      END
C     -------------------------------------------------------
      SUBROUTINE CALCUL
C     *************************************
C     this is the iteration loop
C     *************************************
      PARAMETER(MT=20,MTF=500,MO=10,ML=300,MTL=100)
      COMMON/CONTRO/MST,DT,MDT,TFIN,T,ICT
      COMMON/STATES/ XT(MT,MT,MO),XF(MT),X(MT),XL(MTL,ML),EXIT(MT,MT)
      COMMON/PARAME/P(MT,MT), D(MT,MT), K(MT,MT), AV(MT,MT),
     &              IQ(MT,MT), QD(MT,MT),MARK(MT,MT)

      dimension DEC(MT),DINC(MT),DINCEXIT(MT,MT),DECEXIT(MT,MT)
c     the following loop calculates the exiting fractions out of
c     source state i and going to destination state j according to 
c     imbedded probabilities
      DO 98 I=1,MST
        DO 98 J=1,MST
        DINCEXIT(J,I)=P(J,I)*XF(I)
        EXIT(J,I)=EXIT(J,I)+DINCEXIT(J,I)
 98    CONTINUE 
c 
 500  DO 145 I=1,MST
C        loop for all source states
       DEC(I)=0.0
       DINC(I)=0.0
C        these are the decrement and increment values for each state
C        which are going to be calculated at every iteration         
       DO 120 J=1,MST
c        loop for all destination states  
C            first the decrement of the source state will be equal to
c            the exit term calculated previously (fraction leaving state
c            i and going to state j according to probability j,i) 
           DEC(I)=DEC(I)+DINCEXIT(J,I)
           IF (IQ(J,I).EQ.0) THEN
           DECEXIT(J,I)=DT*D(J,I)*EXIT(J,I)  
c            this is decrement of exit in case there is no fixed latency,
C            note that it is fraction exiting multiplied by the first order rate
C            also dt comes in because this is euler integration method
C            of the chain of first order lags which make up the distributed
C            lag. 
           ELSE
             DECEXIT(J,I)=EXIT(J,I)
c            on the other hand if there is fixed latency the  
C            fraction exiting the state is the decrement for the state 
           ENDIF
c           this endif corresponds to end of establishing the decrement
C            the following term is the increment in state i due to incoming                 
C            recruits from last stage of distributed lag coming from state j 
C            The increment is set to zero if no distributed lag 
c            Note: a nonzero transition proba assumes a distributed
c            lag of at least order 1 
         NK=K(I,J)
         IF(NK.GT.0) THEN
            IF (NK.EQ.1) THEN
                IF (IQ(I,J).EQ.0) THEN
                    DINCD=DT*D(I,J)*EXIT(I,J)
                 ELSE
                    DINCD=DT*D(I,J)*XL(MARK(J,I),IQ(I,J))
                ENDIF
             ELSE
              DINCD=DT*D(I,J)*XT(J,I,NK-1)
            ENDIF
         ELSE
            DINCD=0.0
         ENDIF
         DINC(I)=DINC(I)+DINCD
 120     CONTINUE
c        once the decrement and increment are calculated the values for the state
c        variables are going to be updated: add the inc and substract the dec
       xF(I)=XF(I)-DEC(I)+DINC(I)
  145       continue
       do 110 i=1,mst
c       now the qties in each stage of the fixed delay are going to
c        be promoted to the next stage. First check that there is a
c        non zero fixed delay.
       DO 110 J=1,MST
       FRACT=DT*D(J,I)
           nk1=k(j,i)
c       nk1 is just an aux var for the number of stages of the
c        distributed delay; 
           IF (NK1.GT.2) THEN
c       Now, if the distributed delay is of order >2 then compute other
c        stages, if not then done and go to update
C       All of the stages of the distributed delay
c        will now be computed from the first order process
           DO 118 N=1,NK1-2
                  nq=nk1-n
                  XTIJ1 = XT(I,J,Nq-1)
                  XTIJ1F=XTIJ1*FRACT
                  XTIJN = XT(I,J,Nq)
                  XTIJNF= XTIJN*(1-FRACT)
                  xT(I,J,Nq) = XTIJ1F+XTIJNF

 118        CONTINUE
           ENDIF
           IF (IQ(J,I).NE.0) THEN
C         the first stage of gamma delay is computed using the
C         fraction leaving the last stage of the fixed delay
             NK1=K(J,I)
             IF (NK1.GT.1) THEN
                  XLIQ= FRACT*XL(MARK(I,J),IQ(J,I))
                  XTIJ= XT(I,J,1)*(1-FRACT)
                  xT(I,J,1)=XTIJ+XLIQ
             ENDIF         
c            now the ending or last  position of the fixed delay
c             is modified by promotion from the next to last position
c             and by the rate of entrance to the distributed delay
            IF(IQ(J,I).GE.2) THEN  
                 xL(MARK(I,J),iq(j,i))=XL(MARK(I,J),IQ(J,I)-1)+
     &           XL(MARK(I,J),IQ(J,I))*(1-FRACT)
                ELSE
c             the following corresponds to the special case when
c             the fixed delay is of length 1, in which case there is
c             no promotion from a previous stage, but from the 
c             exiting fraction out of the state i
                 xL(MARK(I,J),iq(j,i))=exit(j,i)+
     &           XL(MARK(I,J),IQ(J,I))*(1-FRACT)
            ENDIF
            IF (IQ(J,I).GE.3) THEN
c        this is done so that the mid entries are obtained by shifting 
c        in the value in the previous stage
                DO 140 L=1,IQ(J,I)-2
                iqml=iq(j,i)-l
                xL(MARK(I,J),iqml) = XL(MARK(I,J),iqml-1)
  140           CONTINUE
            ENDIF
C           and the fraction exiting the satte moves on to first entry
c           of fixed lag except when iq(j,i)=1 (because it was already)
c           calculated)
           if(iq(j,i).ge.2) then
             xL(MARK(I,J),1)=EXIT(J,I)
           endif
C            -------------------------------------
          ELSE
c           or from the fraction exiting the state   
c           when there is no fixed latency
                NK1=K(J,I)
                IF(NK1.GT.1) THEN
                XTIJ= XT(I,J,1)*(1-FRACT) 
                xT(I,J,1) = FRACT*EXIT(J,I) +XTIJ
                ENDIF
          ENDIF
 110      CONTINUE
c       the following loop will update the exiting fractions by        
c       subtracting the decrement
      DO 111 I=1,MST
      DO 111 J=1,MST
      EXIT(J,I)=EXIT(J,I)-DECEXIT(J,I)
 111  CONTINUE
      
      END
C       ----------------------------------------------------          
       SUBROUTINE AGREGA
c     Now the values in the fixed delays and distributed delays for
c     the every transition from state i to state j
c     have to be lumped in one state variable i
      PARAMETER(MT=20,MTF=500,MO=10,ML=300,MTL=100)
      COMMON/CONTRO/MST,DT,MDT,TFIN,T,ICT
      COMMON/STATES/ XT(MT,MT,MO),XF(MT),X(MT),XL(MTL,ML),EXIT(MT,MT)
      COMMON/PARAME/P(MT,MT), D(MT,MT), K(MT,MT), AV(MT,MT),
     &                     IQ(MT,MT), QD(MT,MT),MARK(MT,MT)
      DO 200 I=1,MST
       X(I)=XF(I)
c       add the source variable
      DO 200 J=1,MST
c       now add the exiting fractions
         X(I)=x(i)+EXIT(J,I)
         IF(IQ(J,I).EQ.0) GOTO 221
         DO 220 L=1,IQ(J,I)
           X(I)=X(I) +XL(MARK(I,J),L)
c        add the fixed lag values 
 220       CONTINUE
 221         DO 230 N=1,K(J,I)-1
c        add the distributed delay values
        X(I)=X(I) +XT(I,J,N)
 230       CONTINUE
200   CONTINUE
      END
