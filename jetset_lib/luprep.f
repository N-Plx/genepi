 
C*********************************************************************

C...LUPREP
C...Rearranges partons along strings. 
C...Allows small systems to collapse into one or two particles. 
C...Checks flavours and colour singlet invarient masses.
 
      SUBROUTINE LUPREP(IP)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(D)
      INTEGER LUK,LUCHGE,LUCOMP
C...Commonblocks.
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      SAVE /LUJETS/,/LUDAT1/,/LUDAT2/,/LUDAT3/
C...Local arrays.
      DIMENSION DPS(5),DPC(5),UE(3),PG(5),
     &E1(3),E2(3),E3(3),E4(3),ECL(3)
 
C...Function to give four-product.
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)
 
C...Rearrange parton shower product listing along strings: begin loop.
      I1=N
      DO 130 MQGST=1,2
        DO 120 I=MAX(1,IP),N
          IF(K(I,1).NE.3) GOTO 120
          KC=LUCOMP(K(I,2))
          IF(KC.EQ.0) GOTO 120
          KQ=KCHG(KC,2)
          IF(KQ.EQ.0.OR.(MQGST.EQ.1.AND.KQ.EQ.2)) GOTO 120
 
C...Pick up loose string end.
          KCS=4
          IF(KQ*ISIGN(1,K(I,2)).LT.0) KCS=5
          IA=I
          NSTP=0
  100     NSTP=NSTP+1
          IF(NSTP.GT.4*N) THEN
            CALL LUERRM(14,'(LUPREP:) caught in infinite loop')
            RETURN
          ENDIF
 
C...Copy undecayed parton.
          IF(K(IA,1).EQ.3) THEN
            IF(I1.GE.MSTU(4)-MSTU(32)-5) THEN
              CALL LUERRM(11,'(LUPREP:) no more memory left in LUJETS')
              RETURN
            ENDIF
            I1=I1+1
            K(I1,1)=2
            IF(NSTP.GE.2.AND.KCHG(LUCOMP(K(IA,2)),2).NE.2) K(I1,1)=1
            K(I1,2)=K(IA,2)
            K(I1,3)=IA
            K(I1,4)=0
            K(I1,5)=0
            DO 110 J=1,5
              P(I1,J)=P(IA,J)
              V(I1,J)=V(IA,J)
  110       CONTINUE
            K(IA,1)=K(IA,1)+10
            IF(K(I1,1).EQ.1) GOTO 120
          ENDIF
 
C...Go to next parton in colour space.
          IB=IA
          IF(MOD(K(IB,KCS)/MSTU(5)**2,2).EQ.0.AND.MOD(K(IB,KCS),MSTU(5))
     &    .NE.0) THEN
            IA=MOD(K(IB,KCS),MSTU(5))
            K(IB,KCS)=K(IB,KCS)+MSTU(5)**2
            MREV=0
          ELSE
            IF(K(IB,KCS).GE.2*MSTU(5)**2.OR.MOD(K(IB,KCS)/MSTU(5),
     &      MSTU(5)).EQ.0) KCS=9-KCS
            IA=MOD(K(IB,KCS)/MSTU(5),MSTU(5))
            K(IB,KCS)=K(IB,KCS)+2*MSTU(5)**2
            MREV=1
          ENDIF
          IF(IA.LE.0.OR.IA.GT.N) THEN
            CALL LUERRM(12,'(LUPREP:) colour rearrangement failed')
            RETURN
          ENDIF
          IF(MOD(K(IA,4)/MSTU(5),MSTU(5)).EQ.IB.OR.MOD(K(IA,5)/MSTU(5),
     &    MSTU(5)).EQ.IB) THEN
            IF(MREV.EQ.1) KCS=9-KCS
            IF(MOD(K(IA,KCS)/MSTU(5),MSTU(5)).NE.IB) KCS=9-KCS
            K(IA,KCS)=K(IA,KCS)+2*MSTU(5)**2
          ELSE
            IF(MREV.EQ.0) KCS=9-KCS
            IF(MOD(K(IA,KCS),MSTU(5)).NE.IB) KCS=9-KCS
            K(IA,KCS)=K(IA,KCS)+MSTU(5)**2
          ENDIF
          IF(IA.NE.I) GOTO 100
          K(I1,1)=1
  120   CONTINUE
  130 CONTINUE
      N=I1
 
C...Done if no checks on small-mass systems.
      IF(MSTJ(14).LT.0) RETURN
      IF(MSTJ(14).EQ.0) GOTO 540
 
C...Find lowest-mass colour singlet jet system.
      NS=N
  140 NSIN=N-NS
      PDMIN=1.+PARJ(32)
      IC=0
      DO 190 I=MAX(1,IP),N
        IF(K(I,1).NE.1.AND.K(I,1).NE.2) THEN
        ELSEIF(K(I,1).EQ.2.AND.IC.EQ.0) THEN
          NSIN=NSIN+1
          IC=I
          DO 150 J=1,4
            DPS(J)=P(I,J)
  150     CONTINUE
          MSTJ(93)=1
          DPS(5)=ULMASS(K(I,2))
        ELSEIF(K(I,1).EQ.2) THEN
          DO 160 J=1,4
            DPS(J)=DPS(J)+P(I,J)
  160     CONTINUE
        ELSEIF(IC.NE.0.AND.KCHG(LUCOMP(K(I,2)),2).NE.0) THEN
          DO 170 J=1,4
            DPS(J)=DPS(J)+P(I,J)
  170     CONTINUE
          MSTJ(93)=1
          DPS(5)=DPS(5)+ULMASS(K(I,2))
          PD=SQRT(MAX(0D0,DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2))-
     &    DPS(5)
          IF(PD.LT.PDMIN) THEN
            PDMIN=PD
            DO 180 J=1,5
              DPC(J)=DPS(J)
  180       CONTINUE
            IC1=IC
            IC2=I
          ENDIF
          IC=0
        ELSE
          NSIN=NSIN+1
        ENDIF
  190 CONTINUE
 
C...Done if lowest-mass system above threshold for string frag.
      IF(PDMIN.GE.PARJ(32)) GOTO 540
 
C...Fill small-mass system as cluster.
      NSAV=N
      PECM=SQRT(MAX(0D0,DPC(4)**2-DPC(1)**2-DPC(2)**2-DPC(3)**2))
      K(N+1,1)=11
      K(N+1,2)=91
      K(N+1,3)=IC1
      P(N+1,1)=DPC(1)
      P(N+1,2)=DPC(2)
      P(N+1,3)=DPC(3)
      P(N+1,4)=DPC(4)
      P(N+1,5)=PECM
 
C...Set up history, assuming cluster -> 2 hadrons.
      NBODY=2
      K(N+1,4)=N+2
      K(N+1,5)=N+3
      K(N+2,1)=1
      K(N+3,1)=1
      IF(MSTU(16).NE.2) THEN
        K(N+2,3)=N+1
        K(N+3,3)=N+1
      ELSE
        K(N+2,3)=IC1
        K(N+3,3)=IC2
      ENDIF
      K(N+2,4)=0
      K(N+3,4)=0
      K(N+2,5)=0
      K(N+3,5)=0
      V(N+1,5)=0.
      V(N+2,5)=0.
      V(N+3,5)=0.
 
C...Form two particles from flavours of lowest-mass system, if feasible.
      NTRY = 0
  200 NTRY = NTRY + 1
C...Open string.
      IF(IABS(K(IC1,2)).NE.21) THEN
        KC1=LUCOMP(K(IC1,2))
        KC2=LUCOMP(K(IC2,2))
        IF(KC1.EQ.0.OR.KC2.EQ.0) GOTO 540
        KQ1=KCHG(KC1,2)*ISIGN(1,K(IC1,2))
        KQ2=KCHG(KC2,2)*ISIGN(1,K(IC2,2))
        IF(KQ1+KQ2.NE.0) GOTO 540
C...Start with qq, if there is one. Only allow for rank 1 popcorn meson
  210   K1=K(IC1,2)
        IF(IABS(K(IC2,2)).GT.10) K1=K(IC2,2)
        MSTU(125)=0
        CALL LUDCYK(K1,0,KFLN,K(N+2,2))
        CALL LUDCYK(K(IC1,2)+K(IC2,2)-K1,-KFLN,KFLDMP,K(N+3,2))
        IF(K(N+2,2).EQ.0.OR.K(N+3,2).EQ.0) GOTO 210
C...Closed string.
      ELSE
        IF(IABS(K(IC2,2)).NE.21) GOTO 540
C...No room for popcorn mesons in closed string -> 2 hadrons.
        MSTU(125)=0
  220   CALL LUDCYK(1+INT((2.+PARJ(2))*RLU(0)),0,KFLN,KFDMP)
        CALL LUDCYK(KFLN,0,KFLM,K(N+2,2))
        CALL LUDCYK(-KFLN,-KFLM,KFLDMP,K(N+3,2))
        IF(K(N+2,2).EQ.0.OR.K(N+3,2).EQ.0) GOTO 220
      ENDIF
      P(N+2,5)=ULMASS(K(N+2,2))
      P(N+3,5)=ULMASS(K(N+3,2))
 
C...If it does not work: try again (a number of times), give up
C...(if no place to shuffle momentum), or form one hadron.
      IF(P(N+2,5)+P(N+3,5)+PARJ(64).GE.PECM) THEN
        IF(NTRY.LT.MSTJ(17)) THEN
          GOTO 200
        ELSEIF(NSIN.EQ.1) THEN
          GOTO 540
        ELSE
          GOTO 290
        END IF
      END IF
 
C...Perform two-particle decay of jet system.
C...First step: find reference axis in decaying system rest frame.
C...(Borrow slot N+2 for temporary direction.)
      DO 230 J=1,4
        P(N+2,J)=P(IC1,J)
  230 CONTINUE
      DO 250 I=IC1+1,IC2-1
        IF((K(I,1).EQ.1.OR.K(I,1).EQ.2).AND.
     &  KCHG(LUCOMP(K(I,2)),2).NE.0) THEN
          FRAC1=FOUR(IC2,I)/(FOUR(IC1,I)+FOUR(IC2,I))
          DO 240 J=1,4
            P(N+2,J)=P(N+2,J)+FRAC1*P(I,J)
  240     CONTINUE
        ENDIF
  250 CONTINUE
      CALL PYROBO(N+2,N+2,0.,0.,-SNGL(DPC(1)/DPC(4)),
     &-SNGL(DPC(2)/DPC(4)),
     &-SNGL(DPC(3)/DPC(4)))
      THE1=ULANGL(P(N+2,3),SQRT(P(N+2,1)**2+P(N+2,2)**2))
      PHI1=ULANGL(P(N+2,1),P(N+2,2))
 
C...Second step: generate isotropic/anisotropic decay.
      PA=SQRT((PECM**2-(P(N+2,5)+P(N+3,5))**2)*(PECM**2-
     &(P(N+2,5)-P(N+3,5))**2))/(2.*PECM)
  260 UE(3)=RLU(0)
      PT2=(1.-UE(3)**2)*PA**2
      IF(MSTJ(16).LE.0) THEN
        PREV=0.5
      ELSE
        IF(EXP(-PT2/(2.*PARJ(21)**2)).LT.RLU(0)) GOTO 260
        PR1=P(N+2,5)**2+PT2
        PR2=P(N+3,5)**2+PT2
        ALAMBD=SQRT(MAX(0.,(PECM**2-PR1-PR2)**2-4.*PR1*PR2))
        PREVCF=PARJ(42)
        IF(MSTJ(11).EQ.2) PREVCF=PARJ(39)
        PREV=1./(1.+EXP(MIN(50.,PREVCF*ALAMBD)))
      ENDIF
      IF(RLU(0).LT.PREV) UE(3)=-UE(3)
      PHI=PARU(2)*RLU(0)
      UE(1)=SQRT(1.-UE(3)**2)*COS(PHI)
      UE(2)=SQRT(1.-UE(3)**2)*SIN(PHI)
      DO 270 J=1,3
        P(N+2,J)=PA*UE(J)
        P(N+3,J)=-PA*UE(J)
  270 CONTINUE
      P(N+2,4)=SQRT(PA**2+P(N+2,5)**2)
      P(N+3,4)=SQRT(PA**2+P(N+3,5)**2)
 
C...Third step: move back to event frame and set production vertex.
      CALL PYROBO(N+2,N+3,THE1,PHI1,SNGL(DPC(1)/DPC(4)),
     &SNGL(DPC(2)/DPC(4)),
     &SNGL(DPC(3)/DPC(4)))
      DO 280 J=1,4
        V(N+1,J)=V(IC1,J)
        V(N+2,J)=V(IC1,J)
        V(N+3,J)=V(IC2,J)
  280 CONTINUE
      N=N+3
      GOTO 520
 
C...Else form one particle, if possible.
  290 NBODY=1
      K(N+1,5)=N+2
      DO 300 J=1,4
        V(N+1,J)=V(IC1,J)
        V(N+2,J)=V(IC1,J)
  300 CONTINUE
 
C...Select hadron flavour from available quark flavours.
  310 IF(IABS(K(IC1,2)).GT.100.AND.IABS(K(IC2,2)).GT.100) THEN
        GOTO 540
      ELSEIF(IABS(K(IC1,2)).NE.21) THEN
        CALL LUKFDI(K(IC1,2),K(IC2,2),KFLDMP,K(N+2,2))
      ELSE
        KFLN=1+INT((2.+PARJ(2))*RLU(0))
        CALL LUKFDI(KFLN,-KFLN,KFLDMP,K(N+2,2))
      ENDIF
      IF(K(N+2,2).EQ.0) GOTO 310
      P(N+2,5)=ULMASS(K(N+2,2))
 
C...Use old algorithm for E/p conservation? (EN)
      IF (MSTJ(16).LE.0) GOTO 480
 
C...Find the string piece closest to the cluster by a loop
C...over the undecayed partons not in present cluster. (EN)
      DGLOMI=1D30
      IBEG=0
      I0=0
      DO 340 I1=MAX(1,IP),N-1
        IF(I1.GE.IC1-1.AND.I1.LE.IC2) THEN
          I0=0
        ELSEIF(K(I1,1).EQ.2) THEN
          IF(I0.EQ.0) I0=I1
          I2=I1
  320     I2=I2+1
          IF(K(I2,1).GT.10) GOTO 320  
          IF(KCHG(LUCOMP(K(I2,2)),2).EQ.0) GOTO 320
 
C...Define velocity vectors e1, e2, ecl and differences e3, e4.
          DO 330 J=1,3
            E1(J)=P(I1,J)/P(I1,4)
            E2(J)=P(I2,J)/P(I2,4)
            ECL(J)=P(N+1,J)/P(N+1,4)
            E3(J)=E2(J)-E1(J)
            E4(J)=ECL(J)-E1(J)
  330     CONTINUE
 
C...Calculate minimal D=(e4-alpha*e3)**2 for 0<alpha<1.
          E3S=E3(1)**2+E3(2)**2+E3(3)**2
          E4S=E4(1)**2+E4(2)**2+E4(3)**2
          E34=E3(1)*E4(1)+E3(2)*E4(2)+E3(3)*E4(3)
          IF(E34.LE.0.) THEN
            DDMIN=E4S
          ELSEIF(E34.LT.E3S) THEN
            DDMIN=E4S-E34**2/E3S
          ELSE
            DDMIN=E4S-2.*E34+E3S
          ENDIF
 
C...Is this the smallest so far?
          IF(DDMIN.LT.DGLOMI) THEN
            DGLOMI=DDMIN
            IBEG=I0
            IPCS=I1
          ENDIF
        ELSEIF(K(I1,1).EQ.1.AND.KCHG(LUCOMP(K(I1,2)),2).NE.0) THEN
          I0=0
        ENDIF
  340 CONTINUE
 
C... Check if there are any strings to connect to the new gluon. (EN)
      IF (IBEG.EQ.0) GOTO 480
 
C...Delta_m = m_clus - m_had > 0: emit a 'gluon' (EN)
      IF (P(N+1,5).GE.P(N+2,5)) THEN
 
C...Construct 'gluon' that is needed to put hadron on the mass shell.
        FRAC=P(N+2,5)/P(N+1,5)
        DO 350 J=1,5
          P(N+2,J)=FRAC*P(N+1,J)
          PG(J)=(1.-FRAC)*P(N+1,J)
  350   CONTINUE
 
C... Copy string with new gluon put in.
        N=N+2
        I=IBEG-1
  360   I=I+1
        IF(K(I,1).NE.1.AND.K(I,1).NE.2) GOTO 360
        IF(KCHG(LUCOMP(K(I,2)),2).EQ.0) GOTO 360
        N=N+1
        DO 370 J=1,5
          K(N,J)=K(I,J)
          P(N,J)=P(I,J)
          V(N,J)=V(I,J)
  370   CONTINUE
        K(I,1)=K(I,1)+10
        K(I,4)=N
        K(I,5)=N
        K(N,3)=I
        IF(I.EQ.IPCS) THEN
          N=N+1
          DO 380 J=1,5
            K(N,J)=K(N-1,J)
            P(N,J)=PG(J)
            V(N,J)=V(N-1,J)
  380     CONTINUE
          K(N,2)=21
          K(N,3)=NSAV+1
        ENDIF
        IF(K(I,1).EQ.12) GOTO 360
        GOTO 520
 
C...Delta_m = m_clus - m_had < 0: have to absorb a 'gluon' instead,
C...from string piece endpoints.
      ELSE
 
C...Begin by copying string that should give energy to cluster.
        N=N+2
        I=IBEG-1
  390   I=I+1
        IF(K(I,1).NE.1.AND.K(I,1).NE.2) GOTO 390
        IF(KCHG(LUCOMP(K(I,2)),2).EQ.0) GOTO 390
        N=N+1
        DO 400 J=1,5
          K(N,J)=K(I,J)
          P(N,J)=P(I,J)
          V(N,J)=V(I,J)
  400   CONTINUE
        K(I,1)=K(I,1)+10
        K(I,4)=N
        K(I,5)=N
        K(N,3)=I
        IF(I.EQ.IPCS) I1=N
        IF(K(I,1).EQ.12) GOTO 390
        I2=I1+1
 
C...Set initial Phad.
        DO 410 J=1,4
          P(NSAV+2,J)=P(NSAV+1,J)
  410   CONTINUE
 
C...Calculate Pg, a part of which will be added to Phad later. (EN)
  420   IF(MSTJ(16).EQ.1) THEN
          ALPHA=1.
          BETA=1.
        ELSE
          ALPHA=FOUR(NSAV+1,I2)/FOUR(I1,I2)
          BETA=FOUR(NSAV+1,I1)/FOUR(I1,I2)
        ENDIF
        DO 430 J=1,4
          PG(J)=ALPHA*P(I1,J)+BETA*P(I2,J)
  430   CONTINUE
        PG(5)=SQRT(MAX(1E-20,PG(4)**2-PG(1)**2-PG(2)**2-PG(3)**2))
 
C..Solve 2nd order equation, use the best (smallest) solution. (EN)
        PMSCOL=P(NSAV+2,4)**2-P(NSAV+2,1)**2-P(NSAV+2,2)**2-
     &  P(NSAV+2,3)**2
        PCLPG=(P(NSAV+2,4)*PG(4)-P(NSAV+2,1)*PG(1)-
     &  P(NSAV+2,2)*PG(2)-P(NSAV+2,3)*PG(3))/PG(5)**2
        DELTA=SQRT(PCLPG**2+(P(NSAV+2,5)**2-PMSCOL)/PG(5)**2)-PCLPG
 
C...If all gluon energy eaten, zero it and take a step back.
        ITER=0
        IF(DELTA*ALPHA.GT.1..AND.I1.GT.NSAV+3) THEN
          ITER=1
          DO 440 J=1,4
            P(NSAV+2,J)=P(NSAV+2,J)+P(I1,J)
            P(I1,J)=0.
  440     CONTINUE
          P(I1,5)=0.
          K(I1,1)=K(I1,1)+10
          I1=I1-1
        ENDIF
        IF(DELTA*BETA.GT.1..AND.I2.LT.N) THEN
          ITER=1
          DO 450 J=1,4
            P(NSAV+2,J)=P(NSAV+2,J)+P(I2,J)
            P(I2,J)=0.
  450     CONTINUE
          P(I2,5)=0.
          K(I2,1)=K(I2,1)+10
          I2=I2+1
        ENDIF
        IF(ITER.EQ.1) GOTO 420
 
C...If also all endpoint energy eaten, revert to old procedure.
        IF((1.-DELTA*ALPHA)*P(I1,4).LT.P(I1,5).OR.
     &  (1.-DELTA*BETA)*P(I2,4).LT.P(I2,5)) THEN
          DO 460 I=NSAV+3,N
            IM=K(I,3)
            K(IM,1)=K(IM,1)-10
            K(IM,4)=0
            K(IM,5)=0
  460     CONTINUE
          N=NSAV
          GOTO 480
        ENDIF
 
C... Construct the collapsed hadron and modified string partons.
        DO 470 J=1,4
          P(NSAV+2,J)=P(NSAV+2,J)+DELTA*PG(J)
          P(I1,J)=(1.-DELTA*ALPHA)*P(I1,J)
          P(I2,J)=(1.-DELTA*BETA)*P(I2,J)
  470   CONTINUE
          P(I1,5)=(1.-DELTA*ALPHA)*P(I1,5)
          P(I2,5)=(1.-DELTA*BETA)*P(I2,5)
 
C...Finished with string collapse in new scheme.
        GOTO 520
      ENDIF
 
C... Use old algorithm; by choice or when in trouble.
  480 CONTINUE
C...Find parton/particle which combines to largest extra mass.
      IR=0
      HA=0.
      HSM=0.
      DO 500 MCOMB=1,3
        IF(IR.NE.0) GOTO 500
        DO 490 I=MAX(1,IP),N
          IF(K(I,1).LE.0.OR.K(I,1).GT.10.OR.(I.GE.IC1.AND.I.LE.IC2
     &    .AND.K(I,1).GE.1.AND.K(I,1).LE.2)) GOTO 490
          IF(MCOMB.EQ.1) KCI=LUCOMP(K(I,2))
          IF(MCOMB.EQ.1.AND.KCI.EQ.0) GOTO 490
          IF(MCOMB.EQ.1.AND.KCHG(KCI,2).EQ.0.AND.I.LE.NS) GOTO 490
          IF(MCOMB.EQ.2.AND.IABS(K(I,2)).GT.10.AND.IABS(K(I,2)).LE.100)
     &    GOTO 490
          HCR=DPC(4)*P(I,4)-DPC(1)*P(I,1)-DPC(2)*P(I,2)-DPC(3)*P(I,3)
          HSR=2.*HCR+PECM**2-P(N+2,5)**2-2.*P(N+2,5)*P(I,5)
          IF(HSR.GT.HSM) THEN
            IR=I
            HA=HCR
            HSM=HSR
          ENDIF
  490   CONTINUE
  500 CONTINUE
 
C...Shuffle energy and momentum to put new particle on mass shell.
      IF(IR.NE.0) THEN
        HB=PECM**2+HA
        HC=P(N+2,5)**2+HA
        HD=P(IR,5)**2+HA
        HK2=0.5*(HB*SQRT(MAX(0.,((HB+HC)**2-4.*(HB+HD)*P(N+2,5)**2)/
     &  (HA**2-(PECM*P(IR,5))**2)))-(HB+HC))/(HB+HD)
        HK1=(0.5*(P(N+2,5)**2-PECM**2)+HD*HK2)/HB
        DO 510 J=1,4
          P(N+2,J)=(1.+HK1)*DPC(J)-HK2*P(IR,J)
          P(IR,J)=(1.+HK2)*P(IR,J)-HK1*DPC(J)
  510   CONTINUE
        N=N+2
      ELSE
        CALL LUERRM(3,'(LUPREP:) no match for collapsing cluster')
        RETURN
      ENDIF
 
C...Mark collapsed system and store daughter pointers. Iterate.
  520 DO 530 I=IC1,IC2
        IF((K(I,1).EQ.1.OR.K(I,1).EQ.2).AND.
     &  KCHG(LUCOMP(K(I,2)),2).NE.0) THEN
          K(I,1)=K(I,1)+10
          IF(MSTU(16).NE.2) THEN
            K(I,4)=NSAV+1
            K(I,5)=NSAV+1
          ELSE
            K(I,4)=NSAV+2
            K(I,5)=NSAV+1+NBODY
          ENDIF
        ENDIF
  530 CONTINUE
      IF(N.LT.MSTU(4)-MSTU(32)-5) GOTO 140
 
C...Check flavours and invariant masses in parton systems.
  540 NP=0
      KFN=0
      KQS=0
      DO 550 J=1,5
        DPS(J)=0.
  550 CONTINUE
      DO 580 I=MAX(1,IP),N
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 580
        KC=LUCOMP(K(I,2))
        IF(KC.EQ.0) GOTO 580
        KQ=KCHG(KC,2)*ISIGN(1,K(I,2))
        IF(KQ.EQ.0) GOTO 580
        NP=NP+1
        IF(KQ.NE.2) THEN
          KFN=KFN+1
          KQS=KQS+KQ
          MSTJ(93)=1
          DPS(5)=DPS(5)+ULMASS(K(I,2))
        ENDIF
        DO 560 J=1,4
          DPS(J)=DPS(J)+P(I,J)
  560   CONTINUE
        IF(K(I,1).EQ.1) THEN
          IF(NP.NE.1.AND.(KFN.EQ.1.OR.KFN.GE.3.OR.KQS.NE.0)) CALL
     &    LUERRM(2,'(LUPREP:) unphysical flavour combination')
          IF(NP.NE.1.AND.DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2.LT.
     &    (0.9*PARJ(32)+DPS(5))**2) THEN
            CALL LUERRM(3,'(LUPREP:) too small mass in jet system')
          END IF
          NP=0
          KFN=0
          KQS=0
          DO 570 J=1,5
            DPS(J)=0.
  570     CONTINUE
        ENDIF
  580 CONTINUE
 
      RETURN
      END
