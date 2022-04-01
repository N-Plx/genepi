 

C*********************************************************************
 
C...LUDCYK
C...Handles flavour production in the decay of unstable particles
C...and small string clusters.
 
      SUBROUTINE LUDCYK(KFL1,KFL2,KFL3,KF)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(D)
      IMPLICIT INTEGER(I-N)
      INTEGER LUK,LUCHGE,LUCOMP
C...Commonblocks.
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT1/,/LUDAT2/
 
 
C.. Call LUKFDI directly if no popcorn option is on
      IF(MSTJ(12).LT.2) THEN
         CALL LUKFDI(KFL1,KFL2,KFL3,KF)
         MSTU(124)=KFL3
         RETURN
      ENDIF
 
      KFL3=0
      KF=0
      IF(KFL1.EQ.0) RETURN
      KF1A=IABS(KFL1)
      KF2A=IABS(KFL2)
 
      NSTO=130
      NMAX=MIN(MSTU(125),10)
 
C.. Identify rank 0 cluster qq
      IRANK=1
      IF(KF1A.GT.10.AND.KF1A.LT.10000) IRANK=0
 
      IF(KF2A.GT.0)THEN
C.. Join jets: Fails if store not empty
         IF(MSTU(121).GT.0) THEN
            MSTU(121)=0
            RETURN
         ENDIF
         CALL LUKFDI(KFL1,KFL2,KFL3,KF)
      ELSEIF(KF1A.GT.10.AND.MSTU(121).GT.0)THEN
C.. Pick popcorn meson from store, return same qq, decrease store
         KF=MSTU(NSTO+MSTU(121))
         KFL3=-KFL1
         MSTU(121)=MSTU(121)-1
      ELSE
C.. Generate new flavour. Then done if no diquark is generated
  100    CALL LUKFDI(KFL1,0,KFL3,KF)
         IF(MSTU(121).EQ.-1) GOTO 100
         MSTU(124)=KFL3
         IF(KF.EQ.0.OR.IABS(KFL3).LE.10) RETURN
 
C.. Simple case if no dynamical popcorn suppressions are considered
         IF(MSTJ(12).LT.4) THEN
            IF(MSTU(121).EQ.0) RETURN
            NMES=1
            KFPREV=-KFL3
            CALL LUKFDI(KFPREV,0,KFL3,KFM)
C.. Due to eta+eta' suppr., a qq->M+qq attempt might end as qq->B+q
            IF(IABS(KFL3).LE.10)THEN
               KFL3=-KFPREV
               RETURN
            ENDIF
            GOTO 120
         ENDIF
 
C test output qq against fake Gamma, then return if no popcorn.
         GB=2.
         IF(IRANK.NE.0)THEN
            CALL LUZDIS(1,2103,5.,Z)
            GB=5.*(1.-Z)/Z
            IF(1.-PARF(192)**GB.LT.RLU(0)) THEN
               MSTU(121)=0
               GOTO 100
            ENDIF
         ENDIF
         IF(MSTU(121).EQ.0) RETURN
 
C..Set store size memory. Pick fake dynamical variables of qq.
         NMES=MSTU(121)
         CALL LUPTDI(1,PX3,PY3)
         X=1.
         POPM=0.
         G=GB
         POPG=GB
 
C.. Pick next popcorn meson, test with fake dynamical variables
  110    KFPREV=-KFL3
         PX1=-PX3
         PY1=-PY3
         CALL LUKFDI(KFPREV,0,KFL3,KFM)
         IF(MSTU(121).EQ.-1) GOTO 100
         CALL LUPTDI(KFL3,PX3,PY3)
         PM=ULMASS(KFM)**2+(PX1+PX3)**2+(PY1+PY3)**2
         CALL LUZDIS(KFPREV,KFL3,PM,Z)
         G=(1.-Z)*(G+PM/Z)
         X=(1.-Z)*X
 
         PTST=1.
         GTST=1.
         RTST=RLU(0)
         IF(MSTJ(12).GT.4)THEN
            POPMN=SQRT((1.-X)*(G/X-GB))
            POPM=POPM+PMAS(LUCOMP(KFM),1)-PMAS(LUCOMP(KFM),3)
            PTST=EXP((POPM-POPMN)*PARF(193))
            POPM=POPMN
         ENDIF
         IF(IRANK.NE.0)THEN
            POPGN=X*GB
            GTST=(1.-PARF(192)**POPGN)/(1.-PARF(192)**POPG)
            POPG=POPGN
         ENDIF
         IF(RTST.GT.PTST*GTST)THEN
            MSTU(121)=0
            IF(RTST.GT.PTST) MSTU(121)=-1
            GOTO 100
         ENDIF
 
C.. Store meson
  120    IF(NMES.LE.NMAX) MSTU(NSTO+MSTU(121)+1)=KFM
         IF(MSTU(121).GT.0) GOTO 110
 
C.. Test accepted system size. If OK set global popcorn size variable.
         IF(NMES.GT.NMAX)THEN
            KF=0
            KFL3=0
            RETURN
         ENDIF
         MSTU(121)=NMES
      ENDIF
 
      RETURN
      END
