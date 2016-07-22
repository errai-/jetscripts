 
C*********************************************************************
 
C...PYEVWT
C...Dummy routine, which the user can replace in order to multiply the
C...standard PYTHIA differential cross-section by a process- and
C...kinematics-dependent factor WTXS. For MSTP(142)=1 this corresponds
C...to generation of weighted events, with weight 1/WTXS, while for
C...MSTP(142)=2 it corresponds to a modification of the underlying
C...physics.
 
      SUBROUTINE PYEVWT(WTXS)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      SAVE /PYDAT1/,/PYINT1/,/PYINT2/,/PYSUBS/
C... CSA specific 
      integer CSAMODE
      integer pad
      double precision  MUONRW, GAMMAJRW, ZJRW, ZPRW, HLTRW,
     &  SUSYRW, WWRW, PTPOWER
      common /EXPAR/ pad, CSAMODE,MUONRW, GAMMAJRW, ZJRW, ZPRW, 
     &  HLTRW, SUSYRW, WWRW, PTPOWER
 
C...Set default weight for WTXS.

C      write(*,*) 'Reweighting event ...'
      WTXS=1D0
 
C...Read out subprocess number.
      ISUB=MINT(1)
      ISTSB=ISET(ISUB)
 
C...Read out tau, y*, cos(theta), tau' (where defined, else =0).
      TAU=VINT(21)
      YST=VINT(22)
      CTH=0D0
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) CTH=VINT(23)
      TAUP=0D0
      IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUP=VINT(26)
 
C...Read out x_1, x_2, x_F, shat, that, uhat, p_T^2.
      X1=VINT(41)
      X2=VINT(42)
      XF=X1-X2
      SHAT=VINT(44)
      THAT=VINT(45)
      UHAT=VINT(46)
      PT2=VINT(48)

      PTHAT = SQRT(PT2)  
      PTPOWER = 4.5D0
      WTXS=(PTHAT/CKIN(3))**PTPOWER

      RETURN
      END
