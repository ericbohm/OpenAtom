      PROGRAM PIMD_MAKE_PSEUDO
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(P=0.3614D0,
     #   E1 = 0.2041422096422003D0, E2 = 0.1997535956961481D0,
     #   E3 = 0.2213176596405576D0, E4 = 0.03360430734640255D0,
     #   E5 = 0.4732592578721755D0, E6 =-0.509078520069735D0,
     #   E7 = 0.6772631491947646D0, E8 =-0.369912979092217D0,
     #   E9 = 0.06965131976970335D0,
     #  DE1 = 1.0D0*E1, DE2 = 2.0D0*E2, DE3 = 3.0D0*E3,
     #  DE4 = 4.0D0*E4, DE5 = 5.0D0*E5, DE6 = 6.0D0*E6,
     #  DE7 = 7.0D0*E7, DE8 = 8.0D0*E8, DE9 = 9.0D0*E9)
      CHARACTER*25 PIMD_VPS_FILE
      CHARACTER*25 LGRID_FILE(3)
      COMMON /CSPLIN/ RI(5,1000),VI(5,1000),NPTS(5)
C================================================================
C I) READ THE INSTRUCTION FILE
      OPEN(UNIT=10,FILE='PI_MD.MAKE',STATUS='OLD')
       READ(10,*)NR,RMAX,NANG
       READ(10,*)ZZZ
       READ(10,*)ALPC1,CC1
       READ(10,*)ALPC2,CC2
       READ(10,'(A)')PIMD_VPS_FILE
       DO IANG = 1,NANG+1
        READ(10,'(A)')LGRID_FILE(IANG)
       ENDDO
      CLOSE(10)
C================================================================
C II) CONSTRUCT THE PIMD_VPS_FILE
      OPEN(UNIT=10,FILE=PIMD_VPS_FILE,STATUS='NEW')
C----------------------------------------------------------------
C  A) WRITE THE HEADER
       WRITE(10,*)NR,RMAX,NANG
       ALPC1 = SQRT(ALPC1)
       ALPC2 = SQRT(ALPC2)
       CC1   = CC1*ZZZ
       CC2   = CC2*ZZZ
       ZERO  = 0.D0
       ONE   = 1.D0
       WRITE(10,*)CC1,ALPC1,CC2,ALPC2
       WRITE(10,*)ZERO,ONE
C--------------------------------------------------------------
C B) OPEN THE LOGARITHMIC GRID FILE AND SPLINE THE DATA UP
C    R*PSI_L(R) IS STORED IN VI(2,I)
       DO 100 IANG = 1,NANG+1
        OPEN(UNIT=25,FILE=LGRID_FILE(IANG),STATUS='OLD')
         READ(25,*)NPTS(1)
         NPTS(1) = NPTS(1) + 1
         NPTS(2) = NPTS(1)
         RI(1,1) = 0.0
         RI(2,1) = 0.0
         VI(1,1) = 0.0
         VI(2,1) = 0.0
         DO I = 2,NPTS(1)
          READ(25,*)RI(1,I),VI(1,I),VI(2,I)
          RI(2,I) = RI(1,I)
         ENDDO
C         VI(1,1) = VI(1,2)
        CLOSE(25)
        IF(RI(1,NPTS(1)).LT.RMAX)THEN
          WRITE(*,*)'DISTANCE RANGE ERROR INPUT FILE'
          WRITE(*,*)'RMAX = ',RMAX,' R(NPTS) = ',RI(1,NPTS(1))
          STOP
        ENDIF
        CALL SPLINE(1)
        CALL SPLINE(2)
C--------------------------------------------------------------
C C) PUT THE STUFF ON A EVEN SPACED GRID AND SUBTRACT OUT THE 
C    COULOMB INTERACTION FROM THE POTENTIAL ENERGY 
C    AND WRITE THE RESULTS TO THE PIMD_VPS FILE
        PI = DACOS(-1.D0)
        DR   = RMAX/DFLOAT(NR)
        DO I = 1,NR
         R = DFLOAT(I-1)*DR
         RALP1 = R*ALPC1
         EEE    = EXP(-RALP1*RALP1)
         TT     = 1.0/(1.0+P*RALP1)
         ERF1 = 1.0d0-((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
     #                        +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee
         RALP2 = R*ALPC2
         EEE    = EXP(-RALP2*RALP2)
         TT     = 1.0/(1.0+P*RALP2)
         ERF2 = 1.0d0-((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
     #                        +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee
         IF(I.EQ.1) THEN
           VCOR = -(CC1 + CC2)
           VTOT = VI(1,2)
         ELSE
           VCOR  = -(CC1*ERF1 + CC2*ERF2)/R
           VTOT  = SP(1,R)
         ENDIF
         DVION = vtot-vcor
c         DVION = vtot
         PSI   = SP(2,R)
         WRITE(10,*)DVION,PSI
        ENDDO
 100  CONTINUE
C===============================================================
C CLOSE THE FILE AND YOUR DONE.
      CLOSE(10)
      CALL EXIT
      END
C
      SUBROUTINE SPLINE (IP)
C
C      PEFORMS A SPLINE INTERPOLATION OF 5SERIES OF POINTS IP=1,..5 
C      XI (IP,I) = X(I)  I=1,...NPTS(IP)
C      YI (IP,I) = Y(I)  I=1,...NPTS(IP)
C      INTERPOLATION AT POINT XX IS OBTAINED BY CALLING THE 
C      FUNCTION SP(XX)
C
C     IN MAIN PROGRAM INCLUDE 1ST COMMON BLOCK
C     NPTS(1) = 100 
C     WRITE X VALUES IN INCREASING ORDER INTO XI(1,I)
C     CORRESPONDING Y VALUES INTO Y(1,I)
C     CALL SPLINE(1)
C     SP(1,X) = Y_INTERP
C     SP1(1,X) = Y"_INTERP
C     SP2(1,X) = Y""_INTERP
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /CSPLIN/ XI(5,1000),YI(5,1000),NPTS(5)
      COMMON /CSPSAV/ C(5,3,1000)
      DIMENSION D(1000),DIAG(1000)
      N = NPTS(IP)-1
C     APPROXIMATE INITIAL AND FINAL DERIVATIVES
      C(IP,1,1)   = (YI(IP,2)-YI(IP,1))/(XI(IP,2)-XI(IP,1))
      C(IP,1,N+1) = (YI(IP,N+1)-YI(IP,N))/(XI(IP,N+1)-XI(IP,N))
      C0 = 0.0D00
      C1 = 1.0D00
      C2 = 2.0D00
      C3 = 3.0D00
      DIAG(1) = C1
      D(1)    = C0
      DO 10 M=2,N+1
         MM1 = M-1
         D(M) = XI(IP,M)-XI(IP,MM1)
         DIAG(M) = (YI(IP,M)-YI(IP,MM1))/D(M)
   10 CONTINUE
      DO 20 M=2,N
        MP1 = M+1
        C(IP,1,M) = C3*(D(M)*DIAG(MP1)+D(MP1)*DIAG(M))
        DIAG(M) = C2*(D(M)+D(MP1))
   20 CONTINUE
      DO 30 M=2,N
         MP1 = M+1
         MM1 = M-1
         G = -D(MP1)/DIAG(MM1)
         DIAG(M) = DIAG(M)+G*D(MM1)
         C(IP,1,M) = C(IP,1,M)+G*C(IP,1,MM1)
   30 CONTINUE
      DO 40 M=N,2,-1
         MP1 = M+1
         C(IP,1,M) = (C(IP,1,M)-D(M)*C(IP,1,MP1))/DIAG(M)
   40 CONTINUE
C     CALCULATE ALL OTHER COEFFICIENTS
      DO 50 I=1,N
         IP1 = I+1
         DX = XI(IP,IP1)-XI(IP,I)
         DIVDF1 = (YI(IP,IP1)-YI(IP,I))/DX
         DIVDF3 = C(IP,1,I)+C(IP,1,IP1)-C2*DIVDF1
         C(IP,2,I) = (DIVDF1-C(IP,1,I)-DIVDF3)/DX
         C(IP,3,I) = DIVDF3/(DX*DX)
   50 CONTINUE
      RETURN
      END
C
      double precision FUNCTION SP (IP,XBAR)
C
C      GET THE SPLINE INTERPOLATION
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /CSPLIN/ XI(5,1000),YI(5,1000),NPTS(5)
      COMMON /CSPSAV/ C(5,3,1000)
      N = NPTS(IP)-1
      DO 10 J=1,N
         JP1 = J+1
         IF (XBAR.LT.XI(IP,JP1)) GOTO 30
   10 CONTINUE
      J = N
   30 DX = XBAR-XI(IP,J)
      SP = 
     #   YI(IP,J)+DX*(C(IP,1,J)+DX*(C(IP,2,J)+DX*C(IP,3,J)))
      RETURN
      END

