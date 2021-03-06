! REGE_OWNM_R.F Ales Ziberna, 2006 - ONE WAY, NORMALIZED MATRICES version of REGE (Douglas R. White, 1985)
!  THIS VERSION ALLOWS USER TO SET THE NUMBER OF ITERATIONS 
      subroutine regeownmdiag(R,B,N,NR,ITER)
      DOUBLE PRECISION   R, B, DEG, SUM, SUMM1, SUMM2, XMAX1, XMAX2, CMIKJM1, CMIKJM2
      INTEGER NR, N, ITER, JJ, II !, KR
      DIMENSION  DEG (N), SUM (N,N),R (N,N, NR), B (N,N)

  
!     COMPUTE DEGREE, SUMS FOR I--&gt;K, INITIAL STRUCTURAL EQUIV.
      DO 100 I=1,N
      DEG(I) = 0.0
      DO 100 J=1,N
      SUM(I,J)= R(I,J,1) + R(J,I,2)
  100 DEG(I)=DEG(I) + SUM(I,J)

!     BEGIN ITERATIONS
      DO 700 L=1,ITER
!     INITIALIZE DIFFERENCE IN SUCCESSIVE SE MATRICES
      D = 0.0
!     TAKE POINT I
      DO 520 II = 1, N-1
      I=II
!     IF DEGREE ZERO NEXT I
!     IF(DEG(I).EQ.0.0) GO TO 520
!     TAKE POINT J
      DO 510 JJ= II+1, N
      CM = 0.0
      J=JJ
!     IF DEGREE ZERO NEXT J
      IF((DEG(J)).EQ.0.0) GO TO 506
      I=II
! TAKE EACH OF THE TWO POINTS AS REFERENT
      DO 505 IJ=1,2
      IF (IJ.EQ.1) GOTO 120
      J=II
      I=JJ
!     TAKE POINT K (I--&gt;K, K--&gt;I TO BE EXAMINED)
 120  DO 500 K=1,N
      IF(SUM(I,K).EQ.0.0) GO TO 500
      IF (I.EQ.K) GO TO 500
      XMAX1=0.0
      XMAX2=0.0
!     FIND BEST MATCHING POINT M
      DO 400 M=1,N
      IF(SUM(J,M).EQ.0.0) GO TO 400
      IF(J.EQ.M) GO TO 400
      SUMM1=0.0
      SUMM2=0.0
!     DO 300 KR=1,NR
      SUMM1 = SUMM1 +min (R(I,K,1),r(j,m,1))
      SUMM2 = SUMM2 +min (R(K,I,2),r(m,j,2))
! 300
      CMIKJM1 = SUMM1 * b (max (k,m), min (k,m))
      CMIKJM2 = SUMM2 * b (max (k,m), min (k,m))
!     IF PERFECT MATCH DESIRED, CORRECT MATCH
!     IF(SUMM.NE.SUM(I,K).AND.NOERRS.EQ.1)  CMIKJM=0.0
      IF(CMIKJM1.GT.XMAX1) XMAX1= CMIKJM1
      IF(CMIKJM2.GT.XMAX2) XMAX2= CMIKJM2
 
      IF((XMAX1 + XMAX2).EQ.SUM(I,K)) GO TO 450
 
  400 CONTINUE
!     ADD BEST MATCHES TO REGULAR EQUIVALENCE NUMERATOR FOR I,J
  450 CM=CM+XMAX1 + XMAX2
  500  CONTINUE
      CM=CM + b (max (i,j), min (i,j))*(min(R(I,I,1),r(j,j,1))+min(R(I,I,2),r(j,j,2)))
  505  CONTINUE
!     COMPUTE REGULAR EQUIVALENCE
  506 DM = DEG(II)+DEG(JJ)
      B (II,JJ)= 1.0
      IF(DM.NE.0.0) B (II,JJ)=CM/DM
!     IF(B (II,JJ).LE.CUT) B (II,JJ)=0.0
  510  CONTINUE
  520  CONTINUE

! symmetrize : to lower half matrix
      DO 600 I = 2, N
      DO 600 J = 1, i-1
  600 B(i,j) = B(j,i) 
  700 CONTINUE

      END
