!                                                                       
!     SUBROUTINE        T  R  I  S  E  T       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO SET UP THE POINTER ARRAYS THAT GIVE THE MATRIX STRUCTURE      
! ***  IN "SLAP TRIAD" FORMAT.                                          
!                                                                       
      SUBROUTINE TRISET (ITRI, JTRI) 

      USE CONTRL 
      USE SOLVI 
      USE DIMS
      USE DIMX
      USE SutraMSPrecision

      IMPLICIT NONE

      INTEGER (I4B) :: &
        ITRI (NELT), JTRI (NELT) 
      
      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        IS, IBEG, IEND, &
        JS, &
        KS, &
        M, &
        N, NB1, NBM, NST, NBH1, NBMC, NBMJ, NBMK, NN12, NSTJ, NSTK, NNNBH, NBMNBH

!                                                                       
!.....DEFINE CERTAIN QUANTITIES FOR CONVENIENCE AND EFFICIENCY.         
!                                                                       
      NN12 = NN1 * NN2 
      NNNBH = NN - NBHALF 
      NBH1 = NBHALF + 1 
      NB1 = NB + 1 
!                                                                       
!.....CREATE THE POINTER ARRAYS "ITRI" AND "JTRI", WHICH SPECIFY THE    
!        MATRIX ARRAY STRUCTURE IN "SLAP TRIAD" FORMAT.                 
!                                                                       
      IF (IABS (KTYPE) .EQ.3) THEN 
!.....3-D PROBLEM.                                                      
!                                                                       
         M = 0 
         DO 400 KS = 0, 2 
            NBMK = KS * NN12 
            NSTK = KS * 9 
            DO 400 JS = 0, 2 
               NBMJ = NBMK + JS * NN1 
               NSTJ = NSTK + JS * 3 
               DO 400 IS = 1, 3 
                  NBM = NBMJ + IS 
                  NBMC = NB1 - NBM 
                  NST = NSTJ + IS 
                  IF (NST.LT.14) THEN 
                     IBEG = NBH1 - NBM 
                     IEND = NN 
                  ELSE 
                     IBEG = 1 
                     IEND = NNNBH + NBMC 
                  ENDIF 
                  NBMNBH = NBM - NBHALF 
                  DO 300 N = IBEG, IEND 
                     M = M + 1 
                     ITRI (M) = N 
                     JTRI (M) = N + NBMNBH 
  300             END DO 
  400    CONTINUE 
!                                                                       
      ELSE 
!.....2-D PROBLEM.                                                      
!                                                                       
         M = 0 
         DO 1400 JS = 0, 2 
            NBMJ = JS * NN1 
            NSTJ = JS * 3 
            DO 1400 IS = 1, 3 
               NBM = NBMJ + IS 
               NBMC = NB1 - NBM 
               NST = NSTJ + IS 
               IF (NST.LT.5) THEN 
                  IBEG = NBH1 - NBM 
                  IEND = NN 
               ELSE 
                  IBEG = 1 
                  IEND = NNNBH + NBMC 
               ENDIF 
               NBMNBH = NBM - NBHALF 
               DO 1300 N = IBEG, IEND 
                  M = M + 1 
                  ITRI (M) = N 
                  JTRI (M) = N + NBMNBH 
 1300          END DO 
 1400    CONTINUE 
!                                                                       
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE TRISET                         
