!                                                                       
!     SUBROUTINE        S  O  L  V  E  B       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO SOLVE THE MATRIX EQUATION BY:                                 
! ***   (1) DECOMPOSING THE MATRIX                                      
! ***   (2) MODIFYING THE RIGHT-HAND SIDE                               
! ***   (3) BACK-SUBSTITUTING FOR THE SOLUTION                          
!                                                                       

      SUBROUTINE SOLVEB (KKK, C, R, NNP, IHALFB, MAXNP, MAXBW) 

      USE SutraMSPrecision

      IMPLICIT NONE
      
      INTEGER (I4B) :: &
        KKK, NNP, IHALFB, &
        MAXNP, MAXBW
      REAL (DP) :: &
        C (MAXNP, MAXBW), R (MAXNP) 

      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        IB, IHBP, IBAND, &
        JB, &
        KB, &
        LB, &
        MB, &
        NB, NI, NJ, NL, NK, NR, NU
      REAL (DP) :: &
        A, PIVOTI, SUM

      IHBP = IHALFB + 1 
!                                                                       
!.....DECOMPOSE MATRIX C BY BANDED GAUSSIAN ELIMINATION FOR             
!        NON-SYMMETRIC MATRIX                                           
      IF (KKK - 1) 5, 5, 50 
    5 NU = NNP - IHALFB 
      DO 20 NI = 1, NU 
         PIVOTI = 1.D0 / C (NI, IHBP) 
         NJ = NI + 1 
         IB = IHBP 
         NK = NI + IHALFB 
         DO 10 NL = NJ, NK 
            IB = IB - 1 
            A = - C (NL, IB) * PIVOTI 
            C (NL, IB) = A 
            JB = IB + 1 
            KB = IB + IHALFB 
            LB = IHBP - IB 
            DO 10 MB = JB, KB 
               NB = LB + MB 
   10    C (NL, MB) = C (NL, MB) + A * C (NI, NB) 
   20 END DO 
      NR = NU + 1 
      NU = NNP - 1 
      NK = NNP 
      DO 40 NI = NR, NU 
         PIVOTI = 1.D0 / (C (NI, IHBP) ) 
         NJ = NI + 1 
         IB = IHBP 
         DO 30 NL = NJ, NK 
            IB = IB - 1 
            A = - C (NL, IB) * PIVOTI 
            C (NL, IB) = A 
            JB = IB + 1 
            KB = IB + IHALFB 
            LB = IHBP - IB 
            DO 30 MB = JB, KB 
               NB = LB + MB 
   30    C (NL, MB) = C (NL, MB) + A * C (NI, NB) 
   40 END DO 
      IF (KKK - 1) 50, 44, 50 
   44 RETURN 
!                                                                       
!.....UPDATE RIGHT-HAND SIDE VECTOR, R                                  
   50 NU = NNP + 1 
      IBAND = 2 * IHALFB + 1 
      DO 70 NI = 2, IHBP 
         IB = IHBP - NI + 1 
         NJ = 1 
         SUM = 0.0D0 
         DO 60 JB = IB, IHALFB 
            SUM = SUM + C (NI, JB) * R (NJ) 
   60    NJ = NJ + 1 
   70 R (NI) = R (NI) + SUM 
      IB = 1 
      NL = IHBP + 1 
      DO 90 NI = NL, NNP 
         NJ = NI - IHBP + 1 
         SUM = 0.D0 
         DO 80 JB = IB, IHALFB 
            SUM = SUM + C (NI, JB) * R (NJ) 
   80    NJ = NJ + 1 
   90 R (NI) = R (NI) + SUM 
!                                                                       
!.....BACK SOLVE                                                        
      R (NNP) = R (NNP) / C (NNP, IHBP) 
      DO 110 IB = 2, IHBP 
         NI = NU - IB 
         NJ = NI 
         MB = IHALFB + IB 
         SUM = 0.D0 
         DO 100 JB = NL, MB 
            NJ = NJ + 1 
  100    SUM = SUM + C (NI, JB) * R (NJ) 
  110 R (NI) = (R (NI) - SUM) / C (NI, IHBP) 
      MB = IBAND 
      DO 130 IB = NL, NNP 
         NI = NU - IB 
         NJ = NI 
         SUM = 0.D0 
         DO 120 JB = NL, MB 
            NJ = NJ + 1 
  120    SUM = SUM + C (NI, JB) * R (NJ) 
  130 R (NI) = (R (NI) - SUM) / C (NI, IHBP) 
!                                                                       
!                                                                       
      RETURN 
      END SUBROUTINE SOLVEB                         
