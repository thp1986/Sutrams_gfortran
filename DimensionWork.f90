!     SUBROUTINE        D  I  M  W  R  K       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  RETURN DIMENSIONS FOR SLAP WORK ARRAYS.  THE DIMENSIONS DEPEND   
! ***  ON THE PARTICULAR SOLVER CHOSEN.                                 
!                                                                       
      SUBROUTINE DIMWRK (KSOLVR, NSAVE, NN, NELT, NWI, NWF) 
      use SutraMSPrecision
      IMPLICIT NONE
      INTEGER (I4B) :: &
        KSOLVR, NSAVE, NN, NELT, NWI, NWF
      !LOCALS
      INTEGER (I4B) :: &
        NL
!                                                                       
!.....COMPUTE DIMENSIONS.                                               
!                                                                       
      IF (KSOLVR.EQ.1) THEN 
         NL = (NELT + NN) / 2 
         NWI = 11 + 2 * NL 
         NWF = NL + 5 * NN + 1 
      ELSEIF (KSOLVR.EQ.2) THEN 
         NWI = 31 + 2 * NELT 
         NWF = 2 + NN * (NSAVE+7) + NSAVE * (NSAVE+3) + (NELT - NN) 
      ELSEIF (KSOLVR.EQ.3) THEN 
         NWI = 11 + 2 * NELT 
         NWF = 1 + 3 * NN * (NSAVE+1) + 7 * NN + NSAVE+ (NELT - NN) 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE DIMWRK                         
