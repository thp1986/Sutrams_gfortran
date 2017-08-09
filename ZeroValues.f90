!                                                                       
!     SUBROUTINE        Z  E  R  O             SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO FILL AN ARRAY WITH A CONSTANT VALUE.                          
!                                                                       
      SUBROUTINE ZERO (A, IADIM, FILL) 

      USE SutraMSPrecision

      IMPLICIT NONE

      INTEGER (I4B) :: &
        IADIM
      REAL (DP) :: &
        FILL
      REAL (DP) :: &
        A (IADIM) 

      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        I
!                                                                       
!.....FILL ARRAY A WITH VALUE IN VARIABLE 'FILL'                        
!      FORALL (I = 1:IADIM) 
!        A (I) = FILL 
!      ENDFORALL 
      DO I = 1, IADIM
        A(I) = FILL 
      END DO
!                                                                       
!.....RETURN TO CALLING ROUTINE                                         
      RETURN 
      END SUBROUTINE ZERO                           
