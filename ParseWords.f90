!     SUBROUTINE        P  R  S  W  D  S       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  PARSE A CHARACTER STRING INTO WORDS.  WORDS ARE CONSIDERED TO BE 
! ***  SEPARATED BY SPACES.                                             
!                                                                       
      SUBROUTINE PRSWDS (STRING, NWMAX, WORD, NWORDS) 

      USE SutraMSPrecision

      IMPLICIT NONE
      
      INTEGER (I4B) :: &
        NWMAX, NWORDS
      CHARACTER (len=80) :: WORD(NWMAX) 
      CHARACTER (len=80) :: STRING 

      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        I, &
        M, M1, M2
!                                                                       
      DO 50 I = 1, NWMAX 
         WORD (I) = "" 
   50 END DO 
!                                                                       
      NWORDS = 0 
      M2 = 1 
!                                                                       
  300 CONTINUE 
      DO 350 M = M2, 80 
         IF (STRING (M:M) .NE.' ') THEN 
            M1 = M 
            GOTO 400 
         ENDIF 
  350 END DO 
      RETURN 
!                                                                       
  400 CONTINUE 
      DO 450 M = M1 + 1, 80 
         IF (STRING (M:M) .EQ.' ') THEN 
            M2 = M 
            GOTO 500 
         ENDIF 
  450 END DO 
      M2 = 80 
!                                                                       
  500 CONTINUE 
      NWORDS = NWORDS + 1 
      WORD (NWORDS) = STRING (M1:M2) 
!                                                                       
      IF ( (M2.LT.80) .AND. (NWORDS.LT.NWMAX) ) GOTO 300 
!                                                                       
      RETURN 
      END SUBROUTINE PRSWDS                         
