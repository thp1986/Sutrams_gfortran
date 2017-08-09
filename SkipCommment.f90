!                                                                       
!     SUBROUTINE        S  K  P  C  O  M       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  IDENTIFY AND SKIP OVER COMMENT LINES IN THE UNIT fINP INPUT FILE   
! ***  AND RETURN THE NUMBER OF LINES SKIPPED.                          
!                                                                       

      SUBROUTINE SKPCOM (KU, NLSKIP) 

      USE MSErrorHandler
      USE SutraMSPrecision

      IMPLICIT NONE

      INTEGER (I4B) :: &
        KU, NLSKIP

      !LOCAL VARIABLES
      CHARACTER (LEN=1) :: CDUM

!                                                                       
      NLSKIP = 0 
  100 READ (KU, 111,iostat=ios) CDUM 
      if(ios/=0) call ErrorIO('Comment Line Read')
  111 FORMAT (A1) 
      IF (CDUM.EQ.'#') THEN 
         NLSKIP = NLSKIP + 1 
         GOTO 100 
      ENDIF 
!                                                                       
      BACKSPACE (KU) 
!                                                                       
      RETURN 
      END SUBROUTINE SKPCOM                         
