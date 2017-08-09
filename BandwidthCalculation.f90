!     SUBROUTINE        B  A  N  W  I  D       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO CALCULATE AND CHECK BAND WIDTH OF FINITE ELEMENT MESH.        
!                                                                       
      SUBROUTINE BANWID (IN) 
        USE SOLVI 
        USE FUNITS 
        USE DIMS
        USE DIMX
        use SutraMSPrecision
        USE MSErrorHandler
        IMPLICIT NONE
        INTEGER (I4B) :: &
          IN (NIN) 
        !LOCALS
          INTEGER (I4B) :: &
            NBTEST, &
            NDIF, &
            I, &
            II, &
            L, &
            IELO, &
            IEHI, &
            NDIFF, &
            LEM, &
            NBL
!                                                                       
      NBTEST = 0 
      NDIF = 0 
      II = 0 
      WRITE (fLST, 100) 
  100 FORMAT(////11X,'**** MESH ANALYSIS ****'//) 
!                                                                       
!.....FIND ELEMENT WITH MAXIMUM DIFFERENCE IN NODE NUMBERS              
      DO 2000 L = 1, NE 
         II = II + 1 
         IELO = IN (II) 
         IEHI = IN (II) 
         DO 1000 I = 2, N48 
            II = II + 1 
            IF (IN (II) .LT.IELO) IELO = IN (II) 
 1000    IF (IN (II) .GT.IEHI) IEHI = IN (II) 
         NDIFF = IEHI - IELO 
         IF (NDIFF.GT.NDIF) THEN 
            NDIF = NDIFF 
            LEM = L 
         ENDIF 
         NBL = 2 * NDIFF + 1 
 2000 END DO 
!                                                                       
!.....CALCULATE ACTUAL FULL BAND WIDTH, NB.                             
      NB = 2 * NDIF + 1 
      NBHALF = NDIF + 1 

      !SET PARAMETERS THAT DEPEND ON BANDWIDTH
      NBI  = NB
      NBIX = NB

      WRITE (fLST, 2500) NB, LEM, NBI 
 2500 FORMAT(//13X,'ACTUAL MAXIMUM FULL BANDWIDTH, ',I5,                &
     &   ', WAS CALCULATED IN ELEMENT ',I9/13X,7(1H-),                  &
     &   'INPUT FULL BANDWIDTH IS ',I5)                                 
      IF (NBTEST.EQ.0) GOTO 3000 

 3000 WRITE (fLST, 4000) 
 4000 FORMAT(////////132(1H-)///42X,'E N D   O F   I N P U T   ',       &
                                    'F R O M   U N I T - 5'//132(1H-))                             
      RETURN 
      END SUBROUTINE BANWID                         
