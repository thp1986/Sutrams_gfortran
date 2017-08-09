!     SUBROUTINE        O  B  S  E  R  V       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO MAKE OBSERVATIONS ON PARTICULAR TIME STEPS                    
!                                                                       
      SUBROUTINE OBSERV (IOBS, ITOBS, POBS, UOBS, OBSTIM, PVEC, UVEC, ISTOP)                                                            
      USE CONTRL 
      USE FUNITS 
      USE DIMS
      USE TIMES
      USE OBS
      USE SutraMSPrecision
      IMPLICIT NONE
      CHARACTER(13) UNAME (2) 
      CHARACTER(14) UNDERS 
      INTEGER (I4B) :: &
        ISTOP
      integer (I4B) :: &
        INOB (16) 
      integer (I4B) :: &
        IOBS (NOBSN)
      real (DP) :: &
        POBS (NOBSN, NTOBSN), UOBS (NOBSN, NTOBSN, NSPE), &
        OBSTIM (NTOBSN)
      integer (I4B) :: &
        ITOBS (NTOBSN)
      real (DP) :: &
        PVEC (NN), UVEC (NN, NSPE)
      !LOCALS
      INTEGER (I4B) :: &
        NOBS, NTOBS, &
        I, JJ, K
      !DATA
      DATA UNAME (1)  / 'CONCENTRATION' / , UNAME (2)  / '  TEMPERATURE'  / , &
           UNDERS / '______________' /                                  
      SAVE UNAME, UNDERS 
!                                                                       
!.....NOBS IS ACTUAL NUMBER OF OBSERVATION NODES                        
!.....NTOBS IS MAXIMUM NUMBER OF TIME STEPS WITH OBSERVATIONS           
      NOBS = NOBSN - 1 
      NTOBS = NTOBSN - 2 
!                                                                       
!                                                                       
!.....MAKE OBSERVATIONS EACH NOBCYC TIME STEPS                          
  500 CONTINUE 
      IF (MOD (IT, NOBCYC) .NE.0.AND.IT.GT.1.AND.ISTOP.EQ.0) RETURN 
      IF (IT.EQ.0) RETURN 
      ITCNT = ITCNT + 1 
      ITOBS (ITCNT) = IT 
      OBSTIM (ITCNT) = TSEC 
      DO 1000 JJ = 1, NOBS 
         I = IOBS (JJ) 
         POBS (JJ, ITCNT) = PVEC (I) 
         DO K = 1, NSPE 
           UOBS (JJ, ITCNT, K) = UVEC (I, K) 
         ENDDO 
 1000 END DO 
      RETURN 
!                                                                       
!                                                                       
      END SUBROUTINE OBSERV                         
