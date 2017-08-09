!                                                                       
!     SUBROUTINE        S  O  L  V  E  R       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO CALL THE APPROPRIATE SOLVER.                                  
!                                                                       
      LOGICAL FUNCTION SOLVER (KKK, KPU, KSOLVR, C, R, XITER, B, NNP, IHALFB, &
                               MAXNP, MAXBW)
        USE SOLVI 
        USE PARAMS 
        USE DIMX
        USE SutraMSPrecision
        USE SutraStorage, ONLY : IWK, FWK, IPARM, FPARM, ITRI, JTRI
        
        IMPLICIT NONE

        INTEGER (I4B) :: &
          NNP, &
          MAXNP, MAXBW, &
          IHALFB, &
          KKK, KPU, KSOLVR
        REAL (DP) :: &
          C (MAXNP, MAXBW), R (NNP), XITER (NNP)
        REAL (DP) :: &
          B (NNNX) 
        
        !LOCAL VARIABLES
        LOGICAL :: SLPWRP
        LOGICAL :: LOK

!
        LOK=.true.
!                                                                       
!.......IF KSOLVR=0, CALL BANDED GAUSSIAN (DIRECT) SOLVER.                
!        OTHERWISE, CALL ITERATIVE "SLAP" SOLVER.                       
!                                                                       
        IF (KSOLVR.EQ.0) THEN 
           CALL SOLVEB (KKK, C, R, NNP, IHALFB, MAXNP, MAXBW) 
        ELSEIF (KSOLVR>0.and.KSOLVR<4) then 
           if(.not.SLPWRP(KKK, KPU, KSOLVR, C, R, XITER, B, NNP, IHALFB,    &
                          MAXNP, MAXBW)) LOK=.false.              
        !Additional solvers
        ELSE 
           LOK=.false.
        ENDIF 
!                                                                       
09999&
        SOLVER=LOK
        RETURN 
      END FUNCTION SOLVER                         
