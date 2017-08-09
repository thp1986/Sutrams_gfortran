!     SUBROUTINE        S  T  O  R  E          SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO STORE RESULTS THAT MAY LATER BE USED TO RE-START              
! ***  THE SIMULATION.                                                  
!                                                                       
      SUBROUTINE STORE (PVEC, UVEC, PM1, UM1, CS1, RCIT, SW, PBC) 
      
      USE FUNITS 
      USE DIMS
      USE TIMES
      USE SutraStorage, ONLY : SpecifiedPBC
      USE SutraMSPrecision

      IMPLICIT NONE
      
      REAL (DP) :: &
        PVEC (NN), UVEC (NN, NSPE), PM1 (NN), UM1 (NN, NSPE),   &
        CS1 (NN, NSPE), RCIT (NN), SW (NN), PBC (NBCN)                    

      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        I, IP, &
        K
!
!.....CHECK THAT fRST IS OPEN
      IF(fRST<1) GOTO 9999
!                                                                       
!.....REWIND UNIT-fRST FOR WRITING RESULTS OF CURRENT TIME STEP           
      REWIND (fRST) 
!                                                                       
!.....STORE TIME INFORMATION                                            
      WRITE (fRST, 100) TSEC, DELTP, DELTU 
  100 FORMAT(4D20.10) 
!                                                                       
!.....STORE SOLUTION                                                    
      WRITE (fRST, 105) 
      WRITE (fRST, 110) (PVEC (I), I = 1, NN) 
      WRITE (fRST, 105) 
      DO K = 1, NSPE 
      WRITE (fRST, 110) (UVEC (I, K), I = 1, NN) 
      ENDDO 
      WRITE (fRST, 110) (PM1 (I), I = 1, NN) 
      DO K = 1, NSPE 
      WRITE (fRST, 110) (UM1 (I, K), I = 1, NN) 
      ENDDO 
      DO K = 1, NSPE 
      WRITE (fRST, 110) (CS1 (I, K), I = 1, NN) 
      ENDDO 
      WRITE (fRST, 110) (RCIT (I), I = 1, NN) 
      WRITE (fRST, 110) (SW (I), I = 1, NN) 
      WRITE (fRST, 110) (SpecifiedPBC(IP)%P, IP = 1, NPBC) 
  105 FORMAT("'NONUNIFORM'") 
  110 FORMAT(1PD20.13,1X,1PD20.13,1X,1PD20.13,1X,1PD20.13) 
                                                                       
                                                                       
 9999 RETURN 
      END SUBROUTINE STORE                          
