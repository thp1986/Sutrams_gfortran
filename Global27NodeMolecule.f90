!     SUBROUTINE        G  L  O  2  7          SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO ASSEMBLE RESULTS OF ELEMENTWISE INTEGRATIONS INTO             
! ***  A GLOBAL "27-NODE MOLECULE"-FORMAT MATRIX AND GLOBAL VECTOR      
! ***  FOR BOTH FLOW AND TRANSPORT EQUATIONS.                           
! ***  THIS SUBROUTINE IS NOT CURRENTLY USED BY SUTRA.                  
!                                                                       
      SUBROUTINE GLO27 (L, ML, VOLE, BFLOWE, DFLOWE, BTRANE, DTRANE, IN,&
                        VOL, PMAT, PVEC, UMAT, UVEC, NBI27)                               
      USE CONTRL 
      USE DIMS
      USE DIMX
      use SutraMSPrecision
      IMPLICIT NONE
      !IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      real (DP) :: &
        BFLOWE (8, 8), DFLOWE (8), BTRANE (8, 8), DTRANE (8, 8), VOLE (8)                                                          
      real (DP) :: &
        VOL (NN), PMAT (NELT, NCBI), PVEC (NN), UMAT (NELT, NCBI), UVEC (NN, NSPE)                                            
      integer (I4B) :: &
        IN (NIN), NBI27 (NBIX) 
      INTEGER (I4B) :: &
        L, ML
      !LOCALS
      INTEGER (I4B) :: &
        N1, N8, &
        IB, IE, II, &
        JB, JE, JJ, &
        J27

!                                                                       
      N1 = (L - 1) * N48 + 1 
      N8 = N1 + N48 - 1 
!                                                                       
!.....ADD RESULTS OF INTEGRATIONS OVER ELEMENT L TO GLOBAL              
!        P-MATRIX AND P-VECTOR                                          
      IF (ML - 1) 9050, 9050, 9150 
 9050 IE = 0 
      DO 9100 II = N1, N8 
         IE = IE+1 
         IB = IN (II) 
         VOL (IB) = VOL (IB) + VOLE (IE) 
         PVEC (IB) = PVEC (IB) + DFLOWE (IE) 
         JE = 0 
         DO 9100 JJ = N1, N8 
            JE = JE+1 
            JB = IN (JJ) - IB + NBHALF 
            J27 = NBI27 (JB) 
 9100 PMAT (IB, J27) = PMAT (IB, J27) + BFLOWE (IE, JE) 
      IF (ML - 1) 9150, 9300, 9150 
!                                                                       
!.....ADD RESULTS OF INTEGRATIONS OVER ELEMENT L TO GLOBAL              
!        U-MATRIX                                                       
 9150 IF (NOUMAT.EQ.1) GOTO 9300 
      IE = 0 
      DO 9200 II = N1, N8 
         IE = IE+1 
         IB = IN (II) 
!.....POSITION FOR ADDITION TO U-VECTOR                                 
!     DO K=1,NSPE                                                       
!       UVEC(IB,K)=UVEC(IB,K)+ ((   ))                                  
!     END DO                                                            
         JE = 0 
         DO 9200 JJ = N1, N8 
            JE = JE+1 
            JB = IN (JJ) - IB + NBHALF 
            J27 = NBI27 (JB) 
 9200 UMAT (IB, J27) = UMAT (IB, J27) + DTRANE (IE, JE) + BTRANE (IE,   &
      JE)                                                               
!                                                                       
 9300 CONTINUE 
 9999 CONTINUE 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
      RETURN 
      END SUBROUTINE GLO27                          
