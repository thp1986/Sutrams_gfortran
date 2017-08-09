!     SUBROUTINE        G  L  O  C  O  L           SUTRA VERSION 2D3D.1  
!                                                                        
! *** PURPOSE :                                                          
! ***  TO ASSEMBLE RESULTS OF ELEMENTWISE INTEGRATIONS INTO              
! ***  A GLOBAL "SLAP COLUMN"-FORMAT MATRIX AND GLOBAL VECTOR
! ***  FOR BOTH FLOW AND TRANSPORT EQUATIONS.                            
!                                                                        
!      SUBROUTINE GLOCOL(L,ML,VOLE,BFLOWE,DFLOWE,BTRANE,DTRANE, &
!                        IN,VOL,PMAT,PVEC,UMAT,UVEC,NBI27,MIOFF,IA,JA)
      SUBROUTINE GLOCOL(L,ML,VOLE,BFLOWE,DFLOWE,BTRANE,DTRANE)
      USE CONTRL
      USE DIMS
      USE DIMX
      use SutraMSPrecision
      USE MSErrorHandler
      USE SutraStorage, ONLY : IN,VOL,PMAT,PVEC,UMAT,UVEC,NBI27,MIOFF,IA=>ITRI,JA=>JTRI
      implicit none
      real (DP)     :: BFLOWE(8,8),DFLOWE(8),BTRANE(8,8),DTRANE(8,8),VOLE(8)    
      !locals
      integer (I4B) :: &
        IB, IE, II, &
        JB, JE, JJ, &
        L, &
        M, MBEG, MEND, ML, MM, &
        N1, N8
             
      N1=(L-1)*N48+1                                                     
      N8=N1+N48-1                                                        
!                                                                        
!.....ADD RESULTS OF INTEGRATIONS OVER ELEMENT L TO GLOBAL               
!        P-MATRIX AND P-VECTOR                                           
      IF(ML-1) 9050,9050,9150                                            
 9050 IE=0                                                               
      DO 9100 II=N1,N8                                                   
      IE=IE+1                                                            
      IB=IN(II)                                                          
      VOL(IB)=VOL(IB)+VOLE(IE)                                           
      PVEC(IB)=PVEC(IB)+DFLOWE(IE)                                       
      JE=0                                                               
      DO 9100 JJ=N1,N8                                                   
      JE=JE+1                                                            
      JB = IN(JJ)                                                        
      MBEG = JA(JB)
      MEND = JA(JB + 1) - 1
      DO 9060 MM=MBEG,MEND                  !gm kluge brute force (use bisection later)
         IF (IB.EQ.IA(MM)) THEN
            M = MM
            GOTO 9100
         END IF
 9060 CONTINUE
         call ErrorIO('GloCol:: Could not find match for '//Val2Char(ib)//' or '//Val2Char(jb)//' in P Assembly')
!         print *, 'Problem -- match not found for ', ib, jb !gm kluge
 9100 PMAT(M,1)=PMAT(M,1)+BFLOWE(IE,JE)                                  
      IF(ML-1) 9150,9300,9150                                            
!                                                                        
!.....ADD RESULTS OF INTEGRATIONS OVER ELEMENT L TO GLOBAL               
!        U-MATRIX                                                        
 9150 IF(NOUMAT.EQ.1) GOTO 9300                                          
      IE=0                                                               
      DO 9200 II=N1,N8                                                   
      IE=IE+1                                                            
      IB=IN(II)                                                          
!.....POSITION FOR ADDITION TO U-VECTOR                                  
!        UVEC(IB)=UVEC(IB)+ ((   ))                                      
      JE=0                                                               
      DO 9200 JJ=N1,N8                                                   
      JE=JE+1                                                            
      JB = IN(JJ)                                                        !gm
      MBEG = JA(JB)
      MEND = JA(JB + 1) - 1
      DO 9160 MM=MBEG,MEND                  !gm kluge brute force (use bisection later)
         IF (IB.EQ.IA(MM)) THEN
            M = MM
            GOTO 9200
         END IF
 9160 CONTINUE
         call ErrorIO('GloCol:: Could not find match for '//Val2Char(ib)//' or '//Val2Char(jb)//' in U Assembly')
!         print *, 'Problem -- match not found for ', ib, jb !gm kluge
 9200 UMAT(M,1)=UMAT(M,1)+DTRANE(IE,JE)+BTRANE(IE,JE)                    
!                                                                        
 9300 CONTINUE                                                           
 9999 CONTINUE                                                           
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
!                                                                        
!                                                                        
      RETURN                                                             
      END SUBROUTINE GLOCOL                                                                
