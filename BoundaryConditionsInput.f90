!     SUBROUTINE        B  O  U  N  D          SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO READ AND ORGANIZE SPECIFIED PRESSURE DATA AND                 
! ***  SPECIFIED TEMPERATURE OR CONCENTRATION DATA.                     
!                                                                       
      SUBROUTINE BOUND (IPBC, PBC, IUBC, UBC, IPBCT, IUBCT) 
      USE CONTRL 
      USE PARAMS 
      USE FUNITS 
      USE DIMS
      USE SutraStorage, ONLY : SpecifiedPBC
      USE SutraMSPrecision
      IMPLICIT NONE
      integer (I4B) :: &
        IPBC (NBCN)
      real (DP) :: &
        PBC (NBCN)
      integer (I4B) :: &
        IUBC (NBCN, NSPE)
      real (DP) :: &
        UBC (NBCN, NSPE)
      !LOCALS                                                             
      integer (I4B) :: inode
      integer (I4B) :: &
        IPBCT, IUBCT, ISTOPP, ISTOPU, IP, IPU, &
        I, K, KK, &
        NLSKIP
!                                                                       
!                                                                       
      IPBCT = 1 
      IUBCT = 1 
      ISTOPP = 0 
      ISTOPU = 0 
      IP = 0 
      IPU = 0 
      WRITE (fLST, 50) 
   50 FORMAT(1H1////11X,'B O U N D A R Y   C O N D I T I O N S') 
      IF (NPBC.EQ.0) GOTO 400 
      WRITE (fLST, 100) 
  100 FORMAT(//11X,'**** NODES AT WHICH PRESSURES ARE',                 &
     &   ' SPECIFIED ****'/)                                            
      IF (ME) 107, 110, 114 
!.....SOLUTE TRANSPORT                                                  
  107 WRITE (fLST, 108) 1, trim(adjustl(SPNAME(1)))
  108 FORMAT(11X,'     (AS WELL AS SOLUTE CONCENTRATION OF ANY'         &
     &   /16X,' FLUID INFLOW WHICH MAY OCCUR AT THE POINT'              &
     &   /16X,' OF SPECIFIED PRESSURE)'//55X,'SPECIES ['I3,']'/57X,A,   &
     &   /12X,'NODE',18X,'PRESSURE',                                    &
     &   13X,'CONCENTRATION'//)                                         
      GOTO 120 
!.....SOLUTE AND ENERGY TRANSPORT                                       
  110 WRITE (fLST, 111) 1, trim(adjustl(SPNAME(1)))
  111 FORMAT(11X,'     (AS WELL AS SOLUTE CONCENTRATION OF ANY'         &
     &   /16X,' FLUID INFLOW WHICH MAY OCCUR AT THE POINT'              &
     &   /16X,' OF SPECIFIED PRESSURE)'//55X,'SPECIES ['I3,']',         &
     &   /12X,'NODE',18X,'PRESSURE',                                    &
     &   15X,A//)                                                     
      GOTO 120 
!.....ENERGY TRANSPORT                                                  
  114 WRITE (fLST, 115) 
  115 FORMAT(11X,'     (AS WELL AS TEMPERATURE {DEGREES CELSIUS} OF ANY'&
     &   /16X,' FLUID INFLOW WHICH MAY OCCUR AT THE POINT'              &
     &   /16X,' OF SPECIFIED PRESSURE)'//12X,'NODE',18X,                &
     &   'PRESSURE',13X,'  TEMPERATURE'//)                              
!                                                                       
!.....INPUT DATASET 19                                                  
  120 CALL SKPCOM (fINP, NLSKIP) 
  125 IPU = IPU + 1 
      READ (fINP, * ) inode
      IF (inode > 0) THEN 
         BACKSPACE (fINP) 
         READ (fINP, * ) &
            SpecifiedPBC(IPU)%node, SpecifiedPBC(IPU)%P, ( SpecifiedPBC(IPU)%U(K), K = 1, NSPE )
         WRITE (fLST, 160) &
            SpecifiedPBC(IPU)%node, SpecifiedPBC(IPU)%P, SpecifiedPBC(IPU)%U(1) 
      ELSEIF (inode < 0) THEN 
         IPBCT = - 1
         !set SpecifiedPBC(IPU)%node for user programmed boundary conditions - v2009.1 bug fix
         SpecifiedPBC(IPU)%node = inode 
         WRITE (fLST, 160) SpecifiedPBC(IPU)%node 
      ELSE 
         GOTO 180 
      ENDIF 
  160 FORMAT(7X,I9,6X,1PD20.13,6X,1PD20.13) 
      GOTO 125 
!.....SET IPU TO ACTUAL NUMBER OF SPECIFIED PRESSURE NODES              
  180 IPU = IPU - 1 
!.....WRITE REMAINING SPECIES FOR SPECIFIED PRESSURE NODES              
!.....TO UNIT fLST                                                        
      IF (NSPE.GT.1) WRITE (fLST, 185) NSPE, (trim(adjustl(SPNAME(K))), K = 2, NSPE) 
  185 FORMAT(/11X,'A D D I T I O N A L  ',                              &
     &            'F L U I D   U  D A T A',                             &
     &       /11X,'FOR SPECIFIED PRESSURE NODES',                       &
     &       /11X,'SPECIES  2 TO',I3,                                   &
     &       /12X,'NODE',13XA10,16XA10,                                 &
     &  100(:/29XA10,16XA10))                                           
!.....WRITE SPECIFIED CONCENTRATIONS FOR NSPE 2->NSPE                   
!.....TO UNIT fLST                                                        
      DO I = 1, IPU 
        DO K = 2, NSPE, 2 
          IF (K.EQ.2) THEN 
             IF (NSPE.EQ.2) THEN 
                WRITE (fLST, 190) &
                  SpecifiedPBC(I)%node, SpecifiedPBC(I)%U(K)
             ELSE 
                WRITE (fLST, 190) &
                  SpecifiedPBC(I)%node, ( SpecifiedPBC(I)%U(KK), KK = K, K + 1 )
             ENDIF 
          ENDIF 
          IF (K.GT.2.AND. (K + 1) .LT. (NSPE+1) ) &
            WRITE (fLST, 195) ( SpecifiedPBC(I)%U(KK), KK = K, K + 1 )
          IF (K.GT.2.AND. (K + 1) .GT. (NSPE) ) &
            WRITE (fLST, 195) SpecifiedPBC(I)%U(K)
        ENDDO 
      ENDDO 
  190 FORMAT(7XI9,6X1PE20.13,6X1PE20.13) 
  195 FORMAT(22X1PE20.13,6X1PE20.13) 
!                                                                       
!.....SAVE NUMBER OF LAST SPECIFIED PRESSURE NODE                       
      IP = IPU 
      IF (IP.EQ.NPBC) GOTO 200 
      ISTOPP = 1 
!.....TEST IF SPECIFIED PRESSURES AND CONCENTRATIONS ARE                
!     TIME-DEPENDENT                                                    
  200 IF (IPBCT.NE. - 1) GOTO 400 
!.....TIME-DEPENDENT SPECIFIED PRESSURE NODES                           
!     WRITE MESSAGE APPROPRIATE FOR SIMULATION TO UNIT fLST               
      IF (ME) 205, 210, 215 
!.....SOLUTE TRANSPORT                                                  
  205 WRITE (fLST, 206) 
  206 FORMAT(//12X,'TIME-DEPENDENT SPECIFIED PRESSURE'/12X,'OR INFLOW ',&
     &   'CONCENTRATION INDICATED'/12X,'BY NEGATIVE NODE NUMBER')       
      GOTO 400 
!.....SOLUTE AND ENERGY TRANSPORT                                       
  210 WRITE (fLST, 211) 
  211 FORMAT(//12X,'TIME-DEPENDENT SPECIFIED PRESSURE'/12X,'OR INFLOW ',&
     &   'CONCENTRATION AND TEMPERATURE INDICATED'/12X,'BY NEGATIVE ',  &
     &   'NODE NUMBER')                                                 
      GOTO 400 
!.....ENERGY TRANSPORT                                                  
  215 WRITE (fLST, 216) 
  216 FORMAT(//11X,'TIME-DEPENDENT SPECIFIED PRESSURE'/12X,'OR INFLOW ',&
     &   'TEMPERATURE INDICATED'/12X,'BY NEGATIVE NODE NUMBER')         
!                                                                       
!.....TEST IF ANY SPECIFIED NODES FOR THIS SPECIES                      
  400 IF (SUM (NUBC) .EQ.0) GOTO 6000 
!                                                                       
!.....CALL APPROPRIATE DATASET 20 SUBROUTINE FOR SINGLE AND             
!     MULTIPLE SPECIES DATASETS                                         
      IF (NSPE.LT.2) THEN 
!.......SINGLE SPECIES                                                  
         CALL BND1SPE (IUBC, UBC, IUBCT) 
      ELSE 
!.......MULTIPLE SPECIES                                                
         CALL BNDMSPE (IUBC, UBC, IUBCT) 
      ENDIF 
!                                                                       
 6000 IF (IPBCT.EQ. - 1.OR.IUBCT.EQ. - 1) WRITE (fLST, 7000) 
 7000 FORMAT(////11X,'THE SPECIFIED TIME VARIATIONS ARE ',              &
     &   'USER-PROGRAMMED IN SUBROUTINE  B C T I M E .')                
!                                                                       
!                                                                       
!.....RETURN TO CALLING ROUTINE                                         
      RETURN 
      END SUBROUTINE BOUND                          
