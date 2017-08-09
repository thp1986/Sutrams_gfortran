!     SUBROUTINE        B  N  D  M  S  P  E    SUTRA-MS VERSION 2004.1  
!                                                                       
! *** PURPOSE :                                                         
! ***  TO READ AND ORGANIZE SPECIFIED TEMPERATURE OR CONCENTRATION      
!      DATA FOR SINGLE SPECIES DATASETS.                                
!                                                                       
      SUBROUTINE BNDMSPE (IUBC, UBC, IUBCT) 
      USE CONTRL 
      USE PARAMS 
      USE FUNITS 
      USE DIMS
      USE SutraStorage, ONLY : MultiSpeciesBC
      use SutraMSPrecision
      USE MSErrorHandler
      implicit none
      integer (I4B) :: &
        IUBC (NBCN, NSPE)
      real (DP) :: &
        UBC (NBCN, NSPE) 
      integer (I4B) :: &
        IU (NSPE), IPU (NSPE) 
      !locals
      integer (I4B) :: &
        I, K, &
        IUT, &
        IUBCT, ISTOPU, &
        NLSKIP
        
!
!
      MSErrorValue%cDataSet='BND'
!                                                                       
!.....INITIALIZE FLAGS                                                  
      IUBCT = 1 
      ISTOPU = 0 
!                                                                       
!.....RESET COUNTERS FOR THE NUMBER OF SPECIFIED NODES FOR              
!     SINGLE SPECIES                                                    
      IU = 0 
      IPU = NPBC 
!                                                                       
!.....INPUT DATASET 20                                                  
 1000 CALL SKPCOM (fINP, NLSKIP) 
 1125 READ (fINP, * ) K 
!.....CONTINUE READING DATASET 20 IF K EXCEEDS THE NUMBER OF            
!     SIMULATED SPECIES AND K EXCEEDS ZERO                              
      IF (K.GT.NSPE) GOTO 1125 
!.....TERMINATE READING DATASET 20 IF K EQUALS ZERO                     
      IF (K.EQ.0) GOTO 1180 
!.....INCREMENT IPU AND IU COUNTERS                                     
      IPU (K) = IPU (K) + 1 
      IU (K) = IU (K) + 1 
!                                                                       
!.....BACKSPACE AND REREAD DATASET 20 DATA                              
      BACKSPACE (fINP) 
!.....READ NODE NUMBER OF SPECIFIED CONCENTRATION OR TEMPERATURE        
!     NODE                                                              
      IUT = IPU (K) 
      READ (fINP, * ) K, IUBC (IUT, K) 
      IF (IUBC (IUT, K) .GT.0) THEN 
         BACKSPACE (fINP) 
!         READ (fINP, * ) K, IUBC (IUT, K), UBC (IUT, K) 
         READ (fINP, * ) K, MultiSpeciesBC(K)%SpecifiedU(IU(K))%node, MultiSpeciesBC(K)%SpecifiedU(IU(K))%U
      ELSEIF (IUBC (IUT, K) .LT.0) THEN 
         IUBCT = - 1 
      ENDIF 
!.....CONTINUE READING DATASET 20 UNTIL K EQUALS ZERO                   
      GOTO 1125 
!                                                                       
 1180 CONTINUE 
!                                                                       
!.....DETERMINE IF NUBC DIMENSIONS HAVE BEEN EXCEEDED                   
      DO K = 1, NSPE 
      IF (IU (K) .EQ.NUBC (K) ) CYCLE 
      ISTOPU = 1 
!.......WRITE ERROR MESSAGE IF IU DOES NOT EQUAL NUBC(K)                
      WRITE (fLST, 2010) K, trim(adjustl(SPNAME(K))), IU (K), NUBC (K) 
      ENDDO 
!                                                                       
 2010 FORMAT(////11X,'S P E C I E S [',I3,']',/11X,                     &
     & 'ACTUAL NUMBER OF SPECIFIED'1XA,1X,'NODES READ, ',               &
     &   I9,','/11X,'IS NOT EQUAL TO NUMBER SPECIFIED IN INPUT, ',I9)   
!                                                                       
!.....END SIMULATION IF THERE NEED BE CORRECTIONS TO DATASET 20         
      IF (ISTOPU.EQ.0) GOTO 6000 
      IF (ME) 2500, 3000, 3500 
!.....SOLUTE TRANSPORT                                                  
 2500 IF (ISTOPU.EQ.1) WRITE (fLST, 2550) SUM (IU), SUM (NUBC) 
 2550 FORMAT(////11X,'ACTUAL NUMBER OF SPECIFIED CONCENTRATION NODES',  &
     &   ' READ, ',I9,','/11X,'IS NOT EQUAL TO NUMBER SPECIFIED IN',    &
     &   ' INPUT, ',I9)                                                 
      GOTO 4800 
!.....SOLUTE AND ENERGY TRANSPORT                                       
 3000 IF (ISTOPU.EQ.1) WRITE (fLST, 3050) SUM (IU), SUM (NUBC) 
 3050 FORMAT(////11X,'ACTUAL NUMBER OF SPECIFIED CONCENTRATION AND',    &
     &   ' ENERGY NODES READ, ',I9,','/11X,'IS NOT EQUAL TO NUMBER',    &
     &   ' SPECIFIED IN INPUT, ',I9)                                    
      GOTO 4800 
!.....ENERGY TRANSPORT                                                  
 3500 IF (ISTOPU.EQ.1) WRITE (fLST, 3550) SUM (IU), SUM (NUBC) 
 3550 FORMAT(////11X,'ACTUAL NUMBER OF SPECIFIED TEMPERATURE NODES',    &
     &   ' READ, ',I9,', IS NOT EQUAL TO NUMBER SPECIFIED IN',          &
     &   ' INPUT, ',I9)                                                 
 4800 WRITE (fLST, 5000) 
 5000 FORMAT(////11X,'PLEASE CORRECT DATA AND RERUN.'////////           &
     &   22X,'S I M U L A T I O N   H A L T E D'/                       &
     &   22X,'_________________________________')                       
!
      call ErrorIO('BND1SPE: Actual specified concentration/temperature nodes does not equal NUBC &
 ['//trim(adjustl(Val2Char(SUM(NUBC))))//']')
!                                                                       
!.....WRITE DATA SET 20 DATA TO OUTPUT FILE                             
 6000 DO K = 1, NSPE 
      IF (NUBC (K) .EQ.0) CYCLE 
!.......TEST IF ENERGY TRANSPORT SPECIES                                
      IF (K.EQ.NESP) GOTO 6020 
!.......SOLUTE TRANSPORT                                                
 6010 WRITE (fLST, 6015) K, trim(adjustl(SPNAME(K)))
 6015 FORMAT  (////11X,'**** NODES AT WHICH SOLUTE CONCENTRATIONS ARE ',&
     &   'SPECIFIED TO BE INDEPENDENT OF LOCAL FLOWS AND FLUID SOURCES',&
     &   ' ****'/19X,'S P E C I E S [',I3,']'/24X,A,                    &
     &   //12X,'NODE',13X,'CONCENTRATION'//)                            
      GOTO 6030 
 6020 WRITE (fLST, 6025) K, trim(adjustl(SPNAME(K)))
 6025 FORMAT  (////11X,'**** NODES AT WHICH TEMPERATURES ARE ',         &
     &   'SPECIFIED TO BE INDEPENDENT OF LOCAL FLOWS AND FLUID SOURCES',&
     &   ' ****'/19X,'S P E C I E S [',I3,']'/24X,A,                    &
     &   //12X,'NODE',15X,'TEMPERATURE'//)                              
!                                                                       
!.......WRITE SPECIFIED CONCENTRATIONS OR TEMPERATURES TO OUTPUT        
!       FILE (UNIT fLST)                                                  
 6030 IUT = NPBC 
      DO I = 1, NUBC (K) 
        IUT = IUT + 1 
        IF (MultiSpeciesBC(K)%SpecifiedU(I)%node .GT.0) THEN 
           WRITE (fLST, 6100) MultiSpeciesBC(K)%SpecifiedU(I)%node, MultiSpeciesBC(K)%SpecifiedU(I)%U
        ELSEIF (MultiSpeciesBC(K)%SpecifiedU(I)%node .LT.0) THEN 
           WRITE (fLST, 6100) MultiSpeciesBC(K)%SpecifiedU(I)%node
        ENDIF 
      ENDDO 
!                                                                       
!.......TEST IF TIME-DEPENDENT SPECIFIED U NODES ARE PRESENT            
      IF (IUBCT.NE. - 1) CYCLE 
!                                                                       
!.......TEST IF THIS IS A ENERGY TRANSPORT SPECIES                      
      IF (K.EQ.NESP) GOTO 6050 
 6040 WRITE (fLST, 6045) 
 6045 FORMAT  (//12X,'TIME-DEPENDENT SPECIFIED CONCENTRATION'/12X,'IS ',&
     &     'INDICATED BY NEGATIVE NODE NUMBER')                         
      CYCLE 
 6050 WRITE (fLST, 6055) 
 6055 FORMAT  (//11X,'TIME-DEPENDENT SPECIFIED TEMPERATURE'/12X,'IS ',  &
     &     'INDICATED BY NEGATIVE NODE NUMBER')                         
                                                                        
      ENDDO 
!                                                                       
 6100 FORMAT(11X,I9,6X,1PD20.13) 
!                                                                       
!.....RETURN TO CALLING ROUTINE                                         
      RETURN 
      END SUBROUTINE BNDMSPE                        
