!     SUBROUTINE        B  N  D  1  S  P  E    SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO READ AND ORGANIZE SPECIFIED TEMPERATURE OR CONCENTRATION      
!      DATA FOR SINGLE SPECIES DATASETS.                                
!                                                                       
      SUBROUTINE BND1SPE (IUBC, UBC, IUBCT) 
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
      !locals
      integer (I4B) :: &
        IU, IPU, &
        IUBCT, ISTOPU, &
        NLSKIP
      integer (I4B) :: &
        IVAL
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
!.....TEST IF ENERGY TRANSPORT SPECIES                                  
      IF (NESP.EQ.1) GOTO 550 
!.....SOLUTE TRANSPORT                                                  
  500 WRITE (fLST, 510) 1, trim(adjustl(SPNAME(1)))
  510 FORMAT(////11X,'**** NODES AT WHICH SOLUTE CONCENTRATIONS ARE ',  &
     &   'SPECIFIED TO BE INDEPENDENT OF LOCAL FLOWS AND FLUID SOURCES',&
     &   ' ****'/19X,'S P E C I E S [',I3,']'/24X,A,                    &
     &   //12X,'NODE',13X,'CONCENTRATION'//)                            
      GOTO 1000 
  550 WRITE (fLST, 560) 1, trim(adjustl(SPNAME(1)))
  560 FORMAT(////11X,'**** NODES AT WHICH TEMPERATURES ARE ',           &
     &   'SPECIFIED TO BE INDEPENDENT OF LOCAL FLOWS AND FLUID SOURCES',&
     &   ' ****'/19X,'S P E C I E S [',I3,']'/24X,A,                    &
     &   //12X,'NODE',15X,'TEMPERATURE'//)                              
!                                                                       
!.....INPUT DATASET 20                                                  
 1000 CALL SKPCOM (fINP, NLSKIP) 
 1125 IPU = IPU + 1 
      IU = IU + 1 
!      READ (fINP, * ) MultiSpeciesBC(1)%SpecifiedU(IU)%node
      READ (fINP, * ) IVAL
      IF ( IVAL.GT.0 ) THEN 
         BACKSPACE (fINP) 
         READ (fINP, * ) MultiSpeciesBC(1)%SpecifiedU(IU)%node, MultiSpeciesBC(1)%SpecifiedU(IU)%U
         WRITE (fLST, 1150) MultiSpeciesBC(1)%SpecifiedU(IU)%node, MultiSpeciesBC(1)%SpecifiedU(IU)%U
      ELSEIF ( IVAL.LT.0 ) THEN 
         IUBCT = - 1 
         WRITE (fLST, 1150) MultiSpeciesBC(1)%SpecifiedU(IU)%node 
      ELSE 
         GOTO 1180 
      ENDIF 
 1150 FORMAT(11X,I9,6X,1PD20.13) 
      GOTO 1125 
 1180 IPU = IPU - 1 
      IU = IU - 1 
      IF (IU.EQ.NUBC (1) ) GOTO 1200 
      ISTOPU = 1 
!                                                                       
!.....WRITE ERROR MESSAGE IF IU DOES NOT EQUAL NUBC(1)                  
      WRITE (fLST, 1190) 1, trim(adjustl(SPNAME(1))), IU, NUBC (1) 
 1190 FORMAT(////11X,'S P E C I E S [',I3,']'/11X,                      &
     & 'ACTUAL NUMBER OF SPECIFIED'1XA,1X,'NODES READ, ',               &
     &   I9,','/11X,'IS NOT EQUAL TO NUMBER SPECIFIED IN INPUT, ',I9)   
!.....TEST IF TIME-DEPENDENT SPECIFIED U NODES ARE PRESENT              
 1200 IF (IUBCT.NE. - 1) GOTO 2000 
!                                                                       
!.....TEST IF THIS IS A ENERGY TRANSPORT SPECIES                        
      IF (1.EQ.NESP) GOTO 1215 
 1205 WRITE (fLST, 1206) 
 1206 FORMAT(//12X,'TIME-DEPENDENT SPECIFIED CONCENTRATION'/12X,'IS ',  &
     &   'INDICATED BY NEGATIVE NODE NUMBER')                           
      GOTO 2000 
 1215 WRITE (fLST, 1216) 
 1216 FORMAT(//11X,'TIME-DEPENDENT SPECIFIED TEMPERATURE'/12X,'IS ',    &
     &   'INDICATED BY NEGATIVE NODE NUMBER')                           
!                                                                       
!.....END OF 1 TO NSPE SPECIFIED CONCENTRATION/TEMPERATURE NODES        
 2000 CONTINUE 
!                                                                       
!.....END SIMULATION IF THERE NEED BE CORRECTIONS TO DATASET 20         
      IF (ISTOPU.EQ.0) GOTO 6000 
      IF (ME) 3500, 4000, 4500 
!.....SOLUTE TRANSPORT                                                  
 3500 IF (ISTOPU.EQ.1) WRITE (fLST, 3550) IU, SUM (NUBC) 
 3550 FORMAT(////11X,'ACTUAL NUMBER OF SPECIFIED CONCENTRATION NODES',  &
     &   ' READ, ',I9,','/11X,'IS NOT EQUAL TO NUMBER SPECIFIED IN',    &
     &   ' INPUT, ',I9)                                                 
      GOTO 4800 
!.....SOLUTE AND ENERGY TRANSPORT                                       
 4000 IF (ISTOPU.EQ.1) WRITE (fLST, 4050) IU, SUM (NUBC) 
 4050 FORMAT(////11X,'ACTUAL NUMBER OF SPECIFIED CONCENTRATION AND',    &
     &   ' ENERGY NODES READ, ',I9,','/11X,'IS NOT EQUAL TO NUMBER',    &
     &   ' SPECIFIED IN INPUT, ',I9)                                    
      GOTO 4800 
!.....ENERGY TRANSPORT                                                  
 4500 IF (ISTOPU.EQ.1) WRITE (fLST, 4550) IU, SUM (NUBC) 
 4550 FORMAT(////11X,'ACTUAL NUMBER OF SPECIFIED TEMPERATURE NODES',    &
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
!.....RETURN TO CALLING ROUTINE                                         
 6000 RETURN 
      END SUBROUTINE BND1SPE                        
