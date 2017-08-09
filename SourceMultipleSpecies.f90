
!     SUBROUTINE        S  R  C  M  S  P  E    SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO READ AND ORGANIZE ENERGY OR SOLUTE MASS SOURCE DATA.          
! ***  FOR MULTIPLE SPECIES DATASETS                                    
                                                                        
      SUBROUTINE SRCMSPE (QUIN, IQSOU, IQSOUT) 

      USE ITERAT 
      USE PARAMS 
      USE CONTRL 
      USE FUNITS 
      USE DIMS
      USE SutraMSPrecision
      USE MSErrorHandler

      IMPLICIT NONE

      REAL (DP) :: &
        QUIN (NN, NSPE)
      INTEGER (I4B) :: &
        IQSOU (MNSOU, NSPE), & 
        NIQU (NSPE)

      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        I, K, & 
        NSOUI, IQSOUT, IQCU, IQU, &
        NLSKIP
      REAL (DP) :: &
        QUINC

!
      MSErrorValue%cDataSet='SRC'
!                                                                       
!.....NSOUI IS ACTUAL NUMBER OF SOLUTE MASS OR ENERGY SOURCE NODES      
      NSOUI = MNSOU - 1 
      IQSOUT = 1 
!.....INITIALICE SOURCE NODE COUNTER FOR EACH SPECIES                   
      NIQU = 0 
!                                                                       
!.......INPUT DATASET 18 - MULTI-SPECIES                                
      CALL SKPCOM (fINP, NLSKIP) 
 1000 READ (fINP, * ) K 
!.....IF K EXCEEDS NSPE CONTINUE READING INPUT DATASET 18               
      IF (K.GT.NSPE) GOTO 1000 
!.....TERMINATE READING DATASET 18 IF K EQUALS ZERO                     
      IF (K.EQ.0) GOTO 1150 
!.....BACKSPACE FILE AND READ ALL DATA ON THE LINE                      
      BACKSPACE (fINP) 
      READ (fINP, * ) K, IQCU, QUINC 
      NIQU (K) = NIQU (K) + 1 
!.....TEST IF SPECIFIED DIMENSIONS HAVE BEEN EXCEEDED FOR               
!     THIS SPECIES                                                      
      IF (NIQU (K) .GT.NSOU (K) ) GOTO 1150 
!.....VALID ENTRY - PUSH TO APPROPRIATE POSITION                        
!     IN IQSOU AND QUIN ARRAYS                                          
      IQSOU (NIQU (K), K) = IQCU 
      IF (IQCU.LT.0) IQSOUT = - 1 
      IQU = IABS (IQCU) 
      QUIN (IQU, K) = QUINC 
      GOTO 1000 
!                                                                       
!.....OUTSIDE DATA READ LOOP                                            
 1150 CONTINUE 
!.....TEST IF THE NUMBER OF SOURCE NODES EQUALS THE SPECIFIED           
!.....NUMBER OF SOURCE NODES (NSOU(K))                                  
      DO K = 1, NSPE 
        IF (NIQU (K) .EQ.NSOU (K) ) CYCLE 
!.......END SIMULATION IF THERE NEED BE CORRECTIONS TO DATASET 18       
!       TEST IF ENERGY TRANSPORT SPECIES                                
        IF (K.EQ.NESP) GOTO 2010 
!.......ERROR IN SOLUTE TRANSPORT SOURCES/SINKS                         
 2000   WRITE (fLST, 2005) NIQU (K), NSOU (K), K, trim(adjustl(SPNAME(K)))
 2005   FORMAT  (////11X,'THE NUMBER OF SOLUTE SOURCE NODES READ, ',I9,/  &
     &       11X,'IS NOT EQUAL TO THE NUMBER SPECIFIED,   ',I9,/          &
     &       11X,'S P E C I E S [',I3,']'/16X,A,////                      &
     &       11X,'PLEASE CORRECT DATA AND RERUN'////////                  &
     &       22X,'S I M U L A T I O N   H A L T E D'/                     &
     &       22X,'_________________________________')                     
       call ErrorIO('SRCSPE: Error in solute transport sources/sinks')
!.......ERROR IN ENERGY TRANSPORT SOURCES/SINKS                         
 2010   WRITE (fLST, 2015) NIQU (K), NSOU (K), K, trim(adjustl(SPNAME(K)))
 2015   FORMAT  (////11X,'THE NUMBER OF ENERGY SOURCE NODES READ, ',I9,/  &
     &       11X,'IS NOT EQUAL TO THE NUMBER SPECIFIED,   ',I9,/          &
     &       11X,'S P E C I E S [',I3,']'/16X,A,////                      &
     &       11X,'PLEASE CORRECT DATA AND RERUN'////////                  &
     &       22X,'S I M U L A T I O N   H A L T E D'/                     &
     &       22X,'_________________________________')                     
        call ErrorIO('SRCSPE: Error in energy transport sources/sinks')
      ENDDO 
!.....WRITE MESSAGE TO OUTPUT FILE (UNIT fLST) IF USER PROGRAMMED         
!     BOUNDARY CONDITIONS HAVE BEEN SPECIFIED                           
      IF (IQSOUT.EQ. - 1) WRITE (fLST, 3000) 
!                                                                       
!.....OUTPUT DATASET 18 DATA TO OUTPUT FILE (UNIT fLST)                   
      DO K = 1, NSPE 
        IF (NSOU (K) .EQ.0) CYCLE 
!.......TEST IF ON ENERGY TRANSPORT SPECIES                             
        IF (K.NE.NESP) THEN 
!.........SOLUTE TRANSPORT                                              
         WRITE (fLST, 3010) K, trim(adjustl(SPNAME(K)))
        ELSE 
!.........ENERGY TRANSPORT                                              
         WRITE (fLST, 3020) K, trim(adjustl(SPNAME(K)))
        ENDIF 
        DO I = 1, NSOU (K) 
          IQCU = IQSOU (I, K) 
          IQU = ABS (IQCU) 
          QUINC = QUIN (IQU, K) 
!.........WRITE CONSTANT U SOURCE/SINK                                  
          IF (IQCU.GT.0) WRITE (fLST, 3030) IQCU, QUINC 
        ENDDO 
      ENDDO 
!                                                                       
!.....FORMAT STATEMENTS                                                 
 3000 FORMAT(////11X,'THE SPECIFIED TIME VARIATIONS ARE ',              &
     &   'USER-PROGRAMMED IN SUBROUTINE  B C T I M E .')                
 3010 FORMAT(////////11X,'S O L U T E   S O U R C E   D A T A'/         &
     &   19X,'S P E C I E S [',I3,']'/24X,A,//                          &
     &   11X,'**** NODES AT WHICH SOURCES OR SINKS OF SOLUTE ',         &
     &   'MASS ARE SPECIFIED ****'//11X,'NODE NUMBER',10X,              &
     &   'SOLUTE SOURCE(+)/SINK(-)'/11X,'(MINUS INDICATES',5X,          &
     &   '(SOLUTE MASS/SECOND)'/12X,'TIME-VARYING'/12X,                 &
     &   'SOURCE OR SINK)'//)                                           
 3020 FORMAT(////////11X,'E N E R G Y   S O U R C E   D A T A',/        &
     &   19X,'S P E C I E S [',I3,']'/24X,A,//                          &
     &   11X,'**** NODES AT WHICH SOURCES OR SINKS OF ',                &
     &   'ENERGY ARE SPECIFIED ****'//11X,'NODE NUMBER',10X,            &
     &   'ENERGY SOURCE(+)/SINK(-)'/11X,'(MINUS INDICATES',5X,          &
     &   '(ENERGY/SECOND)'/12X,'TIME-VARYING'/12X,                      &
     &   'SOURCE OR SINK)'//)                                           
 3030 FORMAT(11X,I10,13X,1PE14.7) 
!                                                                       
!.....RETURN TO CALLING ROUTINE                                         
      RETURN 
      END SUBROUTINE SRCMSPE                        
