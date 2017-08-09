!     SUBROUTINE        S  R  C  1  S  P  E    SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO READ AND ORGANIZE ENERGY OR SOLUTE MASS SOURCE DATA.          
! ***  FOR SINGLE SPECIES DATASETS                                      
                                                                        
      SUBROUTINE SRC1SPE (QUIN, IQSOU, IQSOUT) 

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
        IQSOU (MNSOU, NSPE) 
      
      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        IQSOUT, IQCU, IQU, &
        NSOUI, NIQU, &
        NLSKIP
      REAL (DP) :: &
        QUINC 
!
!
      MSErrorValue%cDataSet='SRC'
!                                                                       
!.....NSOUI IS ACTUAL NUMBER OF SOLUTE MASS OR ENERGY SOURCE NODES      
      NSOUI = MNSOU - 1 
      IQSOUT = 1 
!                                                                       
  950 FORMAT(////11X,'THE SPECIFIED TIME VARIATIONS ARE ',              &
     &   'USER-PROGRAMMED IN SUBROUTINE  B C T I M E .')                
!                                                                       
!.....TEST IF ENERGY TRANSPORT                                          
      IF (ME.EQ.1) GOTO 1150 
!.....SOLUTE TRANSPORT                                                  
 1050 WRITE (fLST, 1060) 1, trim(adjustl(SPNAME(1)))
 1060 FORMAT(////////11X,'S O L U T E   S O U R C E   D A T A'/         &
     &   19X,'S P E C I E S [',I3,']'/24X,A,//                          &
     &   11X,'**** NODES AT WHICH SOURCES OR SINKS OF SOLUTE ',         &
     &   'MASS ARE SPECIFIED ****'//11X,'NODE NUMBER',10X,              &
     &   'SOLUTE SOURCE(+)/SINK(-)'/11X,'(MINUS INDICATES',5X,          &
     &   '(SOLUTE MASS/SECOND)'/12X,'TIME-VARYING'/12X,                 &
     &   'SOURCE OR SINK)'//)                                           
      GOTO 1300 
!.....ENERGY TRANSPORT                                                  
 1150 WRITE (fLST, 1160) 1, trim(adjustl(SPNAME(1)))
 1160 FORMAT(////////11X,'E N E R G Y   S O U R C E   D A T A',/        &
     &   19X,'S P E C I E S [',I3,']'/24X,A,//                          &
     &   11X,'**** NODES AT WHICH SOURCES OR SINKS OF ',                &
     &   'ENERGY ARE SPECIFIED ****'//11X,'NODE NUMBER',10X,            &
     &   'ENERGY SOURCE(+)/SINK(-)'/11X,'(MINUS INDICATES',5X,          &
     &   '(ENERGY/SECOND)'/12X,'TIME-VARYING'/12X,                      &
     &   'SOURCE OR SINK)'//)                                           
!                                                                       
!                                                                       
!.....INPUT DATASET 18                                                  
 1300 CALL SKPCOM (fINP, NLSKIP) 
!.....RESET SOURCE NODE COUNTER FOR EACH SPECIES                        
      NIQU = 0 
 1305 CONTINUE 
      READ (fINP, * ) IQCU 
      IF (IQCU.GT.0) THEN 
         BACKSPACE (fINP) 
         READ (fINP, * ) IQCU, QUINC 
      ELSEIF (IQCU.EQ.0) THEN 
         GOTO 1700 
      ENDIF 
      NIQU = NIQU + 1 
      IQSOU (NIQU, 1) = IQCU 
      IF (IQCU.LT.0) IQSOUT = - 1 
      IQU = IABS (IQCU) 
      QUIN (IQU, 1) = QUINC 
!.....TEST FOR CONSTANT U SOURCE/SINK WITH TIME                         
      IF (IQCU.GT.0) GOTO 1450 
!.....TIME VARYING U SOURCE/SINK                                        
      WRITE (fLST, 1500) IQCU 
      GOTO 1600 
!.....CONSTANT U SOURCE/SINK                                            
 1450 WRITE (fLST, 1500) IQCU, QUINC 
 1500 FORMAT(11X,I10,13X,1PE14.7) 
 1600 GOTO 1305 
!.....TEST IF THE NUMBER OF SOURCE NODES EQUALS THE SPECIFIED           
!.....NUMBER                                                            
 1700 IF (NIQU.EQ.NSOU (1) ) GOTO 1890 
!.....END SIMULATION IF THERE NEED BE CORRECTIONS TO DATASET 18         
!.....TEST IF ENERGY TRANSPORT SPECIES                                  
      IF (ME.EQ.1) GOTO 1760 
!.....ERROR IN SOLUTE TRANSPORT SOURCES/SINKS                           
 1740 WRITE (fLST, 1750) NIQU, NSOU (1), 1, trim(adjustl(SPNAME(1)))
 1750 FORMAT(////11X,'THE NUMBER OF SOLUTE SOURCE NODES READ, ',I9,/    &
     &   11X,'IS NOT EQUAL TO THE NUMBER SPECIFIED,   ',I9,/            &
     &   11X,'S P E C I E S [',I3,']'/16X,A,////                        &
     &   11X,'PLEASE CORRECT DATA AND RERUN'////////                    &
     &   22X,'S I M U L A T I O N   H A L T E D'/                       &
     &   22X,'_________________________________')                       
!
      call ErrorIO('SRC1SPE: Error the number of solute sources, &
NIQU ['//trim(adjustl(Val2Char(NIQU)))//'], does not equal the number specified, NSOU ['//trim(adjustl(Val2Char(NSOU(1))))//']')
!
!.....ERROR IN ENERGY TRANSPORT SOURCES/SINKS                           
 1760 WRITE (fLST, 1770) NIQU, NSOU (1), 1, trim(adjustl(SPNAME(1)))
 1770 FORMAT(////11X,'THE NUMBER OF ENERGY SOURCE NODES READ, ',I9,/    &
     &   11X,'IS NOT EQUAL TO THE NUMBER SPECIFIED,   ',I9,/            &
     &   11X,'S P E C I E S [',I3,']'/16X,A,////                        &
     &   11X,'PLEASE CORRECT DATA AND RERUN'////////                    &
     &   22X,'S I M U L A T I O N   H A L T E D'/                       &
     &   22X,'_________________________________')                       
!
      call ErrorIO('SRC1SPE: Error the number of energy sources, &
NIQU ['//trim(adjustl(Val2Char(NIQU)))//'], does not equal the number specified, NSOU ['//trim(adjustl(Val2Char(NSOU(1))))//']')
!
 1890 IF (IQSOUT.EQ. - 1) WRITE (fLST, 950) 
!.....RETURN TO CALLING ROUTINE                                         
      RETURN 
      END SUBROUTINE SRC1SPE                        
