!     SUBROUTINE        S  O  U  R  C  E       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO READ AND ORGANIZE FLUID MASS SOURCE DATA AND ENERGY OR        
! ***  SOLUTE MASS SOURCE DATA.                                         
!                                                                       
      SUBROUTINE SOURCE (QIN, UIN, IQSOP, QUIN, IQSOU, IQSOPT, IQSOUT) 
      USE ITERAT 
      USE PARAMS 
      USE CONTRL 
      USE FUNITS 
      USE DIMS
      use SutraMSPrecision
      USE MSErrorHandler
      USE MSErrorHandler

      IMPLICIT NONE

      REAL (DP) :: &
        QIN (NN), UIN (NN, NSPE)
      INTEGER (I4B) :: &
        IQSOP (NSOP)
      REAL (DP) :: &
        QUIN (NN, NSPE)
      INTEGER (I4B) :: &
        IQSOU (MNSOU, NSPE)                                               
      REAL (DP) :: &
        UINC (NSPE) 
      !locals
      INTEGER (I4B) :: &
        K, KK, &
        NSOPI, NSOUI, NIQP, &
        NLSKIP, &
        IQSOPT, IQSOUT, IQCP, IQP
      REAL (DP) :: &
        QINC
!
!
      MSErrorValue%cDataSet='SRC'
!                                                                       
!.....NSOPI IS ACTUAL NUMBER OF FLUID SOURCE NODES                      
!.....NSOUI IS ACTUAL NUMBER OF SOLUTE MASS OR ENERGY SOURCE NODES      
      NSOPI = NSOP - 1 
      NSOUI = MNSOU - 1 
      IQSOPT = 1 
      IQSOUT = 1 
      NIQP = 0 
      IF (NSOPI.EQ.0) GOTO 1000
      MSErrorValue%cDataSet='SRC'
      select case (ME)
        !solute transport
        case(-1)
          goto 50
        !solute and heat transport
        case(0)
          goto 150
        !heat transport
        case(1)
          goto 200
        case default
          call ErrorIO('Error: ME ['//trim(Val2Char(ME))//'] must be -1, 0, or 1')
      end select
!     SOLUTE TRANSPORT                                                  
   50 WRITE (fLST, 100) 
  100 FORMAT(1H1////11X,'F L U I D   S O U R C E   D A T A'             &
     &   ////11X,'**** NODES AT WHICH FLUID INFLOWS OR OUTFLOWS ARE ',  &
     &   'SPECIFIED ****'//11X,'NODE NUMBER',10X,                       &
     &   'FLUID INFLOW(+)/OUTFLOW(-)',5X,'SOLUTE CONCENTRATION OF'      &
     &   /11X,'(MINUS INDICATES',5X,'(FLUID MASS/SECOND)',              &
     &   12X,'INFLOWING FLUID'/12X,'TIME-VARYING',39X,                  &
     &   '(MASS SOLUTE/MASS WATER)'/12X,'FLOW RATE OR'/12X,             &
     &   'CONCENTRATION)'//)                                            
      GOTO 300 
      !SOLUTE AND HEAT TRANSPORT                                        
  150 WRITE (fLST, 160) trim(adjustl(SPNAME(1)))
  160 FORMAT(1H1////11X,'F L U I D   S O U R C E   D A T A'             &
     &   ////11X,'**** NODES AT WHICH FLUID INFLOWS OR OUTFLOWS ARE ',  &
     &   'SPECIFIED ****'//11X,'NODE NUMBER',10X,                       &
     &   'FLUID INFLOW(+)/OUTFLOW(-)',5X,A,                             &
     &   /11X,'(MINUS INDICATES',5X,'(FLUID MASS/SECOND)',12X,          &
     &   'OF INFLOWING FLUID'/12X,'TIME-VARYING'/12X,'FLOW OR'/12X,     &
     &   'TEMPERATURE)'//)                                              
      GOTO 300 
      !HEAT TRANSPORT                                                   
  200 WRITE (fLST, 210) 
  210 FORMAT(1H1////11X,'F L U I D   S O U R C E   D A T A'             &
     &   ////11X,'**** NODES AT WHICH FLUID INFLOWS OR OUTFLOWS ARE ',  &
     &   'SPECIFIED ****'//11X,'NODE NUMBER',10X,                       &
     &   'FLUID INFLOW(+)/OUTFLOW(-)',5X,'TEMPERATURE {DEGREES CELSIUS}'&
     &   /11X,'(MINUS INDICATES',5X,'(FLUID MASS/SECOND)',12X,          &
     &   'OF INFLOWING FLUID'/12X,'TIME-VARYING'/12X,'FLOW OR'/12X,     &
     &   'TEMPERATURE)'//)                                              
!                                                                       
!.....INPUT DATASET 17                                                  
  300 CONTINUE 
      MSErrorValue%cDataSet=' 17'
      CALL SKPCOM (fINP, NLSKIP) 
      READ (fINP, * ) IQCP 
      IF (IQCP.GT.0) THEN 
         BACKSPACE (fINP) 
         READ (fINP, * ) IQCP, QINC 
!........FLUID SOURCE WITH SPECIFIED CONCENTRATIONS                     
         IF (QINC.GT.0D0) THEN 
            BACKSPACE (fINP) 
            READ (fINP, * ) IQCP, QINC, (UINC (K), K = 1, NSPE) 
         ENDIF 
      ELSEIF (IQCP.EQ.0) THEN 
         GOTO 700 
      ENDIF 
      NIQP = NIQP + 1 
      IQSOP (NIQP) = IQCP 
      IF (IQCP.LT.0) IQSOPT = - 1 
      IQP = IABS (IQCP) 
      QIN (IQP) = QINC 
!.....SPECIFIED CONCENTRATIONS FOR POSITIVE SOURCE [INJECTION]          
      IF (QINC.GT.0D0) THEN 
         DO K = 1, NSPE 
         UIN (IQP, K) = UINC (K) 
         ENDDO 
      ENDIF 
!.....TEST FOR CONSTANT SOURCE FUNCTION WITH TIME AT NODE               
      IF (IQCP.GT.0) GOTO 450 
      WRITE (fLST, 500) IQCP 
      GOTO 600 
!.....TEST IF POSITIVE SOURCE [INJECTION]                               
  450 IF (QINC.GT.0) GOTO 460 
!.....WRITE NODE, AND RATE FOR FLUID SINK                               
      WRITE (fLST, 500) IQCP, QINC 
      GOTO 600 
!.....WRITE NODE,RATE,AND SPECIFIED CONCENTRATION OF FIRST SPECIES      
  460 WRITE (fLST, 500) IQCP, QINC, UINC (1) 
  500 FORMAT(11X,I10,13X,1PE14.7,16X,1PE14.7) 
!.....CONTINUE READING FLUID SOURCES UNTIL IQCP = 0                     
  600 GOTO 300 
!.....TEST IF THE NUMBER OF FLUID SOURCES IS EQUAL TO DIMENSIONS        
  700 IF (NIQP.EQ.NSOPI) GOTO 890 
!.....END SIMULATION IF THERE NEED BE CORRECTIONS TO DATASET 17         
      WRITE (fLST, 750) NIQP, NSOPI 
  750 FORMAT(////11X,'THE NUMBER OF FLUID SOURCE NODES READ, ',I9,      &
     &   ' IS NOT EQUAL TO THE NUMBER SPECIFIED, ',I9////               &
     &   11X,'PLEASE CORRECT DATA AND RERUN'////////                    &
     &   22X,'S I M U L A T I O N   H A L T E D'/                       &
     &   22X,'_________________________________')                       
!
!
      call ErrorIO('SOURCE: Error the number of fluid sources, &
NIQP ['//trim(adjustl(Val2Char(NIQP)))//'], does not equal the number specified, NSOPI ['//trim(adjustl(Val2Char(NSOPI)))//']')
!
!.....WRITE REMAINING SPECIES TO UNIT fLST                                
  890 IF (NSPE.GT.1) WRITE (fLST, 900) NSPE, (trim(adjustl(SPNAME(K))), K = 2, NSPE) 
  900 FORMAT(/11X,'A D D I T I O N A L  ',                              &
     &            'F L U I D   S O U R C E   D A T A',                  &
     &       /11X,'SPECIES  2 TO',I3,                                   &
     &       /11X,'      NODE',13XA10,20XA10,                           &
     &  100(:/34XA10,20XA10))                                            
!.....LOAD SPECIFIED CONCENTRATIONS FOR NSPE 2->NSPE INTO UINC          
!.....IN ORDER TO WRITE TO UNIT fLST                                      
      DO IQP = 1, NSOPI 
      IQCP = IQSOP (IQP) 
      IF (QIN (IQCP) .GT.0D0) THEN 
         UINC = 0D0 
         DO K = 2, NSPE 
         UINC (K) = UIN (IQCP, K) 
         ENDDO 
         DO K = 2, NSPE, 2 
         IF (K.EQ.2) THEN 
            IF (NSPE.EQ.2) THEN 
               WRITE (fLST, 910) IQCP, UINC (K) 
            ELSE 
               WRITE (fLST, 910) IQCP, (UINC (KK), KK = K, K + 1) 
            ENDIF 
         ENDIF 
         IF (K.GT.2.AND. (K + 1) .LT. (NSPE+1) ) WRITE (fLST, 915)        &
           (UINC (KK), KK = K, K + 1)                                     
         IF (K.GT.2.AND. (K + 1) .GT. (NSPE) ) WRITE (fLST, 915) UINC (K) 
         ENDDO 
      ENDIF 
      ENDDO 
  910 FORMAT(11XI10,13X1PE14.7,16X1PE14.7,100(:/34X1PE14.7,16X1PE14.7)) 
  915 FORMAT(34X1PE14.7,16X1PE14.7,100(:/34X1PE14.7,16X1PE14.7)) 
!                                                                       
      IF (IQSOPT.EQ. - 1) WRITE (fLST, 950) 
  950 FORMAT(////11X,'THE SPECIFIED TIME VARIATIONS ARE ',              &
     &   'USER-PROGRAMMED IN SUBROUTINE  B C T I M E .')                
!                                                                       
!.....SPECIFIED CONCENTRATION/HEAT SOURCE DATA FOR EACH SPECIES         
 1000 IF (NSOUI.EQ.0) GOTO 9000 
!.....READ SPECIFIED CONCENTRATION/HEAT SOURCE DATA FOR EACH SPECIES    
      IF (NSPE.LT.2) THEN 
         CALL SRC1SPE (QUIN, IQSOU, IQSOUT) 
      ELSE 
         CALL SRCMSPE (QUIN, IQSOU, IQSOUT) 
      ENDIF 
!.....RETURN TO CALLING ROUTINE                                         
 9000 RETURN 
!                                                                       
      END SUBROUTINE SOURCE                         
