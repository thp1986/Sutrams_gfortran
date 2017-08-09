!     SUBROUTINE        I  N  D  A  T  0       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO INPUT ,OUTPUT, AND ORGANIZE A PORTION OF                      
! ***  UNIT-fINP INPUT DATA (DATASETS 5 THROUGH 7)                        
!                                                                       
      SUBROUTINE INDAT0 () 
      USE ITERAT 
      USE PARAMS 
      USE CONTRL 
      USE SOLVI 
      USE ITSOLR 
      USE ITSOLI 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE GRAVEC
      USE MODSOR
      USE SOLVC
      USE SOLVN
      USE ATSDATA
      use SutraMSPrecision 
      USE MSErrorHandler
      implicit none
      LOGICAL ::lrc 
      CHARACTER (len = 128) aline 
      CHARACTER(10) CDUM10 
      CHARACTER(14) UTYPE (2) 
      CHARACTER(6) STYPE (2) 
      CHARACTER(10) CSOLVP, CSOLVU (NSPE)
      integer (I4B) :: &
        NLSKIP 
      integer (I4B) :: &
        K, &
        KK, &
        M
      DATA UTYPE (1) / ' TEMPERATURES ' / , UTYPE (2) /                 &
      'CONCENTRATIONS' /                                                
      DATA STYPE (1) / 'ENERGY' / , STYPE (2) / 'SOLUTE' / 
      SAVE UTYPE, STYPE 
!                                                                       
!.....INPUT DATASET 5: NUMERICAL CONTROL PARAMETERS                     
      MSErrorValue%cDataSet='  5'
      CALL SKPCOM (fINP, NLSKIP) 
      READ (fINP,*,iostat=ios) UP, GNUP, (GNUU (K), K = 1, NSPE) 
      if(ios/=0) call ErrorIO('Error specifying numerical control parameters')
      WRITE (fLST, 60) UP, GNUP 
      DO K = 1, NSPE 
      WRITE (fLST, 70) GNUU (K), TRIM(ADJUSTL(SPNAME(K))) 
      ENDDO 
   60 FORMAT(////11X,'N U M E R I C A L   C O N T R O L   D A T A',//,  &
     &   11X,F15.5,5X,'"UPSTREAM WEIGHTING" FACTOR'/                    &
     &   11X,1PD15.4,5X,'SPECIFIED PRESSURE',                           &
     &   ' BOUNDARY CONDITION FACTOR')                                  
   70 FORMAT(11X,1PD15.4,5X,'SPECIFIED ',A,                             &
     &   ' BOUNDARY CONDITION FACTOR')                                  
!                                                                       
!.....INPUT DATASET 6: TEMPORAL CONTROL AND SOLUTION CYCLING DATA       
      MSErrorValue%cDataSet='  6'
      CALL SKPCOM (fINP, NLSKIP) 
      READ (fINP,*,iostat=ios) &
        ITMAX, DELT, TMAX, ITCYC, DTMULT, DTMAX, NPCYC, NUCYC                                                             
      if(ios/=0) call ErrorIO('Error specifying temporal control and solution cycling data')
      !Set dtmin to DELT (initial time step)
        dtmin=DELT
                                                                        
      WRITE (fLST, 120) &
        ITMAX, DELT, TMAX, ITCYC, DTMULT, DTMAX, NPCYC, NUCYC                                                             

  120 FORMAT(1H1////11X,'T E M P O R A L   C O N T R O L   A N D   ',   &
     &   'S O L U T I O N   C Y C L I N G   D A T A',                   &
     &   //11X,I15,5X,'MAXIMUM ALLOWED NUMBER OF TIME STEPS'            &
     &   /11X,1PD15.4,5X,'INITIAL TIME STEP (IN SECONDS)'               &
     &   /11X,1PD15.4,5X,'MAXIMUM ALLOWED SIMULATION TIME (IN SECONDS)' &
     &   //11X,I15,5X,'TIME STEP MULTIPLIER CYCLE (IN TIME STEPS)'      &
     &   /11X,0PF15.5,5X,'MULTIPLICATION FACTOR FOR TIME STEP CHANGE'   &
     &   /11X,1PD15.4,5X,'MAXIMUM ALLOWED TIME STEP (IN SECONDS)'       &
     &   //11X,I15,5X,'FLOW SOLUTION CYCLE (IN TIME STEPS)'             &
     &   /11X,I15,5X,'TRANSPORT SOLUTION CYCLE (IN TIME STEPS)')        
      IF (NPCYC.GE.1.AND.NUCYC.GE.1) GOTO 140 

  130 call ErrorIO('NPCYC and NUCYC must be >= 1')
  140 IF (NPCYC.EQ.1.OR.NUCYC.EQ.1) GOTO 160 
      call ErrorIO('NPCYC or NUCYC must = 1')
  160 IF (DELT.LE.DTMAX) GOTO 180 
      call ErrorIO('DELT must be <= DTMAX')
  180 CONTINUE 
!.....SET MAXIMUM ALLOWED TIME STEPS IN SIMULATION FOR                  
!        STEADY-STATE FLOW AND STEADY-STATE TRANSPORT SOLUTION MODES    
      IF (ISSFLO.EQ.1) THEN 
         NPCYC = ITMAX + 1 
         NUCYC = 1 
      ENDIF 
      IF (ISSTRA.EQ.1) ITMAX = 1 
!                                                                       
!.....INPUT DATASET 7A: OUTER (PICARD) ITERATION CONTROLS               
      MSErrorValue%cDataSet=' 7A'
      CALL SKPCOM (fINP, NLSKIP)
      READ (fINP, *,iostat=ios) ITRMAX
      if(ios/=0) call ErrorIO('Error specifying ITRMAX')
      if (ITRMAX<1) ITRMAX=1
      if (ITRMAX.GT.1) then
        BACKSPACE (fINP)
        READ (fINP, *,iostat=ios) ITRMAX, RPMAX, (RUMAX (K), K = 1, NSPE) 
      end if
      if(ios/=0) call ErrorIO('Error specifying outer(picard) iteration controls')
      IF (ITRMAX - 1) 192, 192, 194 
  192 WRITE (fLST, 193) 
  193 FORMAT(////11X,'I T E R A T I O N   C O N T R O L   D A T A',     &
     &   //11X,'  NON-ITERATIVE SOLUTION')                              
      GOTO 196 
  194 WRITE (fLST, 195) ITRMAX, RPMAX, (RUMAX (K), TRIM(ADJUSTL(SPNAME(K))), K = 1, NSPE)                                              
  195 FORMAT(////11X,'I T E R A T I O N   C O N T R O L   D A T A',     &
     &   //11X,I15,5X,'MAXIMUM NUMBER OF ITERATIONS PER TIME STEP',     &
     &   /11X,1PD15.4,5X,'ABSOLUTE CONVERGENCE CRITERION FOR FLOW',     &
     &   ' SOLUTION'50(/11X,1PD15.4,5X,'ABSOLUTE CONVERGENCE CRITERION',&
     &   ' FOR ',A,' TRANSPORT SOLUTION'))                              
  196 CONTINUE 
!                                                                       
!.....INPUT DATASETS 7B AND 7C:  INNER (SOLVER) ITERATION PARAMETERS    
      MSErrorValue%cDataSet=' 7B'
      CALL SKPCOM (fINP, NLSKIP) 
      READ (fINP,*,iostat=ios) CSOLVP 
      if(ios/=0) call ErrorIO('Error specifying parameters for P')
      kk = - 1 
      DO M = 0, NSLVRS - 1 
        IF (CSOLVP.EQ.SOLWRD (M) ) kk = M 
      ENDDO 
      KSOLVP = kk 
      IF (KK>0.AND.KK<4) THEN 
         BACKSPACE (fINP) 
         READ (fINP,*,iostat=ios) CSOLVP, ITRMXP, ITOLP, TOLP, NSAVEP
         if(ios/=0) call ErrorIO('Error specifying parameters for P')
         IF (KSOLVP.EQ.2) THEN
            ITOLP = 0
         END IF
         NSAVEP = 10
      ELSEIF (kk>3) then 
         call ErrorIO('Hook for additional solver - should not execute')
      ENDIF 
!                                                                       
!.....INITIALIZE CSOLVU                                                 
      CSOLVU = '          ' 
!                                                                       
!.....READ SLAP PARAMETERS FOR EACH SPECIES                             
      MSErrorValue%cDataSet=' 7C'
      CALL SKPCOM (fINP, NLSKIP) 
      IF (NSPE.LT.2) THEN 
         READ (fINP,*,iostat=ios) CSOLVU (1) 
         if(ios/=0) call ErrorIO('Error specifying parameters for U')
         kk = - 1 
         DO M = 0, NSLVRS - 1 
           IF (CSOLVU (1) .EQ.SOLWRD (M) ) kk = M 
         ENDDO 
         KSOLVU (1) = kk 
         IF (kk>0.and.kk<4) THEN 
            BACKSPACE (fINP) 
            READ (fINP,*,iostat=ios) CSOLVU(1), ITRMXU(1), ITOLU(1), TOLU(1), NSAVEU(1)
            if(ios/=0) call ErrorIO('Error specifying parameters for U')
            IF (KSOLVU(1).EQ.2) THEN
               ITOLU(1) = 0
            END IF
         ELSEIF (kk>3) then 
	        call ErrorIO('Hook for additional solver - should not execute')
         ENDIF
      ELSE 
  270    READ (fINP, *,iostat=ios) K 
         if(ios/=0) call ErrorIO('Error specifying parameters for U')

!.......CONTINUE READING DATASET 7C IF K EXCEEDS NSPE                   
         IF (K.GT.NSPE) GOTO 270 
!.......TERMINATE DATASET 7C WHEN K EQUALS ZERO                         
         IF (K.EQ.0) GOTO 280 
!.......BACKSPACE AND CONTINUE READING                                  
         BACKSPACE (fINP) 
         READ (fINP, * ) K, CSOLVU (K) 
         kk = - 1 
         DO M = 0, NSLVRS - 1 
         IF (CSOLVU (K) .EQ.SOLWRD (M) ) kk = M 
         enddo 
         KSOLVU (k) = kk 
         IF (kk==0) THEN 
            BACKSPACE (fINP) 
            READ (fINP,'(a128)',iostat=ios) aline 
            if(ios/=0) call ErrorIO('Error specifying parameter for U')
            READ (aline,*,iostat=ios) K, CSOLVU(K), StartUTime(K)                                                
            if(ios/=0) then
              READ (aline,*,iostat=ios) K, CSOLVU(K)
              StartUTime(K)=-1e6
            end if
            if(ios/=0) call ErrorIO('Error specifying parameter for U')
         ELSE IF (kk>0.and.kk<4) THEN 
            BACKSPACE (fINP) 
            READ (fINP,'(a128)',iostat=ios) aline 
            if(ios/=0) call ErrorIO('Error specifying parameter for U')
            READ (aline,*,iostat=ios) K, CSOLVU(K), ITRMXU(K), ITOLU(K), TOLU(K), NSAVEU(K), StartUTime(K)                                                
            if(ios/=0) then
              READ (aline,*,iostat=ios) K, CSOLVU(K), ITRMXU(K), ITOLU(K), TOLU(K), NSAVEU(K)
              StartUTime(K)=-1e6
            end if
            if(ios/=0) call ErrorIO('Error specifying parameter for U')

            IF (KSOLVU(K).EQ.2) THEN
               ITOLU(K) = 0
            END IF
         ELSEIF (kk>3) then 
	        call ErrorIO('Hook for additional solver - should not execute')
         ENDIF 
!.......CONTINUE READING DATASET 7C UNTIL K EQUALS ZERO                 
         GOTO 270 
      ENDIF 
  280 CONTINUE 
!                                                                       
!.....ERROR CHECKING FOR SPECIES SOLVER PARAMETERS                      
      DO K = 1, NSPE 
        IF (LEN_TRIM (CSOLVU(K) ) .LT.1) THEN 
!           INSTOP = INSTOP - 1 
           WRITE (fLST, 290) K, trim(adjustl(SPNAME(K)))
           call ErrorIO(('Error specifying U solver for species'//trim(adjustl(Val2Char(k)))))
        ENDIF 
      ENDDO 
  290 FORMAT(//11X,'SOLVER PARAMETERS NOT DEFINED FOR SPECIES',         &
     &             I3,2X,'[',A,']')                                   

 1000 RETURN 
      END SUBROUTINE INDAT0                         
