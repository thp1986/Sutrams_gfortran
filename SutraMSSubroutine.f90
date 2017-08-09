
!     SUBROUTINE        S  U  T  R  A          SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  MAIN ROUTINE THAT ASSEMBLES AND SOLVES SUTRA SIMULATION.         
!                                                                       

      SUBROUTINE SUTRA (TITLE1, TITLE2) 

      USE ITERAT 
      USE PARAMS 
      USE CONTRL 
      USE SOLVI 
      USE JCOLS 
      USE MODSOR 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE OBS
      USE SOURCEITEMS
      !MS specific modules
      USE TBCS 
      USE ATSDATA
      USE UserSpecifiedOutputTime 
      USE SutraStorage
      USE MSErrorHandler
      USE SutraMSPrecision
      USE TECPLOT 

      IMPLICIT NONE

      LOGICAL :: &
        PRNALL, PRNLST, PRNNOD, PRNELE, PRNOBS 
!.....VERSION 1.1
      LOGICAL :: &
        PRNTPN, PRNTPE 
      LOGICAL :: &
        SOLVER
      LOGICAL :: lrc 
      LOGICAL :: lDELTReduced=.false.
      CHARACTER (len=80) :: TITLE1, TITLE2

        REAL (DP) :: &
            DELTP_, &
            DELTU_

      !LOCAL VARIABLES
        INTEGER (I4B) :: &
          ML, &
          MLSAVE, &
          MATDIM, &
          IHALFB, &
          ISTOP, &
          IGOI, &
          ISTOP_OUTER, &
          IGOI_OUTER

         INTEGER (I4B) :: &
          I, &
          K, &
          KKK, &
          IP, &
          NK, &
          KPU
         INTEGER (I4B) :: &
          itemp
         REAL (DP) :: &
          DELTSAVE, &
          TSECP0, &
          TSECU0, &
          CST, &
          CWT, &
          BDELP1, &
          BDELU1, &
          BDELP, &
          BDELU, &
          DELTM1, &
          RPTIME0, &
          RPTIME1, &
          RP, &
          do1
        REAL (DP) :: &
          ObsTimeDiff, &
          EPSILON=1E-10

!                                                                       
!.....INITIALIZE DNTIME AND DPTIME FOR TRANSIENT                        
!     BOUNDARY CONDITIONS AND OUTPUT PRINTING AT SPECIFIC TIMES         
      DNTIME = 0.0D0 
      DPTIME = 0.0D0 
!                                                                       
!.....INITIALIZE TIME STEP NUMBER                                       
      IT = 0 
!!  V2009.1 BUG FIX FOR BCTIME - 11/05/2009                                                                       
!!.....INITIALIZE COUNTERS FOR TRANSIENT B.C.s PROGRAMMED BY USER        
!      IQSOPT = 0 
!      IQSOUT = 0 
!      IPBCT = 0 
!      IUBCT = 0 
!                                                                       
!.....INITIALIZE VMAG, VANG1, AND VANG2                                 
      VMAG = 0.0D0 
      VANG1 = 0.0D0 
      VANG2 = 0.0D0 
!......MOVED TO SutraMSMainProgram - VERSION 1.1
!!                                                                       
!!.....INPUT INITIAL OR RESTART CONDITIONS AND INITIALIZE PARAMETERS     
!!        (READ UNIT-fICS DATA)                                            
!      write (*,*) 'Reading ICS data...'
!      CALL INDAT2 (PVEC, UVEC, PM1, UM1, UM2, CS1, CS2, CS3, SL, SR,    &
!                   RCIT, SW, DSWDP, PBC, IPBC, IPBCT)                          

!                                                                       
!.....GET FIRST PRINT TIME FROM OUTPUT CONTROL                          
      IF (fOTM.GT.0) THEN
         !initialize DPTIME to simulation starting time (TSEC)
         DPTIME=TSEC
         !Scan otm file until DPTIME exceeds simulation starting time (TSEC)
         write (*,*) 'Reading first user specified Output TiMe...'
         do while (DPTIME<TSEC+DELT) 
           if (.not. ReadUserSpecifiedOutputTime() ) call ErrorIO('SUTRA:: Error reading next print time')
         end do
      ENDIF 
                                                                        
!                                                                       
!.....INITIALIZE DELTSAVE TO DELT                                       
      DELTSAVE = DELT 
!                                                                       
!.....SET STARTING TIME OF SIMULATION CLOCK                             
!     TSEC=TSTART                                                       
      TSECP0 = TSEC 
      TSECU0 = TSEC 
      TMIN = TSEC / 60.D0 
      THOUR = TMIN / 60.D0 
      TDAY = THOUR / 24.D0 
      TWEEK = TDAY / 7.D0 
      TMONTH = TDAY / 30.4375D0 
      TYEAR = TDAY / 365.25D0 
!                                                                       
!.....IF TRANSIENT B.C.s ARE ALLOWED READ INITIAL TIME AT WHICH         
!     TO CHANGE B.C.s                                                   
      lRdTBCTime=.TRUE.
      IF (fTBC.GT.0) then
        write (*,*) 'Reading transient BC data...'
        if(.not. ReadTBCData()) call ErrorIO('SUTRA: Error reading transient BC data')
      end if
!                                                                       
!.....SAVE VALUES FOR ENERGY TRANSPORT                                  
      CST = CS 
      CWT = CW 
!                                                                       
!.....set imaxiter to 0 initially                                       
!      IF (fATS.GT.0) THEN 
      IF ( LATS ) THEN 
         imaxiter = 0 
         NMaxCount = 0 
      ENDIF 
!                                                                       
!.....OUTPUT INITIAL CONDITIONS OR STARTING CONDITIONS                  
      IF (ISSTRA.NE.1) THEN 
         write (*,*) 'Outputting initial/starting conditions to fLST...'
         IF (IABS (KTYPE) .EQ.3) THEN 
            CALL OUTPTCLST_3D(0, 0, 0, PVEC, UVEC, VMAG, VANG1, VANG2, SW) 
         ELSE 
            CALL OUTPTCLST_2D(0, 0, 0, PVEC, UVEC, VMAG, VANG1, SW) 
         ENDIF 
      ENDIF 
!                                                                       
      write (*,*) 'Allocation/Initialization complete...Beginning Time Steps'
!
!.....SET SWITCHES AND PARAMETERS FOR SOLUTION WITH STEADY-STATE FLOW   
      IF (ISSFLO.NE.1) GOTO 1000 
      ML = 1 
      NOUMAT = 0 
      ISSFLO = 2 
      ITER = 0 
      DLTPM1 = DELTP 
      DLTUM1 = DELTU 
      BDELP1 = 1D0 
      BDELP = 0.0D0 
      BDELU = 0.0D0 
      GOTO 1100 
!                                                                       
!                                                                       
! **********************************************************************
!.....BEGIN TIME STEP **************************************************
! **********************************************************************
 1000 IT = IT + 1 
      ITER = 0 
      ML = 0 
      NOUMAT = 0 
      lDELTReduced=.false.
!
!.....SET NOUMAT TO OBTAIN U SOLUTION BY SIMPLE BACK SUBSTITUTION       
!        BEGINNING ON SECOND TIME STEP AFTER A PRESSURE SOLUTION        
!        IF THE SOLUTION IS NON-ITERATIVE (ITRMAX=1)
!        BACKSUBSTITUTION ONLY ALLOWED FOR CASES WHERE NSPE.EQ.1 - SUTRA-MS 1.1                    
      IF ( MOD(IT - 1, NPCYC).NE.0 .AND. &
           MOD(IT, NPCYC).NE.0 .AND.     &
           IT.GT.2 .AND.                 &
           ITRMAX.EQ.1 .AND.             &
           NSPE.EQ.1 ) THEN
         NOUMAT = 1                     
      END IF
!.....CHOOSE SOLUTION VARIABLE ON THIS TIME STEP:                       
!        ML=0 FOR P AND U, ML=1 FOR P ONLY, AND ML=2 FOR U ONLY.        
      IF ( IT.EQ.1.AND.ISSFLO.NE.2 ) GOTO 1005 
      IF ( MOD(IT, NPCYC).NE.0 ) ML = 2 
      IF ( MOD(IT, NUCYC).NE.0 ) ML = 1 
!.....MULTIPLY TIME STEP SIZE BY DTMULT EACH ITCYC TIME STEPS           
      DELTM1 = DELT 
      IF (MOD (IT, ITCYC) .EQ.0.AND.IT.GT.1) THEN 
         DELT = DELT * DTMULT 
         !IF USING Automatic TimeStep Module AND ITERATIONS DURING
         !LAST TIMESTEP EXCEEDED 1/2 IMAXTARGET THEN DO NOT INCREASE
         !THE TIMESTEP DURING THIS ITCYC TIME STEP (WAIT TILL NEXT ITCYC TIME STEP)
!         IF(fATS > 0) THEN
         IF( LATS ) THEN
           IF ( .not.LATSAllowTimeStepIncrease ) DELT = DELT / DTMULT
         END IF 
!........SET TIME STEP SIZE TO MAXIMUM ALLOWED SIZE, DTMAX              
         IF (DELT.GT.DTMAX) DELT = DTMAX 
      ENDIF 
!.....RESET DELT IF IT EXCEEDS TMAX                                     
      IF ( (TSEC + DELT) .GT.TMAX) DELT = TMAX - TSEC 
!                                                                       
!.....RESET DELT IF IT EXCEEDS NEXT OBS TIME
      if(allocated(TIMOBS)) then
        IF ( (TSEC + DELT) .GT.TIMOBS(ONOBS)) THEN
          if(.not. lDELTReduced) DELTSAVE = DELT 
          DELT = TIMOBS(ONOBS)-TSEC
          lTIMOBSRESET=.TRUE.
          lDELTReduced=.true.
        ENDIF 
      end if
!                                                                       
!.....RESET DELT IF IT EXCEEDS DNTIME                                   
      IF (fTBC.GT.0) THEN
         IF(lGetTBC) THEN
           DELT=dtmin
           lGetTBC=.FALSE.
         ELSE
           IF ( (TSEC + DELT) .GT.DNTIME) DELT = DNTIME-TSEC
           IF ( (TSEC + DELT) .EQ.DNTIME) lGetTBC=.TRUE.
        END IF
      ENDIF 
!                                                                       
!.....READ SPECIFIC OUTPUT PRINT TIMES, IF REQUIRED.                    
!     DPTIME ONLY CONTROLS OUTPUT TO fLST, fNOD, fELE, fTPN, AND fTPE OUTPUT FILES   
      IF (fOTM.GT.0) THEN 
!.......RESET DELT IF IT EXCEEDS DPTIME                                 
         IF ( (TSEC + DELT) .GT.DPTIME) THEN 
            if(.not. lDELTReduced) DELTSAVE = DELT 
            DELT = DPTIME-TSEC 
            lDELTReduced=.true.
         ENDIF 
      ENDIF 
                                                                        
!.....NO SIMPLE BACK SUBSTITUTION FOR U IF TIME STEP HAS CHANGED        
      IF (DELT.NE.DELTM1) NOUMAT = 0 
!.....INCREMENT SIMULATION CLOCK, TSEC, TO END OF NEW TIME STEP         
 1005 TSEC = TSEC + DELT 
      TMIN = TSEC / 60.D0 
      THOUR = TMIN / 60.D0 
      TDAY = THOUR / 24.D0 
      TWEEK = TDAY / 7.D0 
      TMONTH = TDAY / 30.4375D0 
      TYEAR = TDAY / 365.25D0 
!
!.....Convergence flags for outer PICARD iterations
      ISTOP_OUTER = 0
      IGOI_OUTER  = 0
!                                                                       
!.....SET TIME STEP FOR P AND/OR U, WHICHEVER ARE SOLVED FOR            
!        ON THIS TIME STEP                                              
      SELECT CASE (ML-1)
        CASE(-1)
          GO TO 1010
        CASE(0)
          GO TO 1020
        CASE(1)
          GO TO 1030
        CASE DEFAULT
          MSErrorValue%cDataSet='SUT'
          call ErrorIO('SUTRA solution: Internal error (ML-1 must be -1, 0, 1) ML-1 = ['//trim(adjustl(Val2Char(ML-1)))//']')
      END SELECT
      
 1010 DLTUM1 = DELTU 
      DLTPM1 = DELTP 
      GOTO 1040 
 1020 DLTPM1 = DELTP 
      GOTO 1040 
 1030 DLTUM1 = DELTU 
 1040 CONTINUE 
      DELTP = TSEC - TSECP0 
      DELTU = TSEC - TSECU0 
!.....SET PROJECTION FACTORS USED ON FIRST ITERATION TO EXTRAPOLATE     
!        AHEAD ONE-HALF TIME STEP                                       
      BDELP = (DELTP / DLTPM1) * 0.50D0 
      BDELU = (DELTU / DLTUM1) * 0.50D0 
      BDELP1 = BDELP + 1.0D0 
      BDELU1 = BDELU + 1.0D0 
!.....INCREMENT CLOCK FOR WHICHEVER OF P AND U WILL BE SOLVED FOR       
!        ON THIS TIME STEP                                              
      SELECT CASE (ML-1)
        CASE(-1)
          GO TO 1060
        CASE(0)
          GO TO 1070
        CASE(1)
          GO TO 1080
        CASE DEFAULT
          MSErrorValue%cDataSet='SUT'
          call ErrorIO('SUTRA solution: Internal error (ML-1 must be -1, 0, 1) ML-1 = ['//trim(adjustl(Val2Char(ML-1)))//']')
      END SELECT
 1060 TSECP0 = TSEC 
      TSECU0 = TSEC 
      GOTO 1090 
 1070 TSECP0 = TSEC 
      GOTO 1090 
 1080 TSECU0 = TSEC 
 1090 CONTINUE 
!                                                                       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....BEGIN ITERATION - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 1100 ITER = ITER + 1 
      if (fOTM>0) then
        WRITE (  *, 1102) IT, ITMAX, ITER, &
                          trim(adjustl(cTime(TSEC))), &
                          trim(adjustl(cTime(DELT))), &
                          trim(adjustl(cTime(DPTIME)))
      else
        WRITE (  *, 1103) IT, ITMAX, ITER, &
                          trim(adjustl(cTime(TSEC))), &
                          trim(adjustl(cTime(DELT)))
      end if
      WRITE (fSMY, 1104) IT, ITMAX, ITER, &
                       trim(adjustl(cTime(TSEC))), &
                       trim(adjustl(cTime(DELT)))
 1102 FORMAT (/1X, 57('-'), /                                           &
     &         1X, 'TIME STEP ', I10, ' OF ', I10, ':  ITERATION ', I10,&
     &        /1X, 'TOTAL TIME: ', A,                                   &
     &        /1X, 'DELT      : ', A,                                   &
     &        /1X, 'print time: ', A,                                   &
     &        /1X, 57('-'))                                          
 1103 FORMAT (/1X, 57('-'), /                                           &
     &         1X, 'TIME STEP ', I10, ' OF ', I10, ':  ITERATION ', I10,&
     &        /1X, 'TOTAL TIME: ', A,                                   &
     &        /1X, 'DELT      : ', A,                                   &
     &        /1X, 57('-'))                                          
 1104 FORMAT (                                                          &
     &         1X, 'TIME STEP ', I10, ' OF ', I10, ':  ITERATION ', I10,&
     &         1X, 'TOTAL TIME: ', A,                                   &
     &         1X, 'DELT      : ', A)                        
!.....VERSION 1.1
!.....DETERMINE IF TECPLOT NODAL OUTPUT WILL BE GENERATED
!     AT THIS TIME
      IF (fTPN>0) THEN
        IF (.NOT. PrintTpVelocity() ) &
          call ErrorIO('SUTRA: Could not determine if Tecplot Nodal Output time')
      END IF
!                                                                       
      SELECT CASE (ML-1)
        CASE(-1)
          GO TO 2000
        CASE(0)
          GO TO 2200
        CASE(1)
          GO TO 2400
        CASE DEFAULT
          MSErrorValue%cDataSet='SUT'
          call ErrorIO('SUTRA solution: Internal error (ML-1 must be -1, 0, 1) ML-1 = ['//trim(adjustl(Val2Char(ML-1)))//']')
      END SELECT

!.....SHIFT AND SET VECTORS FOR TIME STEP WITH BOTH P AND U SOLUTIONS   
      !forall and whole array operations                                
! 2000 forall (i = 1:nn) 
!        PITER (i) = PVEC (i) 
!        PVEL (i) = PVEC (i) 
!        RCITM1 (i) = RCIT (i) 
 2000 DO I = 1, NN
        PITER(I)  = PVEC(I) 
        PVEL(I)   = PVEC(I) 
        RCITM1(I) = RCIT(I) 
!.....RESET RCIT TO BASE DENSITY                                      
        RCIT(I)   = RHOW0 
      END DO 
!      forall (i = 1:nn, k = 1:nspe) 
!          UITER (i, k) = UVEC (i, k) 
!      endforall 
      DO I = 1, NN
        DO K = 1, NSPE
          UITER(I, K) = UVEC(I, K) 
        END DO
      END DO
!.....CALCULATE DENSITY BASED ON UPDATED CONCENTRATIONS                 
!      DO K = 1, NSPE 
!          RCIT (1:NN) = RCIT (1:NN) + DRWDU (K) * (UITER (1:NN, K) - URHOW0 (K) )                                                             
!      ENDDO 
      DO I = 1, NN
        DO K = 1, NSPE 
          RCIT(I) = RCIT(I) + DRWDU(K) * ( UITER(I, K) - URHOW0(K) )                                                             
        END DO 
      END DO
!                                                                       
!.....SPECIFIED PRESSURE FLUXES                                         
!      FORALL (IP = 1:NPBC) 
!        QPLITR (IP) = GNUP * (SpecifiedPBC(IP)%P - PITER ( IABS ( SpecifiedPBC(IP)%node ) ) ) 
!      ENDFORALL
      DO IP = 1, NPBC
        QPLITR(IP) = GNUP * ( SpecifiedPBC(IP)%P - PITER( IABS( SpecifiedPBC(IP)%node ) ) ) 
      END DO
!                                                                       
      IF (ITER.GT.1) GOTO 2600 
!      FORALL (I = 1:NN) 
!        PITER (I) = BDELP1 * PVEC (I) - BDELP * PM1 (I) 
!        PM1 (I) = PVEC (I) 
!      ENDFORALL 
      DO I = 1, NN
        PITER(I) = BDELP1 * PVEC(I) - BDELP * PM1(I) 
        PM1(I)   = PVEC(I) 
      END DO 

!.....SPECIES LOOP TO CALCULATE UITER USING BDELU1 AND BDELU            
!      FORALL (I = 1:NN, K = 1:NSPE) 
!          UITER (I, K) = BDELU1 * UVEC (I, K) - BDELU * UM1 (I, K) 
!          UM2 (I, K) = UM1 (I, K) 
!          UM1 (I, K) = UVEC (I, K) 
!      ENDFORALL 
      DO I = 1, NN
        DO K = 1, NSPE
          UITER(I, K) = BDELU1 * UVEC(I, K) - BDELU * UM1(I, K) 
          UM2(I, K)   = UM1(I, K) 
          UM1(I, K)   = UVEC(I, K) 
        END DO
      END DO 
!.....END SPECIES LOOP                        
      GOTO 2600 
!                                                                       
!.....SHIFT AND SET VECTORS FOR TIME STEP WITH P SOLUTION ONLY          
 2200 DO I = 1, NN
        PVEL(I)  = PVEC(I) 
        PITER(I) = PVEC(I) 
      END DO 
!                                                                       
      IF (ITER.GT.1) GOTO 2600 
!                                                                       
      DO I = 1, NN
        PITER(I) = BDELP1 * PVEC(I) - BDELP * PM1(I) 
        PM1(I)   = PVEC(I) 
!.......SET RCITM1 TO DENSITY FROM LAST ITERATION AND RCIT TO           
!       BASE FLUID DENSITY                                              
        RCITM1(I) = RCIT(I) 
!.......SET RCIT TO THE BASE FLUID DENSITY                              
        RCIT(I) = RHOW0 
      END DO
!.....SPECIES LOOP - reset UVEC to UITER                                
!      FORALL (I = 1:NN, K = 1:NSPE) 
!          UITER (I, K) = UVEC (I, K) 
!      ENDFORALL 
      DO I = 1, NN
        DO K = 1, NSPE 
          UITER(I, K) = UVEC(I, K) 
        END DO
      END DO 
!.....CALCULATE DENSITY BASED ON UPDATED CONCENTRATIONS                 
!      DO K = 1, NSPE 
!        RCIT (1:NN) = RCIT (1:NN) + DRWDU (K) * (UITER (1:NN, K) - URHOW0(K) )                                                             
!      ENDDO 
      DO I = 1, NN
        DO K = 1, NSPE 
          RCIT(I) = RCIT(I) + DRWDU(K) * ( UITER(I, K) - URHOW0(K) )                                                             
        END DO
      END DO 
!                                                                       
      GOTO 2600 
!                                                                       
!.....SHIFT AND SET VECTORS FOR TIME STEP WITH U SOLUTION ONLY          
 2400 IF ( ITER.EQ.1 ) THEN 
!         FORALL (I = 1:NN, K = 1:NSPE) 
!            UITER (I, K) = BDELU1 * UVEC (I, K) - BDELU * UM1 (I, K) 
!         ENDFORALL 
         DO I = 1, NN
           DO K = 1, NSPE
             UITER(I, K) = BDELU1 * UVEC(I, K) - BDELU * UM1(I, K) 
           END DO
         END DO
      ELSE 
!         FORALL (I = 1:NN, K = 1:NSPE) 
!           UITER (I, K) = UVEC (I, K) 
!         ENDFORALL 
         DO I = 1, NN
           DO K = 1, NSPE
             UITER(I, K) = UVEC(I, K) 
           END DO
         END DO 
      ENDIF 
!                                                                       
      IF (NOUMAT.EQ.1) GOTO 2480 
!                                                                       
      IF (ITER.GT.1) GOTO 2600 
!                                                                       
      DO I = 1, NN
        PITER(I)  = PVEC(I) 
        PVEL(I)   = PVEC(I) 
        RCITM1(I) = RCIT(I) 
      END DO
!                                                                       
!.....SPECIFIED PRESSURE FLUXES                                         
!      FORALL (IP = 1:NPBC) 
!        QPLITR (IP) = GNUP * (SpecifiedPBC(IP)%P - PITER ( IABS ( SpecifiedPBC(IP)%node ) ) ) 
!      ENDFORALL 
      DO IP = 1, NPBC
        QPLITR(IP) = GNUP * ( SpecifiedPBC(IP)%P - PITER( IABS( SpecifiedPBC(IP)%node ) ) ) 
      END DO 
!                                                                       
!.....UPDATE UM1 with UM2 and UVEC with UM1                             
! 2480 FORALL (I = 1:NN, K = 1:NSPE) 
!          UM2 (I, K) = UM1 (I, K) 
!          UM1 (I, K) = UVEC (I, K) 
!      ENDFORALL 
 2480 DO I = 1, NN
        DO K = 1, NSPE
          UM2(I, K) = UM1(I, K) 
          UM1(I, K) = UVEC(I, K) 
        END DO
      END DO
!                                                                       
 2600 CONTINUE 
!
!.....AUTOMATIC TIMESTEP
      IF ( LATS ) THEN                                                
        IF ( .not.LATSFailConvergence .and. .not. lContinue ) THEN 
!.........SAVE DATA AT THE BEGINNING OF TIMESTEP FOR USE                  
!.........IF THE IMAXTARGET IS EXCEEDED
           write (*,*) '....Updating ATS data for re-running timestep '//trim(adjustl(Val2Char(IT)))
           DELTP_  = DELTP 
           DELTU_  = DELTU 
           PVEC_   = PVEC 
           UVEC_   = UVEC 
           PM1_    = PM1 
           UM1_    = UM1 
           RCIT_   = RCIT 
           CS1_    = CS1 
           SW_     = SW 
           QPLITR_ = QPLITR
        ENDIF
      END IF
!                                                                       
!                                                                       
!.....INITIALIZE ARRAYS WITH VALUE OF ZERO                              
      MATDIM = NELT * NCBI 
!                                                                       
!     INITIALLY SET U_ML TO ML - WILL RESET
      U_ML = ML                            
!.....SPECIES LOOP FOR SOLUTION OF P AND U  
      DO 5800 K = 1, NSPE 
!.......GLOBAL VARIABLE IN PARAMS MODULE                                
         KSP = K 
!.......SET U_ML FOR SPECIES K
        IF ( NSPE>1 ) THEN
          IF(TSEC-StartUTime(K)<=0.AND.(U_ML(K)-1)==-1) THEN
            !P & U TO P
            U_ML(K)=0
          !U ONLY - CYCLE LOOP
          ELSE IF(TSEC-StartUTime(K)<=0.AND.(U_ML(K)-1)==1) THEN
            write(*,*) 'Note: skipping U('//trim(adjustl(Val2Char(k)))//') solution at time ['//trim(adjustl(Val2Char(TSEC)))//']'
            write(fSMY,'(1x,a)') 'Note: skipping U('//trim(adjustl(Val2Char(k)))//') &
solution at time ['//trim(adjustl(Val2Char(TSEC)))//']'
            CYCLE
          END IF
        END IF
!.......SET TEMPORARY VALUES OF CS AND CW                               
         IF ( KSP.EQ.NESP ) THEN 
            CS = CST 
            CW = CWT 
         ELSE 
            CS = 0.0D0 
            CW = 1.0D0 
         ENDIF 
!                                                                       
!                                                                       
!.....TEST IF P SOLUTION [(ML-1) = -1 OR 0]                             
         SELECT CASE (U_ML(K)-1)
          CASE(-1)
            GO TO 3000
          CASE(0)
            GO TO 3000
          CASE(1)
            GO TO 3300
          CASE DEFAULT
            MSErrorValue%cDataSet='SUT'
            call ErrorIO('SUTRA solution: Internal error testing if P solution &
(ML-1 must be -1, 0, 1) ML-1 = ['//trim(adjustl(Val2Char(ML-1)))//']')
         END SELECT

 3000    CALL ZERO (PMAT, MATDIM, 0.0D0) 
         CALL ZERO (PVEC, NN, 0.0D0) 
         CALL ZERO (VOL, NN, 0.0D0) 
!.....TEST IF U SOLUTION [(ML-1) = -1 OR 1]                             
         SELECT CASE (U_ML(K)-1)
          CASE(-1)
            GO TO 3300
          CASE(0)
            GO TO 3400
          CASE(1)
            GO TO 3300
          CASE DEFAULT
            MSErrorValue%cDataSet='SUT'
            call ErrorIO('SUTRA solution: Internal error testing if U solution &
(ML-1 must be -1, 0, 1) ML-1 = ['//trim(adjustl(Val2Char(ML-1)))//']')
         END SELECT
3300 &
         SELECT CASE (NOUMAT)
          CASE(-1)
            GO TO 3350
          CASE(0)
            GO TO 3350
          CASE(1)
            GO TO 3375
          CASE DEFAULT
            MSErrorValue%cDataSet='SUT'
            call ErrorIO('SUTRA solution: Internal error (NOUMAT must be -1, 0, 1) &
NOUMAT = ['//trim(adjustl(Val2Char(NOUMAT)))//']')
         END SELECT
 3350    CALL ZERO(UMAT, MATDIM, 0.0D0) 
 3375    CALL ZERO(CVEC, NN, 0.0D0) 
!         FORALL (I = 1:NN, NK = KSP:KSP) 
!           UVEC (I, NK) = CVEC (I) 
!         ENDFORALL 
         DO I = 1, NN
           UVEC(I, KSP) = CVEC(I) 
         END DO 
 3400    CONTINUE 
!                                                                       
!.....SET TIME-DEPENDENT BOUNDARY CONDITIONS, SOURCES AND SINKS         
!        FOR THIS TIME STEP                                             
         IF (ITER.EQ.1.AND.IBCT.NE.4) &
           CALL BCTIME (IPBCT, IUBCT, IQSOPT, IQSOUT)                                                          
!                                                                       
!.....SET SORPTION PARAMETERS FOR THIS TIME STEP                        
      IF (U_ML(K).NE.1.AND.ME.LT.1.AND.NOUMAT.EQ.0.AND.TRIM(ADSMOD(K)).NE.'NONE') &
        CALL ADSORB ()                 
!      IF ( U_ML(K).NE.1.AND.ME.LT.1.AND.NOUMAT.EQ.0 ) CALL ADSORB()                 
!                                                                       
!.....DO ELEMENTWISE CALCULATIONS IN MATRIX EQUATION FOR P AND/OR U     
         CALL CPU_TIME (RPTIME0) 
         IF (NOUMAT.EQ.0) THEN 
            IF (IABS (KTYPE) .EQ.3) THEN 
!..... 3-D PROBLEM.                                                     
               CALL ELEMN3 (U_ML(K))                    
            ELSE 
!..... 2-D PROBLEM.                                                     
               CALL ELEMN2 (U_ML(K))        
!......INITIALIZE VANG2 TO ZERO FOR 2-D PROBLEMS                        
               VANG2 = 0.0D0 
            ENDIF 
         ENDIF 
         CALL CPU_TIME (RPTIME1) 
         WRITE ( * , '(1X,A,1X,A,A)') &
                     '....numerical integration     (', trim(adjustl(cTime(RPTIME1 - RPTIME0))) , ')'
!                                                                       
!.....DO NODEWISE CALCULATIONS IN MATRIX EQUATION FOR P AND/OR U        
         CALL NODAL( U_ML(K) )
!                                                                       
!.....SET SPECIFIED P AND U CONDITIONS IN MATRIX EQUATION               
!     FOR P AND/OR U                                                    
         CALL BC( U_ML(K) )
!                                                                       
         CALL CPU_TIME (RPTIME1) 
         WRITE ( * , '(1x,a,1x,a,a)' ) '....matrix assembly completed (',trim(adjustl(cTime(RPTIME1 - RPTIME0))), ')'
         WRITE (fSMY, * ) '....matrix assembly completed' 
!                                                                       
!.....MATRIX EQUATION FOR P AND/OR U ARE COMPLETE, SOLVE EQUATIONS:     
!     WITH DIRECT SOLVER,                                               
!        WHEN KKK=0, DECOMPOSE AND BACK-SUBSTITUTE,                     
!        WHEN KKK=2, BACK-SUBSTITUTE ONLY.                              
!     WITH SLAP SOLVER,                                                 
!        WHEN KKK=0, RESET MATRIX POINTERS TO "TRIAD" FORMAT,           
!        WHEN KKK=2, LEAVE MATRIX POINTERS IN "COLUMN" FORMAT.          
!        KPU=1 WHEN SOLVING FOR P,                                      
!        KPU=2 WHEN SOLVING FOR U.                                      
         IHALFB = NBHALF - 1 
!
!.....SOLVE FOR P                                                       
         SELECT CASE (U_ML(K)-1)
           CASE(-1)
             GO TO 5000
           CASE(0)
             GO TO 5000
           CASE(1)
             GO TO 5500
           CASE DEFAULT
             MSErrorValue%cDataSet='SUT'
             call ErrorIO('SUTRA solution: Internal error initializing P solution &
(ML-1 must be -1, 0, 1) ML-1 = ['//trim(adjustl(Val2Char(ML-1)))//']')
         END SELECT
!
 5000    KKK = 000000 
         KPU = 1 
         WRITE ( * , * ) '....SOLVE FOR P' 
!                                                                       
         CALL CPU_TIME (RPTIME0) 
!                                                                       
         if(.not.SOLVER (KKK, KPU, KSOLVP, PMAT, PVEC, PITER, B, NN,       &
                         IHALFB, NELT, NCBI)) then
           ISTOP_OUTER = -1
           IGOI_OUTER  =  1
         end if
!                                                                       
         CALL CPU_TIME (RPTIME1) 
         WRITE ( * , '(1X,6x,A,1X,A,1X,A)') &
                     'iteration time          (',trim(adjustl(cTime(RPTIME1 - RPTIME0))),')'
!                                                                       
!.....P SOLUTION NOW IN PVEC                                            
!                                                                       
!.....SOLVE FOR U                                                       
         SELECT CASE (U_ML(K)-1)
          CASE(-1)
            GO TO 5500
          CASE(0)
            GO TO 5750
          CASE(1)
            GO TO 5500
          CASE DEFAULT
            MSErrorValue%cDataSet='SUT'
            call ErrorIO('SUTRA solution: Internal error initializing U solution &
(ML-1 must be -1, 0, 1) ML-1 = ['//trim(adjustl(Val2Char(ML-1)))//']')
         END SELECT
!                                                                       
 5500    KKK = 000000 
         KPU = 2 
         WRITE ( * , '(1x,a)' ) '....SOLVE FOR U['//trim(adjustl(Val2Char(K)))//'] - '//trim(adjustl(SPNAME(K)))
!                                                                       
!.....SET CVEC TO UVEC VALUE FOR SPECIES                                
!.....SET CITER TO UITER VALUE FOR SPECIES                              
!         FORALL (I = 1:NN) 
!           CVEC(I)  = UVEC(I, KSP) 
!           CITER(I) = UITER(I, KSP) 
!         ENDFORALL 
         DO I = 1, NN
           CVEC(I)  = UVEC(I, KSP) 
           CITER(I) = UITER(I, KSP) 
         END DO 

         SELECT CASE (NOUMAT)
          CASE(-1)
            KKK=000000
          CASE(0)
            KKK=000000
          CASE(1)
            KKK=2
          CASE DEFAULT
            MSErrorValue%cDataSet='SUT'
            call ErrorIO('SUTRA solution: NOUMAT Error (NOUMAT must be -1, 0, 1) NOUMAT = ['//trim(adjustl(Val2Char(NOUMAT)))//']')
         END SELECT
!                                                                       
         CALL CPU_TIME (RPTIME0) 
!                                                                       
!         if ( ksp.eq.1 .and. IT.le.255 ) then
!           write (1069,'(2i10,g15.6)') IT, ITER, TSEC
!           write (1069,*) 'CVEC'
!           write (1069,'(100(10(g15.6),/))') CVEC
!           write (1069,*) 'CITER'
!           write (1069,'(100(10(g15.6),/))') CITER
!           write (1069,*) 'UMAT'
!           write (1069,'(100(10(g15.6),/))') UMAT
!           write (1069,*) ' '
!         end if
         if ( .not.SOLVER(KKK, KPU, KSOLVU(KSP), UMAT, CVEC, CITER, B, NN, &
                          IHALFB, NELT, NCBI)) then
           ISTOP_OUTER = -1
           if( IGOI_OUTER==1 ) then
             IGOI_OUTER=3
           else
             IGOI_OUTER=2
           end if
         end if
!                                                                       
         CALL CPU_TIME (RPTIME1) 
         WRITE ( * , '(1X,6x,A,1X,A,1X,A)') &
                  'iteration time          (',trim(adjustl(cTime(RPTIME1 - RPTIME0))),')'
!                                                                       
!.....SET UVEC FOR SPECIES TO CVEC VALUE                                
!.....SET UITER FOR SPECIES TO CITER VALUE                              
!         FORALL (I = 1:NN, NK = KSP:KSP) 
!           UVEC (I, NK) = CVEC (I) 
!           UITER (I, NK) = CITER (I) 
!         ENDFORALL 
         DO I = 1, NN
           UVEC(I, KSP)  = CVEC(I) 
           UITER(I, KSP) = CITER(I) 
         END DO
!.....U SOLUTION NOW IN UVEC                                            
!                                                                       
!.....PREVENT MULTIPLE SOLUTIONS OF PRESSURE SOLUTION FOR MULTI SPECIES SIMULATIONS                
!     ML = 0 FOR P AND U                                                
!     ML = 1 FOR P ONLY                                                 
!     ML = 2 FOR U ONLY                                                 
 5750   IF (NSPE.GT.1) THEN 
        !save ML value for use on last species for this iteration       
            IF (K.EQ.1) THEN 
!               IF (U_ML(K).EQ.0) THEN 
               IF (ML.EQ.0) THEN 
                  U_ML(2:NSPE)=2
               ENDIF 
            ENDIF 
         ENDIF 
!                                                                       
 5800 END DO 
!.....END OF SPECIES LOOP                                               
!                                                                       
 5900 CONTINUE 
!                                                                       
!.....CHECK PROGRESS AND CONVERGENCE OF ITERATIONS                      
!        AND SET STOP AND GO FLAGS:                                     
!           ISTOP = -1   NOT CONVERGED - STOP SIMULATION                
!           ISTOP =  0   ITERATIONS LEFT OR CONVERGED - KEEP SIMULATING 
!           ISTOP =  1   LAST TIME STEP REACHED - STOP SIMULATION       
!           ISTOP =  2   MAXIMUM TIME REACHED - STOP SIMULATION         
!           IGOI = 0   P AND U CONVERGED, OR NO ITERATIONS DONE         
!           IGOI = 1   ONLY P HAS NOT YET CONVERGED TO CRITERION        
!           IGOI = 2   ONLY U HAS NOT YET CONVERGED TO CRITERION        
!           IGOI = 3   BOTH P AND U HAVE NOT YET CONVERGED TO CRITERIA  
      ISTOP = 0 
      IGOI = 0 
      SELECT CASE (ITRMAX-1)
        CASE(:-1)
          GO TO 7500
        CASE(0)
          GO TO 7500
        CASE(1:)
          GO TO 7000
       END SELECT
 7000 RPM = 0.D0 
      RP = 0.D0 
      RUM = 0.D0 
      RU = 0.D0 
      IPWORS = 0 
      IUWORS = 0 
      NOCONU = 0 
      SELECT CASE (ML-1)
        CASE(-1)
          GO TO 7050
        CASE(0)
          GO TO 7050
        CASE(1)
          GO TO 7150
        CASE DEFAULT
          MSErrorValue%cDataSet='SUT'
          call ErrorIO('SUTRA solution: Internal error testing for P convergence &
(ML-1 must be -1, 0, 1) ML-1 = ['//trim(adjustl(Val2Char(ML-1)))//']')
       END SELECT
 7050 DO 7100 I = 1, NN 
         RP = DABS( PVEC(I) - PITER(I) ) 
         IF (RP - RPM) 7100, 7060, 7060 
 7060    RPM = RP 
         IPWORS = I 
 7100 END DO 
      IF ( RPM.GT.RPMAX ) IGOI = IGOI + 1 
!.....TEST IF U SOLUTION THIS TIME STEP                                 
7150 &
      SELECT CASE (ML-1)
        CASE(-1)
          GO TO 7200
        CASE(0)
          GO TO 7350
        CASE(1)
          GO TO 7200
        CASE DEFAULT
          MSErrorValue%cDataSet='SUT'
          call ErrorIO('SUTRA solution: Internal error testing for U convergence &
(ML-1 must be -1, 0, 1) ML-1 = ['//trim(adjustl(Val2Char(ML-1)))//']')
       END SELECT
 7200 DO 7300 I = 1, NN 
         DO 7270 K = 1, NSPE 
            RU(K) = DABS( UVEC(I, K) - UITER(I, K) ) 
            IF ( RU(K) - RUM(K) ) 7300, 7260, 7260 
 7260       RUM(K) = RU(K) 
            IUWORS(K) = I 
 7270    END DO 
 7300 END DO 
      !WRITE INFORMATION ABOUT PICARD ITERATIONS
      WRITE ( * , '(/3x,53("*")/,6xa/,6xa20,2(1x,a15))') &
        'PICARD ITERATION INFORMATION', &
        'SPECIES             ', 'CRITERIA       ', 'RESULT         '
      WRITE ( * , '(6x,a10,2(1x,g16.7))' ) &
        'PRESSURE            ', RPMAX, RPM
      DO 7310 K = 1, NSPE 
        WRITE ( * , '(6x,a10,2(1x,g16.7))' ) &
          SPNAME(K), RUMAX(K), RUM(K)
         IF (RUM (K) .GT.RUMAX (K) ) THEN 
            IGOI = IGOI + 2 
            NOCONU (K) = 1 
         ENDIF 
 7310 END DO 
      WRITE ( * , '(3x,53("*"))')
!                                                                       
 7350 CONTINUE 
      IF (IGOI.GT.0.AND.ITER.EQ.ITRMAX) ISTOP = - 1 
      IF (IGOI.GT.0.AND.ISTOP.EQ.0) GOTO 1100 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....END ITERATION - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
 7500 CONTINUE 
!                                                                       
      IF( ISTOP_OUTER == -1) ISTOP = -1
      IF(  IGOI_OUTER  >  0) IGOI  = IGOI_OUTER
      
      IF (ISTOP.NE. - 1.AND.IT.EQ.ITMAX) ISTOP = 1 
      IF (ISTOP.NE. - 1.AND.TSEC.GE.TMAX) ISTOP = 2 
!                                                                       
!.....AUTOMATIC TIMESTEP                                                
!      IF (fATS.GT.0) THEN 
      IF ( LATS ) THEN 
!.......If imaxiter exceeded imaxtarget during the last time step       
!       reduce delt                                                     
!........Initialize logical flag which determines if the timestep length can be increased
         LATSAllowTimeStepIncrease =.true.
!........Initialize logical flag which determines if the ATS convergence criteria has been exceeded
         LATSFailConvergence=.false.
!........SET logical flag to prevent the timestep length from being increased during the next timestep
         IF ( imaxiter > (imaxtarget/2) ) LATSAllowTimeStepIncrease = .false.                                                          
         IF (imaxiter>imaxtarget .and. DELT>dtmin) then
            do1 = DELT 
!...........Reduce DELT by two DTMULT increments                          
            DELT = DELT / DTMULT 
            DELT = DELT / DTMULT 
!...........reset DELT to dtmin if required
            IF (DELT<dtmin) DELT = dtmin 
            IF (fTBC>0.or.fOTM>0) DELTSAVE = DELT 
            WRITE ( * , '(1x,57("*")/,(1xa,1xi5/),3(1x,a,1x,a/),2(1x,a,1x,i5/),1x,57("*"))') &
              'FAILED USER-DEFINED CONVERGENCE CRITERIA ON TIMESTEP', IT, &
              'current total time:', trim(adjustl(cTime(TSEC))),&
              ' previous timestep:', trim(adjustl(cTime(do1))), &
              '      new timestep:', trim(adjustl(cTime(DELT))), &
              '           maxiter:', imaxiter, '        imaxtarget:', imaxtarget              
            WRITE (fLST, '(1x,57("*")/,(1xa,1xi5/),3(1x,a,1x,a/),2(1x,a,1x,i5/),1x,57("*"))') &
              'FAILED USER-DEFINED CONVERGENCE CRITERIA ON TIMESTEP', IT, &
              'current total time:', trim(adjustl(cTime(TSEC))), &
              ' previous timestep:', trim(adjustl(cTime(do1))), &
              '      new timestep:', trim(adjustl(cTime(DELT))), &
              '           maxiter:', imaxiter, '        imaxtarget:', imaxtarget                
!...........RESET imaxiter COUNTER                                        
            imaxiter = 0 
!...........RESET THE APPROPRIATE VALUES FOR P AND U                      
!           AND RESET APPROPRIATE ARRAYS                                  
            IF ( .not.LContinue.and.NMaxCount<NConvMax ) then 
               NMaxCount = NMaxCount + 1 
               !WRITE ADDITIONAL INFORMATION ABOUT RERUN OF TIMESTEP       
               WRITE ( * , '(1x,57("*")/,(1xa,1xi5/),1(1x,a,1x,i5/),1x,57("*"))') &
                      'RERUNNING TIMESTEP', IT, '        NMaxCount:', NMaxCount                                 
               WRITE (fLST, '(1x,57("*")/,(1xa,1xi5/),1(1x,a,1x,i5/),1x,57("*"))') &
                      'RERUNNING TIMESTEP', IT, '        NMaxCount:', NMaxCount                                   
!..............RESET ITER COUNTER                                          
               ITER = 0 
               imaxiter = 0 
!..............REWIND TO INTIAL TIME (TSEC)                                
               TSEC = TSEC - do1 
!..............RESET APPROPRIATE VALUES FOR CURRENT SOLUTION SCHEME        
               selectcase (ML - 1) 
                !P AND U                                                
                 case ( - 1) 
                   TSECP0 = TSEC 
                   TSECU0 = TSEC 
                   DELTP  = DELTP_ 
                   DELTU  = DELTU_ 
                   PVEC   = PVEC_ 
                   UVEC   = UVEC_ 
                   PM1    = PM1_ 
                   UM1    = UM1_ 
                   RCIT   = RCIT_ 
                   CS1    = CS1_ 
                   SW     = SW_ 
                   QPLITR = QPLITR_
                !P ONLY                                                 
                 case (0) 
                   TSECP0 = TSEC 
                   DELTP  = DELTP_ 
                   PVEC   = PVEC_ 
                   PM1    = PM1_ 
                   RCIT   = RCIT_ 
                   SW     = SW_ 
                   QPLITR = QPLITR_
                !U ONLY                                                 
                 case (1) 
                   TSECU0 = TSEC 
                   DELTU  = DELTU_ 
                   UVEC   = UVEC_ 
                   UM1    = UM1_ 
                   CS1    = CS1_ 
               endselect
               LATSFailConvergence=.true.
!..............RERUN THIS TIMESTEP                                         
               GOTO 1005 
            ELSE 
               NMaxCount = 0 
            ENDIF 
         ELSE 
            imaxiter = 0 
            NMaxCount = 0 
         ENDIF 
      ENDIF 
!                                                                       
!.....OUTPUT RESULTS FOR TIME STEP IN ACCORDANCE WITH PRINT CYCLES      
      PRNALL = ( (IT.LE.1) .OR. (ISTOP.NE.0) ) 
      IF(fOTM.LT.1) PRNLST = (PRNALL.OR. (MOD (IT, NPRINT) .EQ.0) ) 
      IF (NCOLPR.GT.0 .AND. fOTM.LT.1) PRNNOD = (PRNALL.OR. (MOD (IT, NCOLPR) .EQ.0) ) 
      IF (LCOLPR.GT.0 .AND. fOTM.LT.1) PRNELE = (PRNALL.OR. (MOD (IT, LCOLPR) .EQ.0) ) 
      ObsTimeDiff=EPSILON
      if(allocated(TIMOBS)) ObsTimeDiff=abs(TIMOBS(ONOBS)-TSEC)
      PRNOBS = (PRNALL.OR. (ObsTimeDiff.LT.EPSILON) ) 
      if(PRNOBS) ONOBS=ONOBS+1
!.....VERSION 1.1
!.....TECPLOT PRINTING                                                  
      IF (NTPNP.LT.1) NTPNP = ITMAX 
      IF (NTPEP.LT.1) NTPEP = ITMAX 
      PRNTPN = (PRNALL.OR. (MOD (IT, NTPNP) .EQ.0) ) 
      PRNTPE = (PRNALL.OR. (MOD (IT, NTPEP) .EQ.0) ) 
!                                                                       
!.....CONTROLLED OUTPUT PRINTING                                        
!     TO fLST, fNOD, AND fELE OUTPUT FILES                               
      IF (fOTM.GT.0) THEN 
         PRNLST  = .FALSE. 
         PRNNOD  = .FALSE. 
         PRNELE  = .FALSE. 
!........VERSION 1.1
         PRNTPN  = .FALSE. 
         PRNTPE  = .FALSE. 
         !OUTPUT ON LAST TIMESTEP OR ITERATION
         ! REGARDLESS OF CONVERGENCE OR WHETHER IT IS SPECIFIED IN CONTROLLED OUTPUT PRINT FILE (OTM)
         IF (ABS(TSEC-DPTIME).LT.EPSILON .OR. ISTOP.NE.0) THEN 
            IF ( fLST.GT.0) PRNLST = .TRUE. 
            IF ( fNOD.GT.0) PRNNOD = .TRUE. 
            IF ( fELE.GT.0) PRNELE = .TRUE. 
            IF ( fTPN.GT.0) PRNTPN = .TRUE. 
            IF ( fTPE.GT.0) PRNTPE = .TRUE. 
         ENDIF 
      ENDIF 
!                                                                       
      IF (PRNLST) THEN 
         itemp=size(U_ML)
         IF (IABS (KTYPE) .EQ.3) THEN 
            CALL OUTPTCLST_3D(ML, ISTOP, IGOI, PVEC, UVEC, VMAG, VANG1, VANG2, SW)                                                  
         ELSE 
            CALL OUTPTCLST_2D(ML, ISTOP, IGOI, PVEC, UVEC, VMAG, VANG1, SW) 
         ENDIF 
!                                                                       
!.......RESET VALUES FOR ENERGY TRANSPORT FOR BUDGET CALCULATIONS       
         CS = CST 
         CW = CWT 
!                                                                       
!.....CALCULATE AND PRINT FLUID MASS AND/OR ENERGY OR                   
!     SOLUTE MASS BUDGET                                                
         IF (KBUDG.EQ.1) &
           CALL BUDGET (ML, IBCT, VOL, SW, DSWDP, RHO,    &
                        QIN, PVEC, PM1, PBC, QPLITR, IPBC, IQSOP, UVEC, UM1, &
                        UM2, UIN, QUIN, IQSOU, UBC, IUBC, CS1, CS2, CS3, SL, SR) 
      ENDIF 
!                                                                       
      IF (PRNNOD) CALL OUTNODE (PVEC, UVEC, SW, IN, X, Y, Z, TITLE1, TITLE2)                                                           
      IF (PRNELE) CALL OUTVELOCITY(VMAG, VANG1, VANG2, IN, X, Y, Z, TITLE1, TITLE2)                                                           
      IF (PRNOBS) CALL OUTOBS(IOBS, X, Y, Z, PVEC, UVEC, SW, TITLE1, TITLE2)                                                           
!.....VERSION 1.1                                                                       
!.....PRINT TECPLOT RESULTS                                             
      IF (PRNTPN .and. fTPN>0) then
        if (.not. PrintTecplotNodeData()) &
          call ErrorIO('SUTRA:: Error writing Tecplot node data')
      end if
      IF (PRNTPE .and. fTPE>0) THEN
        if(.not. PrintTecplotElementData() ) &
          call ErrorIO('SUTRA:: Error writing Tecplot element data')
      END IF
!                                                                       
!.....STORE RESULTS FOR POSSIBLE RESTART OF SIMULATION EACH             
!        ISTORE TIME STEPS AND AFTER LAST TIME STEP                     
      IF (ISTORE.EQ.0) GOTO 8150 
      IF (MOD (IT, ISTORE) .NE.0.AND.ISTOP.EQ.0) GOTO 8150 
      CALL STORE (PVEC, UVEC, PM1, UM1, CS1, RCIT, SW, PBC) 
!                                                                       
!.....RESTORE DELT AFTER MODIFICATION FOR EXACT OUTPUT TIMES                      
 8150 IF (fOTM.GT.0.AND.DELT.LT.DELTSAVE) THEN
        DELT = DELTSAVE 
      END IF
!
!.....RESTORE DELT AFTER MODIFICATION FOR OBS OUTPUT
      IF (lTIMOBSRESET) THEN
        DELT = DELTSAVE
        lTIMOBSRESET=.FALSE.
      END IF
!                                                                       
!.....READ TRANSIENT BOUNDARY CONDITIONS, SOURCES AND SINKS             
!     FOR THIS NEXT TIME STEP, IF REQUIRED.                                  
       IF (lGetTBC) THEN 
!........WRITE TBCs INFO TO SCREEN AND fLST
           WRITE ( *,'(//1X,"READING TRANSIENT B.C.s BEGINNING AT TIME ",A)') &
             trim(adjustl(cTime(DNTIME))) 
           WRITE (fLST,'(//1X,"READING TRANSIENT B.C.s BEGINNING AT TIME ",A)') &
             trim(adjustl(cTime(DNTIME))) 
!........READ REVISED TRANSIENT BOUNDARY CONDITIONS                      
           lRdTBCTime=.FALSE.
           if(.not. ReadTBCData()) call ErrorIO('SUTRA: Error reading transient BC data')
!.......READ NEXT TIME AT WHICH TO CHANGE B.C.s                         
           lRdTBCTime=.TRUE.
           if(.not. ReadTBCData()) call ErrorIO('SUTRA: Error reading transient BC data')
       ENDIF 
!
!.....GET NEXT PRINT TIME
      if (fOTM.GT.0.AND.DPTIME.EQ.TSEC.AND.TSEC.LT.TMAX) then
        if (.not. ReadUserSpecifiedOutputTime() ) call ErrorIO('SUTRA:: Error reading next print time')
      end if
!
      IF (ISTOP.EQ.0) GOTO 1000 
! **********************************************************************
!.....END TIME STEP ****************************************************
! **********************************************************************
!                                                                       
!                                                                       
!.....COMPLETE OUTPUT AND TERMINATE SIMULATION                          
!                                                                       
      WRITE ( * ,  * ) 'S I M U L A T I O N   E N D E D' 
      WRITE (fSMY,  * ) 'S I M U L A T I O N   E N D E D' 
      IF (ISTORE.GT.0 .and. fRST>0) WRITE (fLST, 8100) 
 8100 FORMAT(//////11X,'*** LAST SOLUTION HAS BEEN STORED ON UNIT fRST ***')                                              
!                                                                       
!.....OUTPUT END OF SIMULATION MESSAGE AND RETURN TO MAIN FOR STOP      
      MSErrorValue%cDataSet='SUT'
      IF (ISTOP.GT.0) GOTO 8400 
      write(*,*) IGOI
      SELECT CASE (IGOI - 2)
        CASE(:-1)
          GO TO 8230
        CASE(0)
          GO TO 8260
        CASE(1:)
          GO TO 8290
        CASE DEFAULT
          MSErrorValue%cDataSet='SUT'
          call ErrorIO('SUTRA solution: Internal error setting IGOI-2 ['//trim(adjustl(Val2Char(IGOI-2)))//']')
       END SELECT

!P solution failed 
 8230 call ErrorIO('Error: Simulation terminated due to non-convergent pressure solution')
!P and U solution failed
 8260 select case(ME)
        case(-1:0)
          call ErrorIO('Error: Simulation terminated due to non-convergent solute transport solution')
        case(1)
          call ErrorIO('Error: Simulation terminated due to non-convergent heat transport solution')
      end select
!U solution failed
8290  select case(ME)
        case(-1)
          call ErrorIO('Error: Simulation terminated due to non-convergent pressure and solute transport solution')
        case(0)
          call ErrorIO('Error: Simulation terminated due to non-convergent pressure, solute transport, and heat transport solution')
        case(1)
          call ErrorIO('Error: Simulation terminated due to non-convergent pressure and heat transport solution')
      end select
                                                                       
08400 select case(ISTOP)
        case(2)
          WRITE (fLST, 8550) 
        case default
          WRITE (fLST, 8450) 
      end select

 8450 FORMAT(////////11X,'SUTRA SIMULATION TERMINATED AT COMPLETION OF TIME STEPS'/                                               &
     &               11X,'***** ********** ********** ** ********** ** **** *****')                                               
 8550 FORMAT(////////11X,'SUTRA SIMULATION TERMINATED AT COMPLETION OF TIME PERIOD'/                                              &
     &               11X,'***** ********** ********** ** ********** ** **** ******')                                              

      !RETURN TO MAIN SUTRA PROGRAM
      RETURN 
!                                                                       
      END SUBROUTINE SUTRA                          
