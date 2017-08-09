!     SUBROUTINE        I  N  D  A  T  2       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO READ INITIAL CONDITIONS FROM UNIT-fICS, AND TO                  
! ***  INITIALIZE DATA FOR EITHER WARM OR COLD START OF                 
! ***  THE SIMULATION.                                                  
!                                                                       
      SUBROUTINE INDAT2 (PVEC, UVEC, PM1, UM1, UM2, CS1, CS2, CS3, SL,  &
                         SR, RCIT, SW, DSWDP, PBC, IPBC, IPBCT)                      
      USE PARAMS 
      USE CONTRL 
      USE FUNITS 
      USE DIMS
      USE TIMES
      USE KPRINT
      USE SutraStorage, ONLY : SpecifiedPBC
      use SutraMSPrecision 
      !MS specific modules
      USE MSErrorHandler
      USE SutraZoneModule
      implicit none
      CHARACTER(LEN=10) CPUNI, CUUNI 
      CHARACTER (LEN = 15) :: CTSPE, CSPE 
      real (DP) :: &
        PVEC(NN),&
        UVEC(NN, NSPE),&
        PM1(NN),&
        UM1(NN,NSPE),&
        UM2(NN,NSPE),&
        SL(NN),&
        SR(NN),&
        CS1(NN,NSPE),&
        CS2(NN),&
        CS3(NN),&
        RCIT(NN),&
        SW(NN),&
        DSWDP(NN),&
        PBC(NBCN)                                                       
      integer (I4B) :: &
        IPBC(NBCN)
      !locals
      integer (I4B) :: &
        NLSKIP, &
        IPU, &
        IPBCT, &
        IP
      integer (I4B) :: &
        I, &
        K, &
        KK
      real (DP) :: &
        PUNI, &
        UUNI, &
        RELK
      real (DP) :: &
        tv

!                                                                       
!                                                                       
      MSErrorValue%cDataSet='ICS'
      CALL SKPCOM (fICS, NLSKIP) 
      IF (IREAD) 500, 500, 620 
!.....INPUT INITIAL CONDITIONS FOR WARM START (UNIT-fICS DATA)            
  500 READ (fICS,*,iostat=ios) TSTART, DELTP, DELTU 
      if(ios/=0) call ErrorIO('Error specifying warm start temporal information')
      CALL SKPCOM (fICS, NLSKIP) 
      READ (fICS,*,iostat=ios) CPUNI 
      if(ios/=0) call ErrorIO('Error specifying CPUNI for warm start')
      IF (CPUNI.NE.'NONUNIFORM') THEN 
         call ErrorIO('Error: CPUNI must be "NONUNIFORM" for warm start')
      ENDIF 
      READ (fICS,*,iostat=ios) (PVEC (I), I = 1, NN) 
      if(ios/=0) call ErrorIO('Error specifying PVEC values for warm start')
      CALL SKPCOM (fICS, NLSKIP) 
      READ (fICS,*,iostat=ios) CUUNI 
      if(ios/=0) call ErrorIO('Error specifying CUUNI for warm start')
      IF (CUUNI.NE.'NONUNIFORM') THEN 
         call ErrorIO('Error: CUUNI must be "NONUNIFORM" for warm start')
      ENDIF 
      CALL SKPCOM (fICS, NLSKIP) 
      DO K = 1, NSPE 
        READ (fICS,*,iostat=ios) (UVEC (I, K), I = 1, NN) 
        if(ios/=0) call ErrorIO('Error specifying UVEC('//trim(adjustl(Val2Char(K)))//') values for warm start')
      ENDDO 
      READ (fICS,*,iostat=ios) (PM1 (I), I = 1, NN) 
      if(ios/=0) call ErrorIO('Error specifying PVEC values for warm start')
      DO K = 1, NSPE 
        READ (fICS,*,iostat=ios) (UM1 (I, K), I = 1, NN) 
        if(ios/=0) call ErrorIO('Error specifying UM1('//trim(adjustl(Val2Char(K)))//') values for warm start')
      ENDDO 
      DO K = 1, NSPE 
        READ (fICS,*,iostat=ios) (CS1 (I, K), I = 1, NN) 
        if(ios/=0) call ErrorIO('Error specifying CS1('//trim(adjustl(Val2Char(K)))//') values for warm start')
      ENDDO 
      READ (fICS,*,iostat=ios) (RCIT (I), I = 1, NN) 
      if(ios/=0) call ErrorIO('Error specifying RCIT values for warm start')
      READ (fICS,*,iostat=ios) (SW (I), I = 1, NN) 
      if(ios/=0) call ErrorIO('Error specifying SW values for warm start')
      READ (fICS,*,iostat=ios) (SpecifiedPBC(IPU)%P, IPU = 1, NPBC) 
      if(ios/=0) call ErrorIO('Error specifying PBC values for warm start')
      CALL ZERO (SL, NN, 0.0D0) 
      CALL ZERO (SR, NN, 0.0D0) 
      CALL ZERO (DSWDP, NN, 0.0D0) 
!.....USE CONFORMAL ARRAY PROPERTIES TO SET UM2 EQUAL TO UM1            
      UM2 = UM1 
!                                                                       
      GOTO 1000 
!                                                                       
!.....INPUT INITIAL CONDITIONS FOR COLD START (UNIT-fICS DATA)            
  620 READ (fICS,*,iostat=ios) TSTART 
      if(ios/=0) call ErrorIO('Error specifying TSTART for cold start')
      CALL SKPCOM (fICS, NLSKIP) 
      READ (fICS,*,iostat=ios) CPUNI 
      if(ios/=0) call ErrorIO('Error specifying CPUNI for cold start')
      IF (CPUNI.EQ.'UNIFORM') THEN 
         READ (fICS,*,iostat=ios) PUNI 
         if(ios/=0) call ErrorIO('Error specifying PUNI for cold start')
         DO 625 I = 1, NN 
            PVEC (I) = PUNI 
  625    END DO 
      ELSEIF (CPUNI.EQ.'NONUNIFORM') THEN 
         READ (fICS,*,iostat=ios) (PVEC (I), I = 1, NN) 
         if(ios/=0) call ErrorIO('Error specifying PVEC for cold start')
      ELSE 
         call ErrorIO('Error specifying CPUNI must be "UNIFORM" or "NONUNIFORM" for cold start')
      ENDIF 
      DO 635 K = 1, NSPE 
         CALL SKPCOM (fICS, NLSKIP) 
         READ (fICS,*,iostat=ios) CUUNI 
         if(ios/=0) call ErrorIO('Error specifying CPUNI for cold start')
!.......DETERMINE WHICH SPECIES IF MULTI-SPECIES SIMULATION             
         IF (NSPE.LT.2) THEN 
            KK = K 
         ELSE 
            BACKSPACE (fICS) 
            READ (fICS,*,iostat=ios) CUUNI, CTSPE 
            if(ios/=0) call ErrorIO('Error specifying CUUNI and/or species number for cold start')
            READ (CTSPE,*,iostat=ios) CSPE, KK 
            if(ios/=0) call ErrorIO('Error specifying CSPE or KK for cold start')
         ENDIF 
!.......READ INITIAL U DATA FOR THE SPECIES                             
         IF (CUUNI.EQ.'UNIFORM') THEN 
            READ (fICS,*,iostat=ios) UUNI 
            if(ios/=0) call ErrorIO('Error specifying UUNI for cold start')
            DO 630 I = 1, NN 
               UVEC (I, KK) = UUNI 
  630       END DO 
         ELSEIF (CUUNI.EQ.'NONUNIFORM') THEN 
            READ (fICS,*,iostat=ios) (UVEC (I, KK), I = 1, NN) 
            if(ios/=0) call ErrorIO('Error specifying UVEC('//trim(adjustl(Val2Char(KK)))//') values for cold start')
         ELSE 
            call ErrorIO('Error specifying CPUNI must be "UNIFORM" or "NONUNIFORM" for cold start')
         ENDIF 
  635 END DO 
!.....START-UP WITH NO PROJECTIONS BY SETTING BDELP=BDELU=1.D-16        
!     IN PROJECTION FORMULAE FOUND IN SUBROUTINE SUTRA.                 
      DELTP = DELT * 1.D16 
      DELTU = DELT * 1.D16 
!.....INITIALIZE RCIT TO RHOW0                                          
      CALL ZERO (RCIT, NN, RHOW0) 
!.....INITIALIZE SPECIFIED TIME-VARYING PRESSURES TO INITIAL PRESSURE   
!     VALUES FOR START-UP CALCULATION OF INFLOWS OR OUTFLOWS            
!     (SET QPLITR=0)                                                    
      IF (IPBCT) 680, 740, 740 
  680 DO 730 IP = 1, NPBC 
!         I = IPBC (IP) 
         I = SpecifiedPBC(IP)%node
         IF (I) 700, 700, 730 
  700    SpecifiedPBC(IP)%P = PVEC ( - I) 
  730 END DO 
!.....INITIALIZE P, U, AND CONSISTENT DENSITY                           
  740 DO 800 I = 1, NN 
         PM1 (I) = PVEC (I) 
         DO K = 1, NSPE 
           UM1 (I, K) = UVEC (I, K) 
           UM2 (I, K) = UVEC (I, K) 
           RCIT (I) = RCIT (I) + DRWDU (K) * (UVEC (I, K) - URHOW0 (K) ) 
         ENDDO 
  800 END DO 
!.....INITIALIZE SATURATION, SW(I)                                      
      CALL ZERO (SW, NN, 1.0D0) 
      CALL ZERO (DSWDP, NN, 0.0D0) 
      IF (IUNSAT.NE.1) GOTO 990 
      IUNSAT = 3 
      DO 900 I = 1, NN 
  900 IF (PVEC (I) .LT.0) CALL UNSAT (SW (I), DSWDP (I), RELK, PVEC (I), NodeMap(I) )
  990 CONTINUE
!      CALL ZERO (CS1, NN * NSPE, CS) 
      DO K = 1, NSPE 
        tv = 0.0D0
        IF ( K.EQ.NESP ) THEN
          tv = CS
        END IF
        CALL ZERO (CS1(:,K), NN, tv) 
      END DO
      CALL ZERO (SL, NN, 0.0D0) 
      CALL ZERO (SR, NN, 0.0D0) 
 1000 CONTINUE 
!                                                                       
!.....SET STARTING TIME OF SIMULATION CLOCK, TSEC                       
      TSEC = TSTART 
!                                                                       
!.....RETURN TO CALLING ROUTINE                                         
      RETURN 
      END SUBROUTINE INDAT2                         
