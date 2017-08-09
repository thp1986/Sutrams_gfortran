!     SUBROUTINE        I  N  D  A  T  1       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO INPUT ,OUTPUT, AND ORGANIZE A MAJOR PORTION OF                
! ***  UNIT-fINP INPUT DATA (DATASETS 8 THROUGH 15)                       
!                                                                       
      SUBROUTINE INDAT1()
      USE ITERAT 
      USE PARAMS 
      USE CONTRL 
      USE JCOLS 
      USE MODSOR 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE OBS
      USE GRAVEC
      use SutraMSPrecision
      !MS specific modules
      USE MSErrorHandler
      USE SutraStorage
      USE SutraZoneModule
      USE SobIO  !specified observation nodes
      USE TECPLOT
      implicit none
      CHARACTER (LEN =10) :: CDUM10 
      CHARACTER (LEN =14) :: UTYPE (2) = (/' TEMPERATURES ','CONCENTRATIONS'/)
      CHARACTER (LEN = 6) :: STYPE (2) = (/'ENERGY','SOLUTE'/)
      CHARACTER (LEN = 1) :: CFIND, CSYM (NCOLMX) 
      CHARACTER (LEN = 1) :: CH_K5COL (NCOLMX) 
      CHARACTER (LEN = 1) :: CNODAL, CELMNT, CINCID, CVEL, CBUDG 

      !LOCAL VARIABLE
      LOGICAL :: lModifyLambdas
      character (len=128) :: cViscosityEquation
      character (len=128) :: cViscosity
      integer (I4B) :: &
        NLSKIP, &
        IME, &
        MONSPE, &
        NOBS, &
        NTOBS, &
        JSTOP, &
        NROLD, &
        NRTEST, &
        LRTEST, &
        LROLD, &
        MSTRUC, &
        IMULT, &
        ILAMB
      integer (I4B) :: &
        I, &
        II, &
        JJ, &
        L, &
        LL, &
        K, &
        KK, &
        M, &
        MM
      real (DP) :: &
        BSIGMAW, &
        BURHOW0, &
        BDRWDU, &
        SCALX, &
        SCALY, &
        SCALZ, &
        SCALTH, &
        PORFAC, &
        PMAXFA, &
        PMIDFA, &
        PMINFA, &
        ANG1FA, &
        ANG2FA, &
        ANG3FA, &
        ANGFAC, &
        ALMAXF, &
        ALMIDF, &
        ALMINF, &
        AT1MXF, &
        AT1MDF, &
        AT1MNF, &
        AT2MXF, &
        AT2MDF, &
        AT2MNF, &
        ATMAXF, &
        ATMINF, &
        D2R, &
        ZERO, &
        RADIAX, &
        SINA, &
        SINA2, &
        COSA, &
        COSA2

      real (DP), dimension(NSPE) :: sprod
      real (DP), dimension(NSPE) :: ssorp
!                                                                       
      SAVE UTYPE, STYPE 
!                                                                       
!.....INITIALIZE CSYM                                                   
      CSYM (1) = 'N' 
      CSYM (2) = 'X' 
      CSYM (3) = 'Y' 
      CSYM (4) = 'Z' 
      CSYM (5) = 'P' 
      DO K = 6, (5 + NSPE) 
      CSYM (K) = 'U' 
      ENDDO 
      CSYM (6 + NSPE) = 'S' 
      CSYM (7 + NSPE) = 'R' 
!.....VERSION 1.1
      CSYM (8 + NSPE) = 'H' 
      CSYM (9 + NSPE) = 'V' 
!      CSYM (10 + NSPE) = 'K' 	  		!HYDCON
!                                                                       
!.....INPUT DATASET 8A: OUTPUT CONTROLS AND OPTIONS -- UNIT fLST          
      MSErrorValue%cDataSet=' 8A'
      CALL SKPCOM (fINP, NLSKIP) 
      READ (fINP,*,iostat=ios) NPRINT, CNODAL, CELMNT, CINCID, CVEL, CBUDG 
      if(ios/=0) call ErrorIO('Error specifying output controls and options')
      IF (CNODAL.EQ.'Y') THEN 
         KNODAL = + 1 
      ELSEIF (CNODAL.EQ.'N') THEN 
         KNODAL = 0 
      ELSE 
         call ErrorIO('Error CNODAL must be "Y" or "N"')
      ENDIF 
      IF (CELMNT.EQ.'Y') THEN 
         KELMNT = + 1 
      ELSEIF (CELMNT.EQ.'N') THEN 
         KELMNT = 0 
      ELSE 
         call ErrorIO('Error CELMNT must be "Y" or "N"')
      ENDIF 
      IF (CINCID.EQ.'Y') THEN 
         KINCID = + 1 
      ELSEIF (CINCID.EQ.'N') THEN 
         KINCID = 0 
      ELSE 
         call ErrorIO('Error CINCID must be "Y" or "N"')
      ENDIF 
      IF (CVEL.EQ.'Y') THEN 
         KVEL = + 1 
      ELSEIF (CVEL.EQ.'N') THEN 
         KVEL = 0 
      ELSE 
         call ErrorIO('Error CVEL must be "Y" or "N"')
      ENDIF 
      IF (CBUDG.EQ.'Y') THEN 
         KBUDG = + 1 
      ELSEIF (CBUDG.EQ.'N') THEN 
         KBUDG = 0 
      ELSE 
         call ErrorIO('Error CBUDG must be "Y" or "N"')
      ENDIF 
!.....STANDARD SUTRA WRITE STATEMENTS
      WRITE (fLST, 172) NPRINT 
  172 FORMAT(////11X,'O U T P U T   C O N T R O L S   A N D   ',        &
     &   'O P T I O N S'//11X,'.lst FILE'/11X,'---------'               &
     &   //11X,I6,5X,'PRINTED OUTPUT CYCLE (IN TIME STEPS)')            
      IF (KNODAL.EQ. + 1) WRITE (fLST, 174) 
      IF (KNODAL.EQ.0) WRITE (fLST, 175) 
  174 FORMAT(/11X,'- PRINT NODE COORDINATES, THICKNESSES AND',          &
     &   ' POROSITIES')                                                 
  175 FORMAT(/11X,'- CANCEL PRINT OF NODE COORDINATES, THICKNESSES AND',&
     &   ' POROSITIES')                                                 
      IF (KELMNT.EQ. + 1) WRITE (fLST, 176) 
      IF (KELMNT.EQ.0) WRITE (fLST, 177) 
  176 FORMAT(11X,'- PRINT ELEMENT PERMEABILITIES AND DISPERSIVITIES') 
  177 FORMAT(11X,'- CANCEL PRINT OF ELEMENT PERMEABILITIES AND ',       &
     &   'DISPERSIVITIES')                                              
      IF (KINCID.EQ. + 1) WRITE (fLST, 178) 
      IF (KINCID.EQ.0) WRITE (fLST, 179) 
  178 FORMAT(11X,'- PRINT NODE INCIDENCES IN EACH ELEMENT') 
  179 FORMAT(11X,'- CANCEL PRINT OF NODE INCIDENCES IN EACH ELEMENT') 
      IME = 2 
      IF (ME.EQ. + 1) IME = 1 
      IF (KVEL.EQ. + 1) WRITE (fLST, 184) 
      IF (KVEL.EQ.0) WRITE (fLST, 185) 
  184 FORMAT(/11X,'- CALCULATE AND PRINT VELOCITIES AT ELEMENT ',       &
     &   'CENTROIDS ON EACH TIME STEP WITH OUTPUT')                     
  185 FORMAT(/11X,'- CANCEL PRINT OF VELOCITIES') 
      IF (KBUDG.EQ. + 1) WRITE (fLST, 186) STYPE (IME) 
      IF (KBUDG.EQ.0) WRITE (fLST, 187) 
  186 FORMAT(/11X,'- CALCULATE AND PRINT FLUID AND ',A6,' BUDGETS ',    &
     &   'ON EACH TIME STEP WITH OUTPUT')                               
  187 FORMAT(/11X,'- CANCEL PRINT OF BUDGETS') 
!                                                                       
!.....INPUT DATASET 8B: OUTPUT CONTROLS AND OPTIONS -- UNIT fNOD          
      MSErrorValue%cDataSet=' 8B'
      CALL SKPCOM (fINP, NLSKIP) 
      READ (fINP,*,iostat=ios) NCOLPR 
      if(ios/=0) call ErrorIO('Error specifying nodal output file controls (NCOLPR)')
!.....VERSION 1.1
!.....SET TECPLOT NODAL DATA OUTPUT TO HAVE THE SAME FREQUENCY          
!     A STANDARD SUTRA NODAL DATA OUTPUT.                               
!     TO OUTPUT TO A TECPLOT FILE AND NOT TO THE STANDARD SUTRA         
!     NODAL OUTPUT FILE SET NCOLPR TO THE DESIRED FREQUENCY AND         
!     DO NOT SPECIFY A fNOD FILE IN SUTRA.FIL                             
      IF (fTPN.GT.0) NTPNP = NCOLPR 
                                                                       
      MONSPE = 0 
      CSMSK5 = '                   ' 
!.....VERSION 1.1
      LTpP                     = .FALSE. 
      LTpU                     = .FALSE. 
      LTpSW                    = .FALSE.
      LTpRho                   = .FALSE.
      LTpTotalHead             = .FALSE.
!      LTpHydcon                = .FALSE.
      LTpVelocity              = .FALSE. 
      LTpPrintVelocityAtNodes  = .FALSE.
!      
      DO 1160 M = 1, NCOLMX 
         BACKSPACE (fINP) 
         READ(fINP,*,iostat=ios) NCOLPR,(CH_K5COL(MM), MM = 1, M) 
          if(ios/=0) call ErrorIO('Error specifying nodal output elements')
         MONSPE = MONSPE+1 
         K5COL(MONSPE) = CH_K5COL(M) 
         IF(CH_K5COL(M) .EQ.'U') THEN 
            CSMSK5(MONSPE) = trim(adjustl(SPNAME(1)))
            DO K = 2, NSPE 
              MONSPE = MONSPE+1 
              K5COL(MONSPE) = 'U' 
              CSMSK5(MONSPE) = trim(adjustl(SPNAME(K)))
            ENDDO 
         ENDIF 
!........LOAD APPROPRIATE TEXT INTO TITLE FOR ".nod" FILE(fNOD)          
         IF(CH_K5COL(M) .EQ.'N') CSMSK5(MONSPE) = VARNK5(1) 
         IF(CH_K5COL(M) .EQ.'X') CSMSK5(MONSPE) = VARNK5(2) 
         IF(CH_K5COL(M) .EQ.'Y') CSMSK5(MONSPE) = VARNK5(3) 
         IF(CH_K5COL(M) .EQ.'Z') CSMSK5(MONSPE) = VARNK5(4) 
         IF(CH_K5COL(M) .EQ.'P') CSMSK5(MONSPE) = VARNK5(5) 
         IF(CH_K5COL(M) .EQ.'S') CSMSK5(MONSPE) = VARNK5(7) 
         IF(CH_K5COL(M) .EQ.'R') CSMSK5(MONSPE) = VARNK5(8) 
!........VERSION 1.1
         IF(CH_K5COL(M) .EQ.'H') CSMSK5(MONSPE) = VARNK5(9) 
         IF(CH_K5COL(M) .EQ.'V') CSMSK5(MONSPE) = VARNK5(10) 
!         IF(CH_K5COL(M) .EQ.'K') CSMSK5(MONSPE) = VARNK5(11) 		 !HYDCON 01/6/2017
         !TEST IF UNSATURATED SOLUTION...IF NOT ELIMINATE SW FROM OUTPUT
!         IF(CH_K5COL(M) .EQ.'S' .AND. .not.IUNSAT) MONSPE=MONSPE-1
         IF(CH_K5COL(M) .EQ.'S' .AND. IUNSAT.NE.1) MONSPE=MONSPE-1		!MT: revised due to error during compiling with gfortran
!........VERSION 1.1
!........TEST FOR PRESENCE OF P, U, SW, RHO, OR TOTAL HEAD FOR TECPLOT                     
         IF(CH_K5COL(M) .EQ.'P')                   LTpP                     = .TRUE. 
         IF(CH_K5COL(M) .EQ.'U')                   LTpU                     = .TRUE. 
         !ONLY OUTPUT SW TO TECPLOT IF UNSATURATED SOLUTION (IUNSAT=+1)
         IF(CH_K5COL(M) .EQ.'S' .AND. IUNSAT.EQ.1) LTpSW                    = .TRUE. 
         IF(CH_K5COL(M) .EQ.'R')                   LTpRho                   = .TRUE.
         IF(CH_K5COL(M) .EQ.'H')                   LTpTotalHead             = .TRUE.
         IF(CH_K5COL(M) .EQ.'V')                   LTpVelocity              = .TRUE.
         IF(CH_K5COL(M) .EQ.'V')                   LTpPrintVelocityAtNodes  = .TRUE.
!         IF(CH_K5COL(M) .EQ.'K')                   LTpHydcon                = .TRUE.		 
!........TEST FOR DATA TERMINATION CHARACTER                            
         IF(CH_K5COL(M) .EQ.'-') THEN 
            NCOLS5 = MONSPE-1 
            GOTO 1162 
         ENDIF 
                                                                        
 1160 END DO 
      NCOLS5 = NCOLMX 
 1162 CONTINUE 
      IF (fNOD.GT.0) WRITE (fLST, 1172) NCOLPR 
 1172 FORMAT (//11X,'.nod FILE'/11X,'---------'                         &
     &   //11X,I6,5X,'PRINTED OUTPUT CYCLE (IN TIME STEPS)'/)           
!.....VERSION 1.1                                                                       
!.....TECPLOT NODAL DATA OUTPUT                                         
      IF (fTPN.GT.0) WRITE (fLST, 1173) NTPNP 
 1173 FORMAT (//11X,'TECPLOT NODAL DATA OUTPUT FILE'/11X,'---------'    &
     &   //11X,I6,5X,'PRINTED OUTPUT CYCLE (IN TIME STEPS)'/)           
!                                                                       
!.....SPECIES COUNTER                                                   
      KK = 0 
      DO 1250 M = 1, NCOLS5 
         DO 1200 MM = 1, (K5MAX + NSPE-1) 
            IF(K5COL(M) .EQ.CSYM(MM) ) THEN 
               IF((MM.EQ.1) .AND.(M.NE.1) ) THEN 
                  call ErrorIO('Error NODE NUMBER CAN ONLY APPEAR IN COLUMN 1')
               ENDIF 
               IF((MM.EQ.4) .AND.(IABS(KTYPE) .EQ.2) ) THEN 
                  call ErrorIO('Error CANNOT PRINT "Z" FOR A 2-D PROBLEM')
               ENDIF 
               J5COL(M) = MM 
!..............FOR MULTISPECIES SIMULATIONS                             
!              RESET CSYM AFTER IT IS USED                              
               IF(K5COL(M) .EQ.'U') THEN 
                  KK = KK + 1 
                  WRITE(CFIND, '(I1)') KK 
!................RESET CSYM TO SPECIES NUMBER TO                        
!                CORRECTLY ASSIGN COLUMNS                               
                  CSYM(MM) = CFIND 
               ENDIF 
               GOTO 1250 
            ENDIF 
 1200    END DO 
         call ErrorIO('Error K5COL('//trim(Val2Char(M))//'-- UNRECOGNIZED VALUE')
 1250 END DO 
      WRITE (fLST, 1300) (M, CSMSK5 (M), M = 1, NCOLS5) 
 1300 FORMAT (11X,'COLUMN ',I1,':',2X,A) 
!.....VERSION 1.1
!.....RESET NCOLS5 IF OUTPUTTING RHO, TOTAL HEAD, AND/OR VELOCITY
!     SO THEY ARE NOT INCLUDED IN THE STANDARD ASCII NODE FILE
      IF(LTpRho)       NCOLS5=NCOLS5-1
      IF(LTpTotalHead) NCOLS5=NCOLS5-1
      IF(LTpVelocity)  NCOLS5=NCOLS5-1
!      IF(LTpHydcon)  NCOLS5=NCOLS5-1	  		!HYDCON 01/6/2017
!                                                                       
!.....INPUT DATASET 8C: OUTPUT CONTROLS AND OPTIONS -- UNIT fELE          
      MSErrorValue%cDataSet=' 8C'
      CALL SKPCOM (fINP, NLSKIP) 
      READ (fINP,*,iostat=ios) LCOLPR 
      if(ios/=0) call ErrorIO('Error specifying element output file controls (LCOLPR)')
!.....VERSION 1.1
!.....SET TECPLOT ELEMENT DATA OUTPUT TO HAVE THE SAME FREQUENCY        
!     A STANDARD SUTRA ELEMENT DATA OUTPUT.                             
!     TO OUTPUT TO A TECPLOT FILE AND NOT TO THE STANDARD SUTRA         
!     ELEMENT DATA OUTPUT FILE SET LCOLPR TO THE DESIRED FREQUENCY      
!     AND DO NOT SPECIFY A fELE FILE IN SUTRA.FIL                         
      IF (fTPE.GT.0) NTPEP = LCOLPR 
                                                                       
      DO 2160 M = 1, NCOLMX 
         BACKSPACE (fINP) 
         READ (fINP,*,iostat=ios) LCOLPR, (K6COL (MM), MM = 1, M) 
         if(ios/=0) call ErrorIO('Error specifying element output elements')
         IF (K6COL (M) .EQ.'-') THEN 
            NCOLS6 = M - 1 
            GOTO 2162 
         ENDIF 
 2160 END DO 
      NCOLS6 = NCOLMX 
 2162 CONTINUE 
      IF (fELE.GT.0) WRITE (fLST, 2172) LCOLPR 
 2172 FORMAT (//11X,'.ele FILE'/11X,'---------'                         &
     &   //11X,I6,5X,'PRINTED OUTPUT CYCLE (IN TIME STEPS)'/)           
!.....VERSION 1.1                                                                       
!.....TECPLOT ELEMENT DATA OUTPUT                                       
      IF (fTPE.GT.0) WRITE (fLST, 2173) NTPEP 
 2173 FORMAT (//11X,'TECPLOT ELEMENT DATA OUTPUT FILE'/11X,'---------'  &
     &   //11X,I6,5X,'PRINTED OUTPUT CYCLE (IN TIME STEPS)'/)           
                                                                       
      DO 2250 M = 1, NCOLS6 
         DO 2200 MM = 1, 7 
            IF (K6COL (M) .EQ.K6SYM (MM) ) THEN 
               IF ( (MM.EQ.1) .AND. (M.NE.1) ) THEN 
                  call ErrorIO('Error ELEMENT NUMBER CAN ONLY APPEAR IN COLUMN 1')
               ENDIF 
               IF ( (MM.EQ.4) .AND. (IABS (KTYPE) .EQ.2) ) THEN 
                  call ErrorIO('Error CANNOT PRINT "Z" FOR A 2-D')
               ENDIF 
               IF ( (MM.EQ.7) .AND. (IABS (KTYPE) .EQ.2) ) THEN 
                  call ErrorIO('Error CANNOT PRINT "VZ" FOR A 2-D PROBLEM')
               ENDIF 
               J6COL (M) = MM 
               GOTO 2250 
            ENDIF 
 2200    END DO 
         call ErrorIO('Error K6COL('//trim(Val2Char(M))//'-- UNRECOGNIZED VALUE')
 2250 END DO 
      WRITE (fLST, 2300) (M, VARNK6 (J6COL (M) ), M = 1, NCOLS6) 
 2300 FORMAT (11X,'COLUMN ',I1,':',2X,A) 
!                                                                       
!.....INPUT DATASET 8D: OUTPUT CONTROLS AND OPTIONS -- UNIT fOBS
      if(fSOB>0) goto 4399          
      MSErrorValue%cDataSet=' 8D'
      NOBCYC = ITMAX + 1 
      IF (NOBSN - 1.EQ.0) GOTO 4399 
      CALL SKPCOM (fINP, NLSKIP) 
!.....NOBS IS ACTUAL NUMBER OF OBSERVATION NODES                        
!     NTOBS IS MAXIMUM NUMBER OF TIME STEPS WITH OBSERVATIONS           
      NOBS = NOBSN - 1 
      NTOBS = NTOBSN - 2 
      if (fSOB==0) then
        READ (fINP,*,iostat=ios) NOBCYC, (IOBS (JJ), JJ = 1, NOBSN) 
      end if
      if(ios/=0) call ErrorIO('Error specifying observation point output elements')
      if (fSOB/=0) goto 4399
      WRITE (fLST, 3172) NOBCYC 
 3172 FORMAT (//11X,'.obs FILE'/11X,'---------'                         &
     &   //11X,I6,5X,'PRINTED OUTPUT CYCLE (IN TIME STEPS)'/)           
      ITCNT = 0 
      JSTOP = 0 
      WRITE (fLST, 4060) 
 4060 FORMAT(////11X,'O B S E R V A T I O N   N O D E S') 
      WRITE (fLST, 4065) NOBCYC, NTOBS 
 4065 FORMAT(//11X,'OBSERVATIONS WILL BE MADE EVERY ',I5,' TIME STEPS,' &
     &   /11X,'AS WELL AS ON THE FIRST AND LAST TIME STEP,'             &
     &   /11X,'FOR A TOTAL OF ',I5,' TIME STEPS.')                      
      WRITE (fLST, 4070) 
 4070 FORMAT(//11X,'**** NODES AT WHICH OBSERVATIONS WILL BE MADE',     &
     &   ' ****'//)                                                     
      IF (IOBS (NOBSN) .NE.0) JSTOP = 1 
      DO 4200 JJ = 1, NOBS 
         IF (IOBS (JJ) .LE.0) JSTOP = 1 
 4200 END DO 
      WRITE (fLST, 4300) (IOBS (JJ), JJ = 1, NOBS) 
 4300 FORMAT((11X,10(1X,I9))) 
!.....END SIMULATION IF CORRECTIONS ARE NECESSARY OBSERVATION NODE LIST 
      IF (JSTOP.EQ.1) THEN 
         call ErrorIO('Error in observation node list')
      ENDIF 
 4399 CONTINUE 
!                                                                       
!.....INPUT DATASET 9A: BASE FLUID PROPERTIES                           
      MSErrorValue%cDataSet=' 9A'
      CALL SKPCOM (fINP, NLSKIP) 
      READ (fINP,*,iostat=ios) COMPFL, CW, BSIGMAW, RHOW0, BURHOW0, BDRWDU, BVISC0 
      if(ios/=0) call ErrorIO('Error specifying base fluid properties')
!.....SET SIGMAW, URHOW0, DRWDU, AND VISC0 TO BASE VALUES               
      SIGMAW = BSIGMAW 
      URHOW0 = BURHOW0 
      DRWDU = BDRWDU 
      BVISC0 = BVISC0 
!.....IF ENERGY TRANSPORT ONLY OR MULTISPECIES INCLUDING ENERGY TRANSPORT
      IF (ME.GE.0) THEN 													
         VISC0 (NESP) = 1                                                   
         BVISC0 = 0.0D0                                                     
      ENDIF                                                                 
!      
!.....INPUT DATASET 9B: MULTI SPECIES FLUID PROPERTIES                  
!.....READ ADDITIONAL SIGMAW, URHOW0, DRWDU, AND VISC0 FOR              
!     MULTIPLE SPECIES                                                  
      MSErrorValue%cDataSet=' 9B'
      IF (NSPE.GT.1) THEN 
         CALL SKPCOM (fINP, NLSKIP) 
!.......READ SPECIES NUMBER                                             
 2000    READ (fINP,*,iostat=ios) K 
         if(ios/=0) call ErrorIO('Error specifying species number for multi-species fluid properties')
!.......CONTINUE READING UNTIL K EQUALS ZERO IF IT EXCEEDS NSPE         
         IF (K.GT.NSPE) GOTO 2000 
!.......TERMINATE READ IF K EQUALS ZERO                                 
         IF (K.LT.1) GOTO 2010 
!.......BACKSPACE AND RE-READ LINE                                      
         BACKSPACE (fINP) 
         READ (fINP,*,iostat=ios) K, SIGMAW (K), URHOW0 (K), DRWDU (K), VISC0 (K) 
         if(ios/=0) call ErrorIO('Error specifying multi-species fluid properties')
!.......CONTINUE READING DATASET 9B UNTIL K EQUALS ZERO                 
         GOTO 2000 
      ENDIF 
                                                                    
!.....INPUT DATASET 10: SOLID MATRIX PROPERTIES                         
 2010 MSErrorValue%cDataSet=' 10'
      CALL SKPCOM (fINP, NLSKIP) 
      READ (fINP,*,iostat=ios) COMPMA, CS, SIGMAS, RHOS 
      if(ios/=0) call ErrorIO('Error specifying solid matrix properties')
      !initialize NodeData()%rhos and NodeData()%compma
      !to default value if zone file is not being used - v2009.1 bug fix
      if( fZON == 0 ) then
        NodeData(:)%rhos   = RHOS
        NodeData(:)%compma = COMPMA
      end if
!.....ENERGY TRANSPORT                                                  
      IF (ME.EQ. + 1) WRITE (fLST, 210) COMPFL, COMPMA, CW, CS, VISC0(NESP), &
      RHOS, RHOW0, DRWDU, URHOW0, SIGMAW, SIGMAS                        
  210 FORMAT(1H1////11X,'C O N S T A N T   P R O P E R T I E S   O F',  &
     &   '   F L U I D   A N D   S O L I D   M A T R I X'               &
     &   //11X,1PD15.4,5X,'COMPRESSIBILITY OF FLUID'/11X,1PD15.4,5X,    &
     &   'COMPRESSIBILITY OF POROUS MATRIX'//11X,1PD15.4,5X,            &
     &   'SPECIFIC HEAT CAPACITY OF FLUID',/11X,1PD15.4,5X,             &
     &   'SPECIFIC HEAT CAPACITY OF SOLID GRAIN'//13X,'FLUID VISCOSITY',&
     &   ' IS CALCULATED BY SUTRA AS A FUNCTION OF TEMPERATURE IN ',    &
     &   'UNITS OF {kg/(m*s)}'//11X,1PD15.4,5X,'VISC0, CONVERSION ',    &
     &   'FACTOR FOR VISCOSITY UNITS,  {desired units} = VISC0*',       &
     &   '{kg/(m*s)}'//11X,1PD15.4,5X,'DENSITY OF A SOLID GRAIN'        &
     &   //13X,'FLUID DENSITY, RHOW'/13X,'CALCULATED BY ',              &
     &   'SUTRA IN TERMS OF TEMPERATURE, U, AS:'/13X,'RHOW = RHOW0 + ', &
     &   'DRWDU*(U-URHOW0)'//11X,1PD15.4,5X,'FLUID BASE DENSITY, RHOW0' &
     &   /11X,1PD15.4,5X,'COEFFICIENT OF DENSITY CHANGE WITH ',         &
     &   'TEMPERATURE, DRWDU'/11X,1PD15.4,5X,'TEMPERATURE, URHOW0, ',   &
     &   'AT WHICH FLUID DENSITY IS AT BASE VALUE, RHOW0'               &
     &   //11X,1PD15.4,5X,'THERMAL CONDUCTIVITY OF FLUID'               &
     &   /11X,1PD15.4,5X,'THERMAL CONDUCTIVITY OF SOLID GRAIN, SIGMAS'  &
     &   /31X,'A NEGATIVE VALUE FOR SIGMAS INDICATES SLAMBDA IS ',      &
     &   'IS DEFINED FOR EACH ELEMENT IN DATASET 15D'/31X,              &
     &   'A POSITIVE VALUE FOR SIGMAS INDICATES SLAMBDA IS CONSTANT')   
!.....SOLUTE TRANSPORT OR SOLUTE AND ENERGY TRANSPORT                   
!.....NESP=0 IF ME=-1 {SOLUTE TRANSPORT}                                
!.....NESP=# OF ENERGY TRANSPORT SPECIES FOR ME=0                       
!.....                {SOLUTE & ENERY TRANSPORT}                        
      IF (ME.LT. + 1) THEN 
!.......TEST IF SOLUTE ONLY SIMULATION                                  
!       IF SO MAKE SURE SIGMAS IS NON-NEGATIVE                          
         IF (NESP.EQ.0) SIGMAS = DABS (SIGMAS) 
!.......WRITE DATA TO OUTPUT FILE                                       
         WRITE (fLST, 220) COMPFL, COMPMA, RHOS, RHOW0 
         IF (NESP.LT.1.AND.NSPE.GT.1) WRITE (fLST, 221) BVISC0 
         DO K = 1, NSPE 
         WRITE (fLST, 222) K, trim(adjustl(SPNAME(K)))
!.........WRITE ENERGY TRANSPORT PARAMETERS                             
         IF (K.EQ.NESP) THEN 
            WRITE (fLST, 223) CW, CS, VISC0 (K), DRWDU (K), URHOW0 (K),   &
            SIGMAW (K), SIGMAS                                          
         ELSE 
            WRITE (fLST, 224) DRWDU (K), URHOW0 (K), SIGMAW (K), VISC0 (K) 
         ENDIF 
         ENDDO 
      ENDIF 
  220 FORMAT(1H1////11X,'C O N S T A N T   P R O P E R T I E S   O F',  &
     &   '   F L U I D   A N D   S O L I D   M A T R I X'               &
     &   //11X,1PD15.4,5X,'COMPRESSIBILITY OF FLUID, COMPFL',           &
     &   /11X,1PD15.4,5X,'COMPRESSIBILITY OF POROUS MATRIX, COMPMA '    &
     &   /11X,1PD15.4,5X,'DENSITY OF A SOLID GRAIN, RHOS'               &
     &   //13X,'FLUID DENSITY, RHOW'/13X,'CALCULATED BY ',              &
     &   'SUTRA IN TERMS OF SOLUTE CONCENTRATION/TEMPERTURE, U, AS:',   &
     &   /13X,'RHOW = RHOW0 + DRWDU*(U-URHOW0)'                         &
     &   //11X,1PD15.4,5X,'FLUID BASE DENSITY, RHOW0')                  
!                                                                       
  221 FORMAT(//11X,1PD15.4,5X,'FLUID VISCOSITY',/) 
!                                                                       
  222 FORMAT(//,17X,'S P E C I E S [',I3,']'2X,A,//) 
!                                                                       
  223 FORMAT(11X,1PD15.4,5X,                                            &
     &   'SPECIFIC HEAT CAPACITY OF FLUID, CW',/11X,1PD15.4,5X,         &
     &   'SPECIFIC HEAT CAPACITY OF SOLID GRAIN, CS'//13X,              &
     &   'FLUID VISCOSITY IS CALCULATED BY SUTRA AS A FUNCTION OF ',    &
     &   'TEMPERATURE IN UNITS OF {kg/(m*s)}'//11X,1PD15.4,5X,          &
     &   'VISC0, CONVERSION FACTOR FOR VISCOSITY UNITS,  ',             &
     &   '{desired units} = VISC0*{kg/(m*s)}',                          &
     &   /11X,1PD15.4,5X,'COEFFICIENT OF DENSITY CHANGE WITH ',         &
     &   'TEMPERATURE, DRWDU'/11X,1PD15.4,5X,'TEMPERATURE, URHOW0, ',   &
     &   'AT WHICH FLUID DENSITY IS AT BASE VALUE, RHOW0',              &
     &   /11X,1PD15.4,5X,'THERMAL CONDUCTIVITY OF FLUID, SIGMAW',       &
     &   /11X,1PD15.4,5X,'THERMAL CONDUCTIVITY OF SOLID GRAIN, SIGMAS'  &
     &   /31X,'A NEGATIVE VALUE FOR SIGMAS INDICATES SLAMBDA IS ',      &
     &   'IS DEFINED FOR EACH ELEMENT IN DATASET 15D'/31X,              &
     &   'A POSITIVE VALUE FOR SIGMAS INDICATES SLAMBDA IS CONSTANT')   
!                                                                       
  224 FORMAT(11X,1PD15.4,5X,'COEFFICIENT OF DENSITY CHANGE WITH ',      &
     &   'SOLUTE CONCENTRATION, DRWDU'                                  &
     &   /11X,1PD15.4,5X,'SOLUTE CONCENTRATION, URHOW0, ',              &
     &   'AT WHICH FLUID DENSITY IS AT BASE VALUE, RHOW0'               &
     &   /11X,1PD15.4,5X,'MOLECULAR DIFFUSIVITY OF SOLUTE IN FLUID, ',  &
     &   'SIGMAW'/11X,1PD15.4,5X,'COEFFICIENT OF VISCOSITY CHANGE '     &
     &   'SOLUTE CONCENTRATION, VISC0')                                 
!
!.....Calculate summary of viscosity equation
      cViscosityEquation='Viscosity = Base Viscosity'
      cViscosity='Viscosity = '//trim(adjustl(Val2Char(BVISC0)))
      do k=1,NSPE
        if(k==NESP) then
          cViscosityEquation=trim(adjustl(cViscosityEquation))//' + (Viscosity Multiplier * SUTRA Temperature Function)'
          cViscosity=trim(adjustl(cViscosity))//' + ('//trim(adjustl(Val2Char(VISC0(k))))//' * SUTRA Temperature Function)'
        else
          cViscosityEquation=trim(adjustl(cViscosityEquation))//' &
+ (dViscosity/d'//trim(adjustl(SPNAME(K)))//' * '//trim(adjustl(SPNAME(K)))//')'
          cViscosity=trim(adjustl(cViscosity))//' + ('//trim(adjustl(Val2Char(VISC0(k))))//' * '//trim(adjustl(SPNAME(K)))//')'
        end if
      end do
      write (fLST,'(///132("-")/,a,/,a,/,a,/,132("-")/)') &
        'SUTRA-MS Viscosity Equation', cViscosityEquation, cViscosity

!
!.....Set SOP from zones
      if(.not.MkZoneSOP(COMPFL)) call ErrorIO('Error specifying SOP for nodal zones')
!                                                                       
!.....INITIALIZE ADSMOD TO NONE AND ADSORPTION COEFFICIENTS TO ZERO     
      ADSMOD = 'NONE      ' 
      CHI1 = 0.0D0 
      CHI2 = 0.0D0 
!.....INPUT DATASET 11: ADSORPTION PARAMETERS                           
      MSErrorValue%cDataSet=' 11'
      CALL SKPCOM (fINP, NLSKIP) 
      IF (NSPE.LT.2) THEN 
         READ (fINP,*,iostat=ios) ADSMOD(1) 
         if(ios/=0) call ErrorIO('Error specifying ADSMOD for species 1')
!.......NON-SORBING SOLUTE OR ENERGY                                    
      IF (ADSMOD (1) .NE.'NONE      '.AND.NESP.NE.1) THEN 
            BACKSPACE (fINP) 
            READ (fINP,*,iostat=ios) ADSMOD(1), CHI1(1), CHI2(1) 
            if(ios/=0) call ErrorIO('Error specifying adsorption parameters for species 1')
         ENDIF 
      ELSE 
 2015    READ (fINP,*,iostat=ios) K 
         if(ios/=0) call ErrorIO('Error specifying adsorption species number')
!.......CONTINUE READING UNTIL K EQUALS ZERO IF IT EXCEEDS NSPE         
         IF (K.GT.NSPE) GOTO 2015 
!.......TERMINATE READ IF K EQUALS ZERO                                 
         IF (K.LT.1) GOTO 2020 
!.......BACKSPACE AND RE-READ LINE                                      
         BACKSPACE (fINP) 
         READ (fINP,*,iostat=ios) K, ADSMOD(K) 
         if(ios/=0) call ErrorIO('Error specifying ADSMOD for species '//trim(Val2Char(K)))
!.......NON-SORBING SOLUTE OR ENERGY                                    
         IF (ADSMOD(K).NE.'NONE      '.AND.K.NE.NESP) THEN 
            BACKSPACE (fINP) 
            READ (fINP,*,iostat=ios) K, ADSMOD (K), CHI1 (K), CHI2 (K) 
            if(ios/=0) call ErrorIO('Error specifying sorption parameters for species '//trim(Val2Char(K)))
         ENDIF 
!.......CONTINUE READING DATASET 9B UNTIL K EQUALS ZERO                 
         GOTO 2015 
      ENDIF 
!.....ENERGY TRANSPORT - NO SORPTION WITH ENERGY ONLY SIMULATION        
 2020 IF (ME.EQ. + 1) GOTO 249 
      DO 248 K = 1, NSPE 
!.....CHECK IF NON-SORBING SOLUTE OR ENERGY                             
        IF (ADSMOD (K) .EQ.'NONE      '.OR.K.EQ.NESP) GOTO 234 
        WRITE (fLST, 232) ADSMOD (K), K, trim(adjustl(SPNAME(K)))
  232   FORMAT(////11X,'A D S O R P T I O N   P A R A M E T E R S'        &
     &           //17X,A10,5X,'EQUILIBRIUM SORPTION ISOTHERM',            &
     &           //16X,' S P E C I E S [',I3,']'/17X,A)                       
        GOTO 236 
!.....NON-SORBING SOLUTE OR ENERGY                                      
  234   WRITE (fLST, 235) K, trim(adjustl(SPNAME(K)))
  235   FORMAT(////11X,'A D S O R P T I O N   P A R A M E T E R S'        &
     &           //16X,' S P E C I E S [',I3,']'/17X,A,                   &
     &            /17X,'NON-SORBING SOLUTE')                                    
        GOTO 248 
!.....CHECK FOR CORRECT KEYWORDS                                        
  236   IF((trim(ADSMOD(K)) .EQ.'NONE') .OR.(trim(ADSMOD(K)) .EQ.'LINEAR') .OR. &
     &     (trim(ADSMOD(K)) .EQ.'FREUNDLICH') .OR.(trim(ADSMOD(K)) .EQ.'LANGMUIR') ) GOTO 238                                                     
        !Correct adsorption keywords not found
        call ErrorIO('Error: unknown sorption type ['//trim(ADSMOD(K))//'specifed for species '//trim(Val2Char(K))//' &
['//trim(adjustl(SPNAME(K)))//']')

  238   IF (trim(ADSMOD(K)) .EQ.'LINEAR') WRITE (fLST, 242) CHI1 (K) 
  242   FORMAT(11X,1PD15.4,5X,'LINEAR DISTRIBUTION COEFFICIENT') 

        IF (trim(ADSMOD(K)) .EQ.'FREUNDLICH') WRITE (fLST, 244) CHI1 (K) ,CHI2 (K)                                                       
  244   FORMAT(11X,1PD15.4,5X,'FREUNDLICH DISTRIBUTION COEFFICIENT'       &
     &        /11X,1PD15.4,5X,'SECOND FREUNDLICH COEFFICIENT')               

        IF (trim(ADSMOD(K)) .EQ.'FREUNDLICH'.AND.CHI2 (K) .LE.0.D0) THEN 
          call ErrorIO('Error: second sorption coefficient for species '//trim(Val2Char(K))//' &
['//trim(adjustl(SPNAME(K)))//'] must be greater than 0.0')
        ENDIF 
      IF (trim(ADSMOD(K)) .EQ.'LANGMUIR') WRITE (fLST, 246) CHI1(K) , CHI2(K)                                                                
  246 FORMAT(11X,1PD15.4,5X,'LANGMUIR DISTRIBUTION COEFFICIENT'         &
     &   /11X,1PD15.4,5X,'SECOND LANGMUIR COEFFICIENT')                 

      IF ( (trim(ADSMOD(K)) .EQ.'LANGMUIR' .OR. trim(ADSMOD(K)) .EQ.'FREUNDLICH') .AND. ITRMAX.EQ.1) &
        WRITE (fLST, 247) K, SPNAME (K)                      
  247 FORMAT(/11X,'W A R N I N G - SORPTION MODEL IS NON-LINEAR, ',     &
     &   'BUT ITERATIVE SOLUTION WAS NOT REQUESTED.'/11X,'PLEASE CHECK',&
     &   ' ITERATION CONTROLS!',                                        &
     &   /16X,' S P E C I E S [',I3,'] - ',A)                                

!.....END OF NSPE LOOP FOR PRINTING ADSORPTION PARAMETERS               
  248 END DO 
!                                                                       
!.....INITIALIZE PRODUCTION OF ENERGY AND SOLUTE MASS VARIABLES         
      PRODF0 = 0.0D0 
      PRODS0 = 0.0D0 
      PRODF1 = 0.0D0 
      PRODS1 = 0.0D0 
!.....INPUT DATASET 12: PRODUCTION OF ENERGY OR SOLUTE MASS             
  249 MSErrorValue%cDataSet=' 12'
      CALL SKPCOM (fINP, NLSKIP) 
  250 IF (NSPE.LT.2) THEN 
         K = 1 
         READ(fINP,*,iostat=ios) PRODF0(K), PRODS0(K), PRODF1(K), PRODS1(K) 
         if(ios/=0) call ErrorIO('Error specifying energy or solute mass production variables for species 1')
      ELSE 
         READ(fINP,*,iostat=ios) K 
         if(ios/=0) call ErrorIO('Error specifying species number for production variables')
!.......IF K EXCEEDS NSPE CONTINUE READING INPUT DATASET 12             
         IF(K.GT.NSPE) GOTO 250 
!.......TERMINATE READING DATASET 12 IF K EQUALS ZERO                   
         IF(K.EQ.0) GOTO 266 
!.......BACKSPACE UNIT fINP AND READ THE REST OF DATA FOR                 
!       FOR SPECIES K                                                   
         BACKSPACE(fINP) 
         READ(fINP,*,iostat=ios) K, PRODF0(K), PRODS0(K), PRODF1(K), PRODS1(K) 
         if(ios/=0) call ErrorIO('Error specifying energy or solute mass production variables for species '//trim(Val2Char(K)))
      ENDIF
!.....WRITE DATASET 12 TO OUTPUT FILE (UNIT fLST)                         
!.....SOLUTE TRANSPORT                                                  
      IF(K.NE.NESP) WRITE(fLST, 251) K, trim(adjustl(SPNAME(K))), PRODF0(K), PRODS0(K), PRODF1(K), PRODS1(K)                                        
  251 FORMAT(////11X,'P R O D U C T I O N   A N D   D E C A Y   O F   ',&
     &   'S P E C I E S   M A S S'/12X,' S P E C I E S [',I3,'] - ',    &
     &   A,//13X,'PRODUCTION RATE (+)'/13X,'DECAY RATE (-)',          &
     &   //11X,1PD15.4,5X,'ZERO-ORDER RATE OF SOLUTE '                  &
     &   'MASS PRODUCTION/DECAY IN FLUID'/11X,1PD15.4,5X,               &
     &   'ZERO-ORDER RATE OF ADSORBATE MASS PRODUCTION/DECAY IN ',      &
     &   'IMMOBILE PHASE'/11X,1PD15.4,5X,'FIRST-ORDER RATE OF SOLUTE ', &
     &   'MASS PRODUCTION/DECAY IN FLUID'/11X,1PD15.4,5X,               &
     &   'FIRST-ORDER RATE OF ADSORBATE MASS PRODUCTION/DECAY IN ',     &
     &   'IMMOBILE PHASE')                                              
!.....ENERGY TRANSPORT                                                  
      IF(K.EQ.NESP) WRITE(fLST, 255) K, trim(adjustl(SPNAME(K))), PRODF0(K), PRODS0(K)                                                                
  255 FORMAT(////11X,'P R O D U C T I O N   A N D   L O S S   O F   ',  &
     &   'E N E R G Y'/12X,' S P E C I E S [',I3,'] - ',                &
     &   A,//13X,'PRODUCTION RATE (+)'/13X,'LOSS RATE (-)',           &
     &   //11X,1PD15.4,5X,'ZERO-ORDER RATE OF ENERGY ',                 &
     &   'PRODUCTION/LOSS IN FLUID'/11X,1PD15.4,5X,                     &
     &   'ZERO-ORDER RATE OF ENERGY PRODUCTION/LOSS IN ',               &
     &   'SOLID GRAINS')                                                
!                                                                       
!.....TERMINATE READING OF DATASET 12 IF THIS IS A SINGLE               
!     SPECIES TRANSPORT SIMULATION                                      
      IF (NSPE.LT.2) GOTO 266 
!                                                                       
  265 GOTO 250 
  266 CONTINUE 
!                                                                       
!.....SET PARAMETER SWITCHES FOR EITHER SOLUTE, ENERGY/SOLUTE, OR       
!.....SOLUTE TRANSPORT
      select case (ME)
        case(-1)
          goto 272
        case(0)
          goto 273
        case(1)
          goto 274
        case default
          call ErrorIO('Error:  ME must be -1, 0, or 1 ['//trim(Val2Char(ME))//']')
      end select
!     FOR SOLUTE ONLY TRANSPORT:                                        
  272 CS = 0.0D0 
      CW = 1.D00 
      SIGMAS = 0.0D0 
      GOTO 278 
!.....FOR SOLUTE AND ENERGY TRANSPORT                                   
  273 ADSMOD(NESP)  = 'NONE      ' 
      CHI1(NESP)    = 0.0D0 
      CHI2(NESP)    = 0.0D0 
      PRODF1(NESP)  = 0.0D0 
      PRODS1(NESP)  = 0.0D0 
      GOTO 278 
!     FOR ENERGY TRANSPORT:                                             
  274 ADSMOD(1)  = 'NONE      ' 
      CHI1(1)    = 0.0D0 
      CHI2(1)    = 0.0D0 
      PRODF1(1)  = 0.0D0 
      PRODS1(1)  = 0.0D0 
  278 CONTINUE 
      !initialize production and sorption variables of ProdSorp
      !to default values if zone file is not being used - VERSION 1.1
      if( fZON == 0 ) then
        do I = 1, NN
          do K = 1, NSPE
            ProdSorp(I)%prodf0(K) = PRODF0(K)
            ProdSorp(I)%prods0(K) = PRODS0(K)
            ProdSorp(I)%prodf1(K) = PRODF1(K)
            ProdSorp(I)%prods1(K) = PRODS1(K)
            ProdSorp(I)%chi1(K)   = CHI1(K)
            ProdSorp(I)%chi2(K)   = CHI2(K)
          end do
        end do
      end if
!                                                                       
      IF (IABS (KTYPE) .EQ.3) THEN 
!.....READ 3-D INPUT FROM DATASETS 13 - 15.                             
!                                                                       
!.....INPUT DATASET 13: ORIENTATION OF COORDINATES TO GRAVITY           
         MSErrorValue%cDataSet=' 13'
         CALL SKPCOM (fINP, NLSKIP) 
         READ (fINP,*,iostat=ios) GRAVX, GRAVY, GRAVZ 
         if(ios/=0) call ErrorIO('Error specifying gravity for 3D problem')
         WRITE (fLST, 320) GRAVX, GRAVY, GRAVZ 
  320 FORMAT(////11X,'C O O R D I N A T E   O R I E N T A T I O N   ',  &
     &   'T O   G R A V I T Y'//13X,'COMPONENT OF GRAVITY VECTOR',      &
     &   /13X,'IN +X DIRECTION, GRAVX'/11X,1PD15.4,5X,                  &
     &   'GRAVX = -GRAV * D(ELEVATION)/DX'//13X,'COMPONENT OF GRAVITY', &
     &   ' VECTOR'/13X,'IN +Y DIRECTION, GRAVY'/11X,1PD15.4,5X,         &
     &   'GRAVY = -GRAV * D(ELEVATION)/DY'//13X,'COMPONENT OF GRAVITY', &
     &   ' VECTOR'/13X,'IN +Z DIRECTION, GRAVZ'/11X,1PD15.4,5X,         &
     &   'GRAVZ = -GRAV * D(ELEVATION)/DZ')                             
!                                                                       
!.....INPUT DATASETS 14A AND 14B: NODEWISE DATA                         
         MSErrorValue%cDataSet='14A'
         CALL SKPCOM (fINP, NLSKIP) 
         READ (fINP,*,iostat=ios) CDUM10, SCALX, SCALY, SCALZ, PORFAC 
         if(ios/=0) call ErrorIO('Error specifying nodewise multiplier data')
         IF (trim(CDUM10).NE.'NODE') THEN 
           call ErrorIO('Error: Dataset 14A must begin with the word "NODE"')
         ENDIF 

         NRTEST = 1 
         MSErrorValue%cDataSet='14B'
         CALL SKPCOM (fINP, NLSKIP) 
         DO I = 1, NN 
            READ(fINP,*,iostat=ios) II
            if(ios/=0) call ErrorIO('Error specifying node number - error on line ['//trim(Val2Char(I))//'] of Data Set 14B')
            if(I>NN)  call ErrorIO('Error: Specified node number ['//trim(Val2Char(II))//'] exceeds NE ['//trim(Val2Char(NN))//']')
            backspace(fINP)
            if(fZON==0) then
              READ(fINP,*,iostat=ios) II, NodeMap(ii), X(II), Y(II), Z(II), NodeData(II)%por
              NodeMap(ii) = II
            end if
            if(fZON/=0) READ(fINP,*,iostat=ios) II, NodeMap(ii), X(II), Y(II), Z(II)
            if(ios/=0) call ErrorIO('Error specifying nodewise data - error on line ['//trim(Val2Char(I))//'] of Data Set 14B')
            X(II) = X(II) * SCALX 
            Y(II) = Y(II) * SCALY 
            Z(II) = Z(II) * SCALZ 
            if(fZON==0) NodeData(NodeMap(ii))%por = NodeData(NodeMap(ii))%por * PORFAC 
!           SET SPECIFIC PRESSURE STORATIVITY, SOP.                           
            if(fZON==0) NodeData(NodeMap(ii))%sop = (1.D0 - NodeData(NodeMap(ii))%por ) &
* NodeData(NodeMap(ii))%compma + NodeData(NodeMap(ii))%por * COMPFL


            IF(I.GT.1.AND.NodeMap(ii) .NE.NROLD) NRTEST = NRTEST + 1 
            NROLD = NodeMap(ii) 
         END DO

  460    IF (KNODAL.EQ.0) WRITE (fLST, 461) SCALX, SCALY, SCALZ, PORFAC 
  461 FORMAT(1H1////11X,'N O D E   I N F O R M A T I O N'//16X,         &
     &   'PRINTOUT OF NODE COORDINATES AND POROSITIES ',                &
     &   'CANCELLED.'//16X,'SCALE FACTORS :'/33X,1PD15.4,5X,'X-SCALE'/  &
     &   33X,1PD15.4,5X,'Y-SCALE'/33X,1PD15.4,5X,'Z-SCALE'/             &
     &   33X,1PD15.4,5X,'POROSITY FACTOR')                              
         IF (IUNSAT.EQ.1.AND.KNODAL.EQ.0.AND.NRTEST.NE.1) WRITE (fLST,463)                                                           
         IF (IUNSAT.EQ.1.AND.KNODAL.EQ.0.AND.NRTEST.EQ.1) WRITE (fLST,465)                                                           
  463 FORMAT(33X,'MORE THAN ONE REGION OF UNSATURATED PROPERTIES HAS ', &
     &   'BEEN SPECIFIED AMONG THE NODES.')                             
  465 FORMAT(33X,'ONLY ONE REGION OF UNSATURATED PROPERTIES HAS ',      &
     &   'BEEN SPECIFIED AMONG THE NODES.')                             
         IF (KNODAL.EQ. + 1.AND.(IUNSAT.NE.1.and.fZON==0)) &
           WRITE (fLST, 470) (I, X (I), Y (I), Z (I), NodeData(NodeMap(i))%por, I = 1, NN)                              
  470 FORMAT(1H1//11X,'N O D E   I N F O R M A T I O N'//13X,           &
     &   'NODE',7X,'X',16X,'Y',21X,'Z',10X,'POROSITY'//                 &
     &   (8X,I9,3(3X,1PD14.5),6X,0PF8.5))                               
         IF (KNODAL.EQ. + 1.AND.(IUNSAT.EQ.1.or.fZON/=0)) &
           WRITE (fLST, 480) (I, NodeMap(i), X (I), Y (I), Z (I), NodeData(NodeMap(i))%por, I = 1, NN)
  480 FORMAT(1H1//11X,'N O D E   I N F O R M A T I O N'//13X,'NODE',3X, &
     &   'REGION',7X,'X',16X,'Y',21X,'Z',10X,'POROSITY'//               &
     &   (8X,I9,3X,I6,3(3X,1PD14.5),6X,0PF8.5))                         
!
!.....INPUT DATASET 14C: SORPTION PARAMETERS SPECIFIED BY NODE
!.....INPUT DATASET 14D: PRODUCTION PARAMETERS SPECIFIED BY NODE
      if( fZON == 0 ) then
         CALL ReadSorption14C()
         CALL ReadSorption14D()
      end if
!                                                                       
!.....INPUT DATASETS 15A AND 15B: ELEMENTWISE DATA                      
         MSErrorValue%cDataSet='15B'
         CALL SKPCOM (fINP, NLSKIP) 
         READ (fINP,*,iostat=ios) CDUM10, PMAXFA, PMIDFA, PMINFA, ANG1FA, ANG2FA,  &
         ANG3FA, ALMAXF, ALMIDF, ALMINF, AT1MXF, AT1MDF, AT1MNF
         if(ios/=0) call ErrorIO('Error specifying elementwise multiplier data')
         IF (trim(CDUM10).NE.'ELEMENT') THEN 
           call ErrorIO('Error: Dataset 15A must begin with the word "ELEMENT"')
         ENDIF 
         !Set AT2 multiplies equal to AT1 multipliers
           AT2MXF=AT1MXF
           AT2MDF=AT1MDF
           AT2MNF=AT1MNF

         IF (KELMNT.EQ. + 1) WRITE (fLST, 500) 
  500 FORMAT(1H1//11X,'E L E M E N T   I N F O R M A T I O N'//         &
     &   11X,'ELEMENT',4X,'MAXIMUM',9X,'MIDDLE',10X,'MINIMUM',9X,       &
     &   'ANGLE1',9X,'ANGLE2',9X,'ANGLE3',7X,                           &
     &   'LONGITUDINAL',5X,'LONGITUDINAL',5X,'LONGITUDINAL',5X,         &
     &   'TRANSVERSE',6X,'TRANSVERSE',6X,'TRANSVERSE',6X,               &
     &   'TRANSVERSE',6X,'TRANSVERSE',6X,'TRANSVERSE'/                  &
     &   22X,'PERMEABILITY',4X,'PERMEABILITY',4X,'PERMEABILITY',2X,     &
     &   '(IN DEGREES)',3X,'(IN DEGREES)',3X,'(IN DEGREES)',3X,         &
     &   'DISPERSIVITY',5X,'DISPERSIVITY',5X,'DISPERSIVITY',5X,         &
     &   'DISPERSIVITY_1',2X,'DISPERSIVITY_1',2X,'DISPERSIVITY_1',2X,   &
     &   'DISPERSIVITY_2',2X,'DISPERSIVITY_2',2X,'DISPERSIVITY_2'/)     
         IF (KELMNT.EQ. + 1.AND.IUNSAT.EQ.1) WRITE (fLST, 508) 
  508 FORMAT(14X,'REGION') 
         IF (KELMNT.EQ. + 1) WRITE (fLST, 509) 
  509 FORMAT(/) 

         LRTEST = 1 
         MSErrorValue%cDataSet='15B'
         CALL SKPCOM (fINP, NLSKIP) 
         DO 550 LL = 1, NE 
            READ(fINP,*,iostat=ios) L
            if(ios/=0) call ErrorIO('Error specifying element number - error on line ['//trim(Val2Char(LL))//'] of Data Set 15B')
            if(L>NE)  call ErrorIO('Error: Specified element number ['//trim(Val2Char(L))//'] &
exceeds NE ['//trim(Val2Char(NE))//']')
            backspace(fINP)
            if(fZON==0) then
              READ(fINP,*,iostat=ios) L, ElemMap(l), &
                                             ElemData(L)%pmax, ElemData(L)%pmid, ElemData(L)%pmin, &
                                             ElemData(L)%angle1, ElemData(L)%angle2, ElemData(L)%angle3, &
                                             ElemData(L)%almax, ElemData(L)%almid, ElemData(L)%almin, &
                                             ElemData(L)%at1max, ElemData(L)%at1mid, ElemData(L)%at1min
              ElemMap(l)=L
            end if
            if(fZON/=0) READ(fINP,*,iostat=ios) L, ElemMap(l)
            if(ios/=0) call ErrorIO('Error specifying elementwise data - error on line ['//trim(Val2Char(LL))//'] of Data Set 15B')

            IF (LL.GT.1.AND.ElemMap(l) .NE.LROLD) LRTEST = LRTEST + 1 
            LROLD = ElemMap(l) 

            if(fZON==0) then
              ElemData(ElemMap(l))%pmax = ElemData(ElemMap(l))%pmax * PMAXFA 
              ElemData(ElemMap(l))%pmid = ElemData(ElemMap(l))%pmid * PMIDFA 
              ElemData(ElemMap(l))%pmin = ElemData(ElemMap(l))%pmin * PMINFA 
              ElemData(ElemMap(l))%angle1 = ElemData(ElemMap(l))%angle1 * ANG1FA 
              ElemData(ElemMap(l))%angle2 = ElemData(ElemMap(l))%angle2 * ANG2FA 
              ElemData(ElemMap(l))%angle3 = ElemData(ElemMap(l))%angle3 * ANG3FA 
            end if

            if(fZON==0) then
              ElemData(ElemMap(l))%almax = ElemData(ElemMap(l))%almax * ALMAXF 
              ElemData(ElemMap(l))%almid = ElemData(ElemMap(l))%almid * ALMIDF 
              ElemData(ElemMap(l))%almin = ElemData(ElemMap(l))%almin * ALMINF 
              ElemData(ElemMap(l))%at1max = ElemData(ElemMap(l))%at1max * AT1MXF 
              ElemData(ElemMap(l))%at1mid = ElemData(ElemMap(l))%at1mid * AT1MDF 
              ElemData(ElemMap(l))%at1min = ElemData(ElemMap(l))%at1min * AT1MNF 
              ElemData(ElemMap(l))%at2max = ElemData(ElemMap(l))%at1max
              ElemData(ElemMap(l))%at2mid = ElemData(ElemMap(l))%at1mid
              ElemData(ElemMap(l))%at2min = ElemData(ElemMap(l))%at1min
            end if

            IF(KELMNT.EQ. + 1.AND.(IUNSAT.NE.1.and.fZON==0)) &
              WRITE(fLST, 520) L, &
                ElemData(ElemMap(l))%pmax, ElemData(ElemMap(l))%pmid, ElemData(ElemMap(l))%pmin, &
                ElemData(ElemMap(l))%angle1, ElemData(ElemMap(l))%angle2, ElemData(ElemMap(l))%angle3, &
                ElemData(ElemMap(l))%almax, ElemData(ElemMap(l))%almid, ElemData(ElemMap(l))%almin, &
                ElemData(ElemMap(l))%at1max, ElemData(ElemMap(l))%at1mid, ElemData(ElemMap(l))%at1min, &
                ElemData(ElemMap(l))%at2max, ElemData(ElemMap(l))%at2mid, ElemData(ElemMap(l))%at2min

            IF(KELMNT.EQ. + 1.AND.(IUNSAT.EQ.1.or.fZON==0)) &
              WRITE(fLST, 530) L, ElemMap(l), &
                ElemData(ElemMap(l))%pmax, ElemData(ElemMap(l))%pmid, ElemData(ElemMap(l))%pmin, &
                ElemData(ElemMap(l))%angle1, ElemData(ElemMap(l))%angle2, ElemData(ElemMap(l))%angle3, &
                ElemData(ElemMap(l))%almax, ElemData(ElemMap(l))%almid, ElemData(ElemMap(l))%almin, &
                ElemData(ElemMap(l))%at1max, ElemData(ElemMap(l))%at1mid, ElemData(ElemMap(l))%at1min, &
                ElemData(ElemMap(l))%at2max, ElemData(ElemMap(l))%at2mid, ElemData(ElemMap(l))%at2min


  520 FORMAT(9X,I9,2X,3(1PD14.5,2X),3(0PF10.3,5X),9(1PD14.5,2X)) 
  530 FORMAT(3X,I9,1X,I5,2X,3(1PD14.5,2X),3(0PF10.3,5X),9(1PD14.5,2X)) 
!                                                                       
!.....ROTATE PERMEABILITY FROM MAX/MID/MIN TO X/Y/Z DIRECTIONS.         
!.....THIS SECTION OF CODE (THROUGH THE CALL TO "TENSOR") WAS LIFTED    
!.....FROM DAVE POLLOCK'S "indat13-dwp.f".                              
            if(fZON==0) then
              D2R = 1.745329252D-2 
              ElemData(ElemMap(l))%pangl1 = D2R * ElemData(ElemMap(l))%angle1
              ElemData(ElemMap(l))%pangl2 = D2R * ElemData(ElemMap(l))%angle2 
              ElemData(ElemMap(l))%pangl3 = D2R * ElemData(ElemMap(l))%angle3 
              ZERO = 0D0 
              MSTRUC = 1
              CALL TENSOR (ElemData(ElemMap(l))%pmax, ZERO, ZERO, ZERO, &
                ElemData(ElemMap(l))%pmid, ZERO, ZERO, ZERO, ElemData(ElemMap(l))%pmin, &
                ElemData(ElemMap(l))%pangl1, ElemData(ElemMap(l))%pangl2, ElemData(ElemMap(l))%pangl3, &
                ElemData(ElemMap(l))%permxx, ElemData(ElemMap(l))%permxy, ElemData(ElemMap(l))%permxz, &
                ElemData(ElemMap(l))%permyx, ElemData(ElemMap(l))%permyy, ElemData(ElemMap(l))%permyz, &
                ElemData(ElemMap(l))%permzx, ElemData(ElemMap(l))%permzy, ElemData(ElemMap(l))%permzz, &
                MSTRUC)                 
            end if

  550    END DO 
         IF (KELMNT.EQ.0.and.fZON==0) &
         WRITE (fLST, 569) PMAXFA, PMIDFA, PMINFA,       &
         ANG1FA, ANG2FA, ANG3FA, ALMAXF, ALMIDF, ALMINF, AT1MXF, AT1MDF,&
         AT1MNF, AT2MXF, AT2MDF, AT2MNF                                 
  569 FORMAT(////11X,'E L E M E N T   I N F O R M A T I O N'//          &
     &   16X,'PRINTOUT OF ELEMENT PERMEABILITIES AND DISPERSIVITIES ',  &
     &   'CANCELLED.'//16X,'SCALE FACTORS :'/33X,1PD15.4,5X,'MAXIMUM ', &
     &   'PERMEABILITY FACTOR'/33X,1PD15.4,5X,'MIDDLE PERMEABILITY ',   &
     &   'FACTOR '/33X,1PD15.4,5X,'MINIMUM PERMEABILITY FACTOR'/        &
     &   33X,1PD15.4,5X,'ANGLE1 FACTOR'/33X,1PD15.4,5X,'ANGLE2 FACTOR'/ &
     &   33X,1PD15.4,5X,'ANGLE3 FACTOR'/                                &
     &   33X,1PD15.4,5X,'FACTOR FOR LONGITUDINAL DISPERSIVITY IN ',     &
     &   'MAX-PERM DIRECTION'/33X,1PD15.4,5X,'FACTOR FOR LONGITUDINAL ',&
     &   'DISPERSIVITY IN MID-PERM DIRECTION'/33X,1PD15.4,5X,'FACTOR ', &
     &   'FOR LONGITUDINAL DISPERSIVITY IN MIN-PERM DIRECTION'/         &
     &   33X,1PD15.4,5X,'FACTOR FOR TRANSVERSE DISPERSIVITY_1 IN ',     &
     &   'MAX-PERM DIRECTION'/33X,1PD15.4,5X,'FACTOR FOR TRANSVERSE ',  &
     &   'DISPERSIVITY_1 IN MID-PERM DIRECTION'/33X,1PD15.4,5X,'FACTOR',&
     &   ' FOR TRANSVERSE DISPERSIVITY_1 IN MIN-PERM DIRECTION'/        &
     &   33X,1PD15.4,5X,'FACTOR FOR TRANSVERSE DISPERSIVITY_2 IN ',     &
     &   'MAX-PERM DIRECTION'/33X,1PD15.4,5X,'FACTOR FOR TRANSVERSE ',  &
     &   'DISPERSIVITY_2 IN MID-PERM DIRECTION'/33X,1PD15.4,5X,'FACTOR',&
     &   ' FOR TRANSVERSE DISPERSIVITY_2 IN MIN-PERM DIRECTION')        

         IF (IUNSAT.EQ.1.AND.KELMNT.EQ.0.AND.LRTEST.NE.1) WRITE (fLST,573)                                                           
         IF (IUNSAT.EQ.1.AND.KELMNT.EQ.0.AND.LRTEST.EQ.1) WRITE (fLST,575)                                                           
  573 FORMAT(33X,'MORE THAN ONE REGION OF UNSATURATED PROPERTIES HAS ', &
     &   'BEEN SPECIFIED AMONG THE ELEMENTS.')                          
  575 FORMAT(33X,'ONLY ONE REGION OF UNSATURATED PROPERTIES HAS ',      &
     &   'BEEN SPECIFIED AMONG THE ELEMENTS.')                          
!                                                                       
      ELSE 
!.....READ 2-D INPUT FROM DATASETS 13 - 15.                             
!.....NOTE THAT Z = THICKNESS; AT1MAX = ATMAX; AT1MIN = ATMIN;          
!.....AND PANGL1 = PANGLE.                                              
!                                                                       
!.....INPUT DATASET 13: ORIENTATION OF COORDINATES TO GRAVITY           
         MSErrorValue%cDataSet=' 13'
         CALL SKPCOM (fINP, NLSKIP) 
         READ (fINP,*,iostat=ios) GRAVX, GRAVY 
         if(ios/=0) call ErrorIO('Error specifying gravity for 2D problem')
         GRAVZ = 0D0 
         WRITE (fLST, 1320) GRAVX, GRAVY 
 1320 FORMAT(////11X,'C O O R D I N A T E   O R I E N T A T I O N   ',  &
     &   'T O   G R A V I T Y'//13X,'COMPONENT OF GRAVITY VECTOR',      &
     &   /13X,'IN +X DIRECTION, GRAVX'/11X,1PD15.4,5X,                  &
     &   'GRAVX = -GRAV * D(ELEVATION)/DX'//13X,'COMPONENT OF GRAVITY', &
     &   ' VECTOR'/13X,'IN +Y DIRECTION, GRAVY'/11X,1PD15.4,5X,         &
     &   'GRAVY = -GRAV * D(ELEVATION)/DY')                             
!                                                                       
!.....INPUT DATASETS 14A AND 14B: NODEWISE DATA                         
         MSErrorValue%cDataSet='14A'
         CALL SKPCOM (fINP, NLSKIP) 
         READ (fINP,*,iostat=ios) CDUM10, SCALX, SCALY, SCALTH, PORFAC 
         if(ios/=0) call ErrorIO('Error specifying nodewise multiplier data')
         IF (trim(CDUM10).NE.'NODE') THEN 
           call ErrorIO('Error: Dataset 14A must begin with the word "NODE"')
         ENDIF 

         NRTEST = 1 
         MSErrorValue%cDataSet='14B'
         CALL SKPCOM (fINP, NLSKIP) 
         DO I = 1, NN 
            READ(fINP,*,iostat=ios) II
            if(ios/=0) call ErrorIO('Error specifying node number - error on line ['//trim(Val2Char(I))//'] of Data Set 15B')
            if(I>NN)  call ErrorIO('Error: Specified node number ['//trim(Val2Char(II))//'] exceeds NE ['//trim(Val2Char(NN))//']')
            backspace(fINP)
            if(fZON==0) then
              READ(fINP,*,iostat=ios) II, NodeMap(ii), X(II), Y(II), Z(II), NodeData(II)%por  
              NodeMap(ii) = II
            end if
            if(fZON/=0) READ(fINP,*,iostat=ios) II, NodeMap(ii), X(II), Y(II), Z(II)
            if(ios/=0) call ErrorIO('Error specifying nodewise data - error on line ['//trim(Val2Char(I))//'] of Data Set 14B')
            X(II) = X(II) * SCALX 
            Y(II) = Y(II) * SCALY 
            Z(II) = Z(II) * SCALTH 
            if(fZON==0) NodeData(NodeMap(ii))%por = NodeData(NodeMap(ii))%por * PORFAC 
!           SET SPECIFIC PRESSURE STORATIVITY, SOP.                           
            if(fZON==0) NodeData(NodeMap(ii))%sop = (1.D0 - NodeData(NodeMap(ii))%por ) * &
NodeData(NodeMap(ii))%compma + NodeData(NodeMap(ii))%por * COMPFL 

            IF(I.GT.1.AND.NodeMap(ii) .NE.NROLD) NRTEST = NRTEST + 1 
            NROLD = NodeMap(ii) 

         end do

 1460    IF (KNODAL.EQ.0) WRITE (fLST, 1461) SCALX, SCALY, SCALTH, PORFAC 
 1461 FORMAT(1H1////11X,'N O D E   I N F O R M A T I O N'//16X,         &
     &   'PRINTOUT OF NODE COORDINATES, THICKNESSES AND POROSITIES ',   &
     &   'CANCELLED.'//16X,'SCALE FACTORS :'/33X,1PD15.4,5X,'X-SCALE'/  &
     &   33X,1PD15.4,5X,'Y-SCALE'/33X,1PD15.4,5X,'THICKNESS FACTOR'/    &
     &   33X,1PD15.4,5X,'POROSITY FACTOR')                              
         IF (IUNSAT.EQ.1.AND.KNODAL.EQ.0.AND.NRTEST.NE.1) WRITE (fLST,1463)                                                          
         IF (IUNSAT.EQ.1.AND.KNODAL.EQ.0.AND.NRTEST.EQ.1) WRITE (fLST,1465)                                                          
 1463 FORMAT(33X,'MORE THAN ONE REGION OF UNSATURATED PROPERTIES HAS ', &
     &   'BEEN SPECIFIED AMONG THE NODES.')                             
 1465 FORMAT(33X,'ONLY ONE REGION OF UNSATURATED PROPERTIES HAS ',      &
     &   'BEEN SPECIFIED AMONG THE NODES.')                             
         IF (KNODAL.EQ. + 1.AND.(IUNSAT.NE.1.and.fZON==0)) &
           WRITE (fLST, 1470) (I, X (I), Y (I), Z (I), NodeData(NodeMap(i))%por, I = 1, NN)                              

 1470 FORMAT(1H1//11X,'N O D E   I N F O R M A T I O N'//13X,           &
     &   'NODE',7X,'X',16X,'Y',17X,'THICKNESS',6X,'POROSITY'//          &
     &   (8X,I9,3(3X,1PD14.5),6X,0PF8.5))                               
         IF (KNODAL.EQ. + 1.AND.(IUNSAT.EQ.1.or.fZON/=0)) &
           WRITE (fLST, 1480) (I, NodeMap(i), X (I), Y (I), Z (I), NodeData(NodeMap(i))%por, I = 1, NN)                   

 1480 FORMAT(1H1//11X,'N O D E   I N F O R M A T I O N'//13X,'NODE',3X, &
     &   'REGION',7X,'X',16X,'Y',17X,'THICKNESS',6X,'POROSITY'//        &
     &   (8X,I9,3X,I6,3(3X,1PD14.5),6X,0PF8.5))                         
!
!.....INPUT DATASET 14C: SORPTION PARAMETERS SPECIFIED BY NODE
!.....INPUT DATASET 14D: PRODUCTION PARAMETERS SPECIFIED BY NODE
      if( fZON == 0 ) then
         CALL ReadSorption14C()
         CALL ReadSorption14D()
      end if
!                                                                       
!.....INPUT DATASETS 15A AND 15B: ELEMENTWISE DATA                      
         MSErrorValue%cDataSet='15A'
         CALL SKPCOM (fINP, NLSKIP) 
         READ (fINP,*,iostat=ios) CDUM10, PMAXFA, PMINFA, ANGFAC, ALMAXF, ALMINF, ATMAXF, ATMINF                                                 
         if(ios/=0) call ErrorIO('Error specifying elementwise multiplier data')
         IF (trim(CDUM10).NE.'ELEMENT') THEN 
           call ErrorIO('Error: Dataset 15A must begin with the word "ELEMENT"')
         ENDIF 

         IF (KELMNT.EQ. + 1) WRITE (fLST, 1500) 
 1500 FORMAT(1H1//11X,'E L E M E N T   I N F O R M A T I O N'//         &
     &   11X,'ELEMENT',4X,'MAXIMUM',9X,'MINIMUM',12X,                   &
     &   'ANGLE BETWEEN',3X,'LONGITUDINAL',3X,'LONGITUDINAL',5X,        &
     &   'TRANSVERSE',5X,'TRANSVERSE'/                                  &
     &   22X,'PERMEABILITY',4X,'PERMEABILITY',4X,                       &
     &   '+X-DIRECTION AND',3X,'DISPERSIVITY',3X,'DISPERSIVITY',3X,     &
     &   'DISPERSIVITY',3X,'DISPERSIVITY'/                              &
     &   50X,'MAXIMUM PERMEABILITY',3X,' IN MAX-PERM',                  &
     &   3X,' IN MIN-PERM',3X,' IN MAX-PERM',3X,' IN MIN-PERM'/         &
     &   58X,'(IN DEGREES)',3X,'   DIRECTION',3X,'   DIRECTION',        &
     &   3X,'   DIRECTION',3X,'   DIRECTION')                           
         IF (KELMNT.EQ. + 1.AND.IUNSAT.EQ.1) WRITE (fLST, 1508) 
 1508 FORMAT(14X,'REGION') 
         IF (KELMNT.EQ. + 1) WRITE (fLST, 1509) 
 1509 FORMAT(/) 
         LRTEST = 1 
         MSErrorValue%cDataSet='15B'
         CALL SKPCOM (fINP, NLSKIP) 
         DO 1550 LL = 1, NE 
            READ(fINP,*,iostat=ios) L
            if(ios/=0) call ErrorIO('Error specifying element number - error on line ['//trim(Val2Char(LL))//'] of Data Set 15B')
            if(L>NE)  call ErrorIO('Error: Specified element number ['//trim(Val2Char(L))//'] &
exceeds NE ['//trim(Val2Char(NE))//']')
            backspace(fINP)

            if(fZON==0) then
              READ(fINP,*,iostat=ios) &
                         L, ElemMap(l), &
                         ElemData(L)%pmax, ElemData(L)%pmin, &
                         ElemData(L)%anglex, &
                         ElemData(L)%almax, ElemData(L)%almin, &
                         ElemData(L)%at1max, ElemData(L)%at1min
              ElemMap(l)=L
            end if
            if(fZON/=0) READ(fINP,*,iostat=ios) &
                         L, ElemMap(l)
            if(ios/=0) call ErrorIO('Error specifying elementwise data - error on line ['//trim(Val2Char(LL))//'] of Data Set 15B')

            IF(LL.GT.1.AND.ElemMap(l) .NE.LROLD) LRTEST = LRTEST + 1 
            LROLD = ElemMap(l) 

            if(fZON==0) then
              ElemData(ElemMap(l))%pmax = ElemData(ElemMap(l))%pmax * PMAXFA 
              ElemData(ElemMap(l))%pmin = ElemData(ElemMap(l))%pmin * PMINFA 
              ElemData(ElemMap(l))%anglex = ElemData(ElemMap(l))%anglex * ANGFAC 

              ElemData(ElemMap(l))%almax = ElemData(ElemMap(l))%almax * ALMAXF 
              ElemData(ElemMap(l))%almin = ElemData(ElemMap(l))%almin * ALMINF 
              ElemData(ElemMap(l))%at1max = ElemData(ElemMap(l))%at1max * ATMAXF 
              ElemData(ElemMap(l))%at1min = ElemData(ElemMap(l))%at1min * ATMINF
            end if

             
            IF(KELMNT.EQ. + 1.AND.(IUNSAT.NE.1.and.fZON==0)) &
              WRITE(fLST, 1520) &
                L, &
                ElemData(ElemMap(l))%pmax, ElemData(ElemMap(l))%pmin, &
                ElemData(ElemMap(l))%anglex, &
                ElemData(ElemMap(l))%almax, ElemData(ElemMap(l))%almin, &
                ElemData(ElemMap(l))%at1max, ElemData(ElemMap(l))%at1min
 1520 FORMAT(9X,I9,2X,2(1PD14.5,2X),5X,5(0PF10.3,5X)) 
            IF(KELMNT.EQ. + 1.AND.(IUNSAT.EQ.1.or.fZON/=0)) &
              WRITE(fLST, 1530) &
                L, ElemMap(l), &
                ElemData(ElemMap(l))%pmax, ElemData(ElemMap(l))%pmin, &
                ElemData(ElemMap(l))%anglex, &
                ElemData(ElemMap(l))%almax, ElemData(ElemMap(l))%almin, &
                ElemData(ElemMap(l))%at1max, ElemData(ElemMap(l))%at1min
 1530 FORMAT(3X,I9,1X,I5,2X,2(1PD14.5,2X),5X,5(0PF10.3,5X)) 

            if(fZON==0) then
              RADIAX = 1.745329D-2 * ElemData(ElemMap(l))%anglex 
              SINA = SIN(RADIAX) 
              COSA = COS(RADIAX) 
              SINA2 = SINA * SINA 
              COSA2 = COSA * COSA 
              ElemData(ElemMap(l))%permxx = ElemData(ElemMap(l))%pmax * COSA2 + ElemData(ElemMap(l))%pmin * SINA2 
              ElemData(ElemMap(l))%permyy = ElemData(ElemMap(l))%pmax * SINA2 + ElemData(ElemMap(l))%pmin * COSA2 
              ElemData(ElemMap(l))%permxy =(ElemData(ElemMap(l))%pmax - ElemData(ElemMap(l))%pmin) * SINA * COSA 
              ElemData(ElemMap(l))%permyx = ElemData(ElemMap(l))%permxy 
              ElemData(ElemMap(l))%pangl1 = RADIAX 
            end if


 1550    END DO 
         IF (KELMNT.EQ.0.and.fZON==0) & 
           WRITE (fLST, 1569) PMAXFA, PMINFA, ANGFAC, ALMAXF, ALMINF, ATMAXF, ATMINF                                 
 1569 FORMAT(////11X,'E L E M E N T   I N F O R M A T I O N'//          &
     &   16X,'PRINTOUT OF ELEMENT PERMEABILITIES AND DISPERSIVITIES ',  &
     &   'CANCELLED.'//16X,'SCALE FACTORS :'/33X,1PD15.4,5X,'MAXIMUM ', &
     &   'PERMEABILITY FACTOR'/33X,1PD15.4,5X,'MINIMUM PERMEABILITY ',  &
     &   'FACTOR'/33X,1PD15.4,5X,'ANGLE FROM +X TO MAXIMUM DIRECTION',  &
     &   ' FACTOR'/33X,1PD15.4,5X,'FACTOR FOR LONGITUDINAL DISPERSIVITY'&
     &  ,' IN MAX-PERM DIRECTION'/33X,1PD15.4,5X,                       &
     &   'FACTOR FOR LONGITUDINAL DISPERSIVITY IN MIN-PERM DIRECTION',  &
     &   /33X,1PD15.4,5X,'FACTOR FOR TRANSVERSE DISPERSIVITY',          &
     &   ' IN MAX-PERM DIRECTION'/33X,1PD15.4,5X,                       &
     &   'FACTOR FOR TRANSVERSE DISPERSIVITY IN MIN-PERM DIRECTION')    
         IF (IUNSAT.EQ.1.AND.KELMNT.EQ.0.AND.LRTEST.NE.1) WRITE (fLST,1573)                                                          
         IF (IUNSAT.EQ.1.AND.KELMNT.EQ.0.AND.LRTEST.EQ.1) WRITE (fLST,1575)                                                          
 1573 FORMAT(33X,'MORE THAN ONE REGION OF UNSATURATED PROPERTIES HAS ', &
     &   'BEEN SPECIFIED AMONG THE ELEMENTS.')                          
 1575 FORMAT(33X,'ONLY ONE REGION OF UNSATURATED PROPERTIES HAS ',      &
     &   'BEEN SPECIFIED AMONG THE ELEMENTS.')                          
!                                                                       
      ENDIF 
!                                                                       
!.....INPUT DATASETS 15C AND 15D [IF REQUIRED]: ELEMENTWISE DATA        
!.....TEST IF MULTIPLE SPECIES ARE BEING SIMULATED                      
!     ATSPMULT IS INITIALIZED TO 1 SO ONLY THE SPECIES WITH             
!     DISPERSIVITY MULTIPLIERS OTHER THAN 1.0 NEED TO BE DEFINED        
      MSErrorValue%cDataSet='15C'
      IF (NSPE.GT.1) THEN 
         CALL SKPCOM (fINP, NLSKIP) 
 1600    READ (fINP,*,iostat=ios) IMULT 
         if(ios/=0) call ErrorIO('Error specifying species number for dispersivity multiplier')
         
!.......TEST FOR TERMINATION OF DATASET 15C                             
         IF (IMULT.EQ.0) GOTO 1605 
!.......TEST IF L EXCEEDS THE NUMBER OF SPECIES                         
         IF (IMULT.GT.NSPE) THEN 
            WRITE (fLST, 1610) IMULT, NSPE 
            GOTO 1600 
         ENDIF 
         BACKSPACE (fINP) 
         READ (fINP,*,iostat=ios) IMULT, ATSPMULT(IMULT) 
         if(ios/=0) call ErrorIO('Error specifying species number for dispersivity multiplier for species'&
//trim(Val2Char(IMULT))//' ['//trim(SPNAME(IMULT))//']')
         GOTO 1600 
!.......WRITE ATSPMULT TO UNIT fLST                                        
 1605    WRITE (fLST, 1615) (IMULT, trim(adjustl(SPNAME(IMULT))), ATSPMULT(IMULT), IMULT = 1, NSPE)                                               
      ENDIF 
 1610 FORMAT(//2X'W A R N I N G'/2X'ERROR IN DATASET 15C'               &
     &        /2X'SPECIES',1XI3,1X'EXCEEDS THE TOTAL NUMBER OF SPECIES',&
     &            ' [',1XI3,1X'] SPECIFIED IN DATASET 3'                &
     &        /2X'I G N O R I N G  E N T R Y')                          
 1615 FORMAT(//2X'D I S P E R S I V I T Y  M U L T I P L I E R'         &
     &       //2X'ALL DISPERSIVITY VALUES FOR A GIVEN SPECIES'          &
     &        /2X'[MIN, MAX, AND MID] SPECIFIED IN DATASET 15B'         &
     &        /2X'FOR EACH ELEMENT ARE MULTIPLIED BY THE'               &
     &        /2X'DISPERSIVITY MULTIPLIER'                              &
     &       //2X'   SPECIES',2X'NAME      ',2X,5X,'MULTIPLER'          &
     &        /48('-')/(2X,I10,2X,A10,2X1PD15.8))                       
!                                                                       
!.....INPUT DATASET 15D [IF SIGMAS IS LESS THAN 0.0D0]                  
!.....INITIALIZE SLAMBDA TO DABS(SIGMAS)                                
!     ONLY THE ELEMENTS THAT ARE NOT EQUAL TO DABS(SIGMAS)              
!     NEED TO BE DEFINED                                                
      if(fZON==0) then
        do ILAMB=1,ElemZones
          ElemData(ilamb)%lambdas = DABS (SIGMAS) 
        end do
      end if
!                                                                       
!.....TEST IF SLAMBDA IS TO BE READ
      MSErrorValue%cDataSet='15D'
      IF ( SIGMAS.LT.0.0D0.AND.ME.GT. - 1 ) THEN 
         CALL SKPCOM (fINP, NLSKIP) 
         read (fINP,*,iostat=ios) CDUM10
         if(ios/=0) call ErrorIO('Error specifying bulk thermal conductivity equation type')
         if(trim(CDUM10).eq.'AVERAGE'.or.trim(CDUM10).eq.'average') then
           LVolAvgLambda=.true.
           write(fLST,'(//a//)') 'Volumetric Average Bulk Thermal Conductivity Equation Used for Heat Transport Equation'
         else if(trim(CDUM10).eq.'GEOMETRIC'.or.trim(CDUM10).eq.'geometric') then
           LVolAvgLambda=.false.
           write(fLST,'(//a//)') 'Geometric Mean Bulk Thermal Conductivity Equation Used for Heat Transport Equation'
         !assumes that only average or geometric keywords are present
         else
           backspace(fINP)
         end if

         lModifyLambdas = .FALSE.
         CALL SKPCOM (fINP, NLSKIP) 
 1650    READ (fINP,*,iostat=ios) ILAMB 
         if(ios/=0) call ErrorIO('Error specifying element number for lambda &
modification.  Confirm that CEQ is either "AVERAGE" or "GEOMETRIC"')
         if(ILAMB>NE)  call ErrorIO('Error: Specified element number ['//trim(Val2Char(ILAMB))//'] &
exceeds NE ['//trim(Val2Char(NE))//']')
!.......TEST FOR TERMINATION OF DATASET 15D                             
         IF (ILAMB.EQ.0) GOTO 1655 
         BACKSPACE (fINP)
         if(fZON==0) READ (fINP,*,iostat=ios) ILAMB, ElemData(ElemMap(ILAMB))%lambdas
         if(ios/=0) call ErrorIO('Error specifying lambda modification for species'//trim(Val2Char(IMULT))//'&
 ['//trim(SPNAME(IMULT))//']')
         lModifyLambdas = .TRUE.

         GOTO 1650 
!.......WRITE SLAMBDA TO UNIT fLST                                        
 1655   if(lModifyLambdas) then
          WRITE (fLST, 1660) 
          WRITE (fLST, 1665) (ILAMB, ElemData(ElemMap(ILAMB))%lambdas, ILAMB = 1, NE) 
        end if

      ENDIF 
!                                                                       
 1660 FORMAT(////52X'  E N E R G Y  T R A N S P O R T  '/               &
     &           55X'  L A M B D A  S O L I D  '/132('-'))              
 1665 FORMAT(//132('-'),                                                &
     &       //2X,5(6X,'NODE',16X)/(2X,5(1X,I9,1X,1PD15.8)))            
!                                                                       
!.....SET SIGMAS TO ZERO IF COUPLED SOLUTE/ENERGY TRANSPORT             
!     SINCE SIGMAS VALUES FOR ENERGY TRANSPORT HAVE BEEN PUSHED         
!     INTO SLAMBDA ABOVE                                                
!     ALLOWS MINIMAL MODIFICATION TO EQUATIONS IN ELEMN3 AND            
!     ELEMN2 SUBROUTINES                                                
      IF (ME.EQ.0) SIGMAS = 0.0D0 
!
!.....Find observation nodes if using Specified Observation Locations (SOB)
      if(fSOB>0) then
          MSErrorValue%cDataSet='SB3'
        if(.not.CalcObsNode()) call ErrorIO('Error calculating observation nodes from SOB x,y,z data')
      end if

!                                                                       
!                                                                       
 1000 RETURN 
      END SUBROUTINE INDAT1                         
