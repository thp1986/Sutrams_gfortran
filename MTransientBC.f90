!     MODULE            T  B  C  S             SUTRA-MS VERSION 2004.1
                                                                       
! ***
! *** PURPOSE :                                                         
! ***  MODULE THAT IS USED TO READ TRANSIENT BOUNDARY CONDITION (TBC OPTION)
! ***  DATA (READTBCDATA).  THE MODULE ALSO GIVES THE PROGRAM ACCESS TO THE
! ***  THE VARIABLES LGETTBC AND LRDTBCTIME WHICH ARE USED TO INSTRUCT THE
! ***  THE FUNCTION READ.  THE GLOBAL VARIABLE DNTIME IS THE TIME NEW TBC
! ***  WILL BE READ.
! ***

      module TBCS
        !MS specific modules
        USE SutraMSPrecision
        USE MSErrorHandler
        implicit none
        logical                           :: lGetTBC=.FALSE.
        logical                           :: lRdTBCTime=.TRUE.
        real  (DP)                        :: DNTIME
        public :: &
          lGetTBC, lRdTBCTime, DNTIME, &
          ReadTBCData

        contains

          !function to read Transient BC data
          logical function ReadTBCData()
            USE SutraStorage, ONLY: QIN, UIN, IQSOP, QUIN, IQSOU, IPBC, PBC, IUBC, UBC, &
                                    SpecifiedPBC, MultiSpeciesBC
            USE ITERAT 
            USE PARAMS 
            USE CONTRL 
            USE FUNITS 
            USE DIMS
            USE TIMES
            implicit none
            !local
            CHARACTER (LEN=20) :: &
              CTERM
            INTEGER (I4B) :: &
              NLSKIP
            INTEGER (I4B) :: &
              NTPBC, &
              NTUBC, &
              NTSOP, &
              NTSOU
            INTEGER (I4B) :: &
              I, &
              K
            INTEGER (I4B) :: &
              IU, &
              IQU, &
              ID, &
              IP, &
              IUT, &
              IPU, &
              IQCP, &
              IQP
            INTEGER (I4B) :: &
              NIQP, &
              NIQU
            REAL (DP) :: &
              UINC(NSPE)
            REAL (DP) :: &
              PP, &
              UU, &
              PUU(NSPE), &
              QINC, &
              QUINC

      !
            lOk=.false.
      !.....FORMAT STATEMENTS                                                 
       2100 FORMAT(////11X,'*** T R A N S I E N T  B O U N D A R Y  ',        &
           &     'C O N D I T I O N S',                                       &
           &/11X,'**** NODES AT WHICH PRESSURE IS ',                          &
           &     'SPECIFIED TO BE INDEPENDENT OF LOCAL FLOWS AND FLUID ',     &
           &     'SOURCES ****')                                              
       2110 FORMAT(//1X,'SPECIFIED PRESSURE NODE ('I10,') '  &
           &            'WAS NOT DEFINED AS A SPECIFIED NODE '                &
           &            'IN DATASET 19')                                      
       2120 FORMAT(11X,I10,13X,99(1PE14.7)) 
       2200 FORMAT(////11X,'*** T R A N S I E N T  B O U N D A R Y  ',        &
           &     'C O N D I T I O N S',                                       &
           &/11X,'**** NODES AT WHICH U CONCENTRATIONS ARE ',                 &
           &     'SPECIFIED TO BE INDEPENDENT OF LOCAL FLOWS AND FLUID ',     &
           &     'SOURCES ****')                                              
       2210 FORMAT(//1X,'SPECIFIED CONCENTRATION/TEMPERATURE NODE ('I10,') '  &
           &            'WAS NOT DEFINED AS A SPECIFIED NODE '                &
           &            'IN DATASET 20')                                      
       2220 FORMAT(11X,I10,I10,3X,1PE14.7) 
       2230 FORMAT(11X,I10,13X,1PE14.7) 
       2300 FORMAT(1H1////11X,'         T R A N S I E N T'                    &
           &   /11X,'F L U I D   S O U R C E   D A T A'                       &
           &////11X,'**** NODES AT WHICH FLUID INFLOWS OR OUTFLOWS ARE ',     &
           &'SPECIFIED ****'//11X,'NODE NUMBER',10X,                          &
           &'FLUID INFLOW(+)/OUTFLOW(-)'//)                                   
       2310 FORMAT(//1X,'SPECIFIED FLUID SOURCE/SINK NODE ('I10,') WAS NOT '  &
           &            'DEFINED AS A SOURCE/SINK NODE IN DATASET 17')        
       2320 FORMAT(11X,I10,13X,1PE14.7,99(:16X,1PE14.7)) 
       2330 FORMAT(//1X,'NUMBER OF FLUID SOURCE/SINKS (',I5,') '              &
           &            'EXCEEDS NSOP (',I5,')',                              &
           &        /1X,'CORRECT fTBC DATASET AND RERUN')                       
       2400 FORMAT(1H1////11X,'         T R A N S I E N T'                    &
           &   /11X,'S O L U T E   M A S S  /  T E M P E R A T U R E'         &
           &   /11X 'S O U R C E  D A T A'                                    &
           &////11X,'**** NODES AT WHICH SOLUTE/ENERGY INFLOWS OR OUTFLOWS ', &
           &'ARE SPECIFIED ****'//11X,'NODE NUMBER',10X,                      &
           &'FLUID INFLOW(+)/OUTFLOW(-)'//)                                   
       2410 FORMAT(//1X,'SPECIFIED SOLUTE/ENERGY SOURCE/SINK NODE ('I10,') '  &
           &            'WAS NOT DEFINED AS A SOURCE/SINK NODE IN DATASET 18')        
      !
            MSErrorValue%cDataSet='TBC'
      !.....READ TIME TO READ NEXT TRANSIENT DATA SET                         
            IF(lRdTBCTime) THEN
              CALL SKPCOM (fTBC, NLSKIP)
              READ (fTBC, *, IOSTAT = ios) CTERM 
              if(ios/=0) call ErrorIO('Error specifying dataset 1 (next starting time) for transient boundary conditions')
              IF(TRIM(CTERM)=='END'.OR.TRIM(CTERM)=='end') THEN
                DNTIME=TMAX+1.0D0
                GOTO 999
              END IF
              READ (CTERM, *, IOSTAT = ios) DNTIME 
              if(ios/=0) call ErrorIO('Internal Error converting CTERM ['//trim(adjustl(CTERM))//'] to DNTIME')
              GOTO 999
            END IF

      !
      !.....READ TRANSIENT BOUNDARY CONDITION DATA
            CALL SKPCOM (fTBC, NLSKIP) 
            READ (fTBC, *, IOSTAT = ios) NTPBC, NTUBC, NTSOP, NTSOU 
            if(ios/=0) call ErrorIO('Error specifying dataset 2 for transient boundary conditions')
      !                                                                       
      !.....TEST FOR TRANSIENT PARAMETERS THAT EXCEED SPECIFIED               
      !     DIMENSIONS                                                        
            IF (NTPBC.GT.NPBC) &
               call ErrorIO('NTPBC ['//trim(Val2Char(NTPBC))//'] EXCEEDS NPBC ['//trim(Val2Char(NPBC))//']')

            IF (NTSOU.GT. (MNSOU - 1) ) &
               call ErrorIO('NTSOU ['//trim(Val2Char(NTSOU))//'] EXCEEDS NSOU ['//trim(Val2Char(MNSOU-1))//']')
      !                                                                       
      !.....READ REVISED TRANSIENT BOUNDARY CONDITIONS                        
      !                                                                       
      !.....SPECIFIED PRESSURE BOUNDARY CONDITIONS                            
100 &
            IF (NTPBC.LT.1) GOTO 200
            IU = 0
            WRITE (fLST, 2100) 
            CALL SKPCOM (fTBC, NLSKIP) 
110 &
            IUT = 0 
            READ (fTBC,*,IOSTAT=ios) ID 
            if(ios/=0) call ErrorIO('Error specifying node number for specified pressure boundary condition')
      !.....TEST THAT I IS NOT EQUAL TO ZERO                                  
            IF (ID.EQ.0) GOTO 140 
            BACKSPACE (fTBC) 
            READ (fTBC,*,IOSTAT=ios) ID, PP, (PUU(K),K=1,NSPE)
            if(ios/=0) call ErrorIO('Error specifying specified pressure boundary condition data')

      !.....TEST THAT NODE WAS SPECIFIED IN DATASET 19                        
            DO IP = 1, NPBC 
               IF (SpecifiedPBC(IP)%node .EQ. ID) THEN 
                  WRITE ( *, * ) ID, SpecifiedPBC(IP)%node, IP, NPBC
                  IUT = IP 
                  EXIT 
               ENDIF 
            END DO 

            IF (IUT.EQ.0) THEN 
               WRITE (fLST, 2110) ABS (ID) 
               GOTO 110 
            ENDIF 

            IU = IU + 1 
            SpecifiedPBC(IUT)%P = PP
            DO K=1,NSPE
             SpecifiedPBC(IUT)%U(K) = PUU(K)
            END DO
            WRITE (fLST, 2120) ID, PP, (PUU(K),K=1,NSPE)
      !.....continue reading specified concentration/temperature nodes until 0 encountered
            GOTO 110 

      !.....TEST FOR ERROR IN NUMBER OF SPECIFIED PRESSURE NODES SPECIFIED IN UNIT fTBC                                        
140 &
            IF (IU.GT.NPBC) & 
               call ErrorIO('Actual NTPBC ['//trim(Val2Char(IU))//'] EXCEEDS NPBC ['//trim(Val2Char(NPBC))//']')


      !                                                                       
      !.....SPECIFIED CONCENTRATION/TEMPERATURE BOUNDARY CONDITIONS           
200 &
            IF (NTUBC.LT.1) GOTO 300 
            WRITE (fLST, 2200) 
            IU = 0 
            CALL SKPCOM (fTBC, NLSKIP) 
210 &
            READ (fTBC,*,IOSTAT=ios) I 
            if(ios/=0) call ErrorIO('Error specifying node number for specified concentration/temperature boundary condition')
      !.....TEST THAT I IS NOT EQUAL TO ZERO                                  
            IF (I.EQ.0) GOTO 240 
            BACKSPACE (fTBC) 
            K = 1 

            IF (NSPE.LT.2) THEN
              READ (fTBC,*,IOSTAT=ios) ID, UU
            ELSE 
              READ (fTBC,*,IOSTAT=ios) K, ID, UU 
            END IF
            if(ios/=0) call ErrorIO('Error specifying specified concentration/temperature boundary condition data')

      !.....TEST THAT NODE WAS SPECIFIED IN DATASET 20                        
            IUT = 0 
            DO IPU = 1, NUBC(K)
               IF (MultiSpeciesBC(K)%SpecifiedU(IPU)%node .EQ.ID) THEN 
                  IUT = IPU 
                  EXIT 
               ENDIF 
            END DO 

            IF (IUT.EQ.0) THEN 
               WRITE (fLST, 2210) ABS (ID) 
               GOTO 210 
            ENDIF 

            IU = IU + 1 
            MultiSpeciesBC(K)%SpecifiedU(IUT)%U = UU 
            IF (NSPE.LT.2) THEN
              WRITE (fLST, 2230) ID, UU 
            ELSE
              WRITE (fLST, 2220) K, ID, UU 
            END IF
      !.....continue reading specified concentration/temperature nodes until 0 encountered
            GOTO 210 

      !.....TEST FOR ERROR IN NUMBER OF SPECIFIED CONCENTRATION/TEMPERATURE   
      !     NODES SPECIFIED IN UNIT fTBC                                        
240 &
            IF (IU.GT.MNUBC) & 
               call ErrorIO('Actual NTUBC ['//trim(Val2Char(IU))//'] EXCEEDS MNUBC ['//trim(Val2Char(MNUBC))//']')
      !                                                                       
      !.....SPECIFIED FLUID SOURCE AND SINK BOUNDARY CONDITIONS               
300 &
            IF (NTSOP.LT.1) GOTO 400 
            WRITE (fLST, 2300) 
            NIQP = 0 
310 &
            CONTINUE 
            CALL SKPCOM (fTBC, NLSKIP) 
320 &
            READ (fTBC,*,IOSTAT=ios) IQCP 
            if(ios/=0) call ErrorIO('Error specifying node number for specified fluid source/sink boundary condition')
            IF (IQCP.GT.0) THEN 
               BACKSPACE (fTBC) 
               READ (fTBC,*,IOSTAT=ios) IQCP, QINC 
               if(ios/=0) call ErrorIO('Error specifying specified [Q] fluid source/sink boundary condition data')
      !........FLUID SOURCE WITH SPECIFIED CONCENTRATIONS                     
               IF (QINC.GT.0D0) THEN 
                  BACKSPACE (fTBC) 
                  READ (fTBC,*,IOSTAT=ios) IQCP, QINC, (UINC (K), K = 1, NSPE) 
                  if(ios/=0) call ErrorIO('Error specifying specified [U] fluid source/sink boundary condition data')
               ENDIF 
            ELSEIF (IQCP.EQ.0) THEN 
               GOTO 350 
            ENDIF 
            IQP = 0 
      !.....TEST THAT NODE WAS SPECIFIED IN DATA SET 17                       
            DO I = 1, (NSOP - 1) 
               IF (IQSOP (I) .EQ.ABS (IQCP) ) THEN 
                  IQP = ABS (IQCP) 
                  EXIT 
               ENDIF 
            END DO 
            IF (IQP.EQ.0) THEN 
               WRITE (fLST, 2310) ABS (IQCP) 
               GOTO 320 
            ENDIF 
            NIQP = NIQP + 1 
            QIN (IQP) = QINC 
      !.....SPECIFIED CONCENTRATIONS FOR POSITIVE SOURCE [INJECTION]          
            IF (QINC.GT.0D0) THEN 
               DO K = 1, NSPE 
                 UIN (IQP, K) = UINC (K) 
               ENDDO 
            ENDIF 
      !.....TEST IF POSITIVE SOURCE [INJECTION]                               
            IF (QINC.GT.0) GOTO 340 
      !.....WRITE NODE, AND RATE FOR FLUID SINK                               
            WRITE (fLST, 2320) IQCP, QINC 
            GOTO 320 
      !.....WRITE NODE,RATE,AND SPECIFIED CONCENTRATION OF FIRST SPECIES      
340 &
            WRITE (fLST, 2320) IQCP, QINC, (UINC (K), K = 1, NSPE) 
      !.....CONTINUE READING FLUID SOURCES UNTIL IQCP = 0                     
            GOTO 320 

      !.....TEST FOR ERROR IN NUMBER OF FLUID SOURCE SINKS SPECIFIED          
      !     IN UNIT fTBC                                                        
350 &
            IF (NIQP.GT. (NSOP - 1) ) & 
              call ErrorIO('Actual NTSOP ['//trim(Val2Char(NIQP))//'] EXCEEDS NSOP ['//trim(Val2Char(NSOP-1))//']')
      !                                                                       
      !.....SPECIFIED SOLUTE/ENERGY SOURCE AND SINK BOUNDARY CONDITIONS       
400 &
            IF (NTSOU.LT.1) GOTO 500 
            WRITE (fLST, 2400) 
            NIQU = 0 
410 &
            CONTINUE 
            CALL SKPCOM (fTBC, NLSKIP) 
420 &
            READ (fTBC,*,IOSTAT=ios) I 
            if(ios/=0) call ErrorIO('Error specifying node number for specified concentration/temperature boundary condition')
      !.....TEST THAT I IS NOT EQUAL TO ZERO                                  
            IF (I.EQ.0) GOTO 440 
            BACKSPACE (fTBC) 
            K = 1 

            IF (NSPE.LT.2) THEN
              READ (fTBC,*,IOSTAT=ios) ID, QUINC
            ELSE 
              READ (fTBC,*,IOSTAT=ios) K, ID, QUINC 
            END IF
            if(ios/=0) call ErrorIO('Error specifying specified concentration/energy source/sink boundary condition data')

      !.....TEST THAT NODE WAS SPECIFIED IN DATASET 18
            IQU = 0 
            DO I = 1, (MNSOU - 1) 
               IF (IQSOU (I, K) .EQ.ID) THEN 
                  IQU = ABS (ID) 
                  EXIT 
               ENDIF 
            END DO 
            IF (IQU.EQ.0) THEN 
               WRITE (fLST, 2410) ABS (ID) 
               GOTO 420 
            ENDIF 
            NIQU = NIQU + 1 
            QUIN (IQU,K) = QUINC 

            IF (NSPE.LT.2) THEN
              WRITE (fLST, 2230) ID, QUINC
            ELSE
              WRITE (fLST, 2220) K, ID, QUINC 
            END IF
      !.....continue reading specified solute/energy source nodes until 0 encountered
            GOTO 420 

      !.....TEST FOR ERROR IN NUMBER OF SPECIFIED CONCENTRATION/ENERGY SOURCE/SINK NODES SPECIFIED IN UNIT fTBC                                        
440 &
            IF (NIQU.GT.(MNSOU-1)) & 
               call ErrorIO('Actual NTSOU ['//trim(Val2Char(NIQU))//'] EXCEEDS MNSOU ['//trim(Val2Char(MNSOU-1))//']')
                                                                        
500 &
            GOTO 999 
      !                                                                       
      !.....RETURN TO CALLING ROUTINE                                         
999 & 
            lOk=.true.
            ReadTBCData=lOk
            RETURN 
            end function ReadTBCData

      end module  TBCS
