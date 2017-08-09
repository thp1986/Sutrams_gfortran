!     SUBROUTINE        O  T  T  O  B  S       SUTRA-MS VERSION 2004.1

! *** PURPOSE :                                                         
! ***  TO PRINT THE SOLUTION AT OBSERVATION NODES.  SPECIFICALLY,       
! ***  TO PRINT PRESSURES OR HEADS, CONCENTRATIONS OR TEMPERATURES,     
! ***  AND SATURATIONS IN A COLUMNWISE FORMAT.  OUTPUT IS TO UNIT fOBS.   
!                                                                       
      SUBROUTINE OUTOBS (IOBS, X, Y, Z, PVEC, UVEC, SW, TITLE1, TITLE2) 

      USE ITERAT 
      USE CONTRL 
      USE SOLVI 
      USE PARAMS 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE OBS
      USE GRAVEC
      USE PLT1
      USE SutraMSPrecision

      IMPLICIT NONE

      real (DP) :: &
        X (NN), Y (NN), Z (NN) 
      real (DP) :: &
        PVEC (NN), UVEC (NN, NSPE), SW (NN) 
      integer (I4B) :: &
        IOBS (NOBSN) 
      CHARACTER(1) TITLE1 (80), TITLE2 (80) 
      CHARACTER(8) HORP 
      CHARACTER (LEN = 13) ::TORC (NSPE), CORC (NSPE) 
      CHARACTER (LEN = 2) ::CNSPE 
      CHARACTER (LEN = 10) ::CNE, CNS 
      CHARACTER (LEN = 516) ::C964 
      CHARACTER (LEN = 256) ::C965, C966, C967, C968, C969, C980 
      real (DP) :: &
        TT (9999)
      integer (I4B) :: &
        ITT (9999), ISTORC (9999), ISHORP (9999), ISSATU (9999)

      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        IO, &
        NSTART, NEND, &
        NOBS, &
        LCHORP, LCTORC, &
        JJ, JT, JTMAX, &
        K, KT, KTMAX
      REAL (DP) :: &
        DKTM2, DELTK, &
        TS

      IF (.NOT.ONCEK7) THEN 
!.....BUILD FORMAT STATEMENTS FOR MULTI SPECIES OBSERVATION OUTPUT      
         NSTART = (NSPE-1) / 2 
         NEND = (NSPE-NSTART - 1) 
         NSTART = NSTART * 15 
         NEND = NEND * 15 
         WRITE (CNSPE, '(I2)') (NSPE-1) 
         WRITE (CNS, '(I10)') NSTART 
         WRITE (CNE, '(I10)') NEND 
!.....C964 FORMAT                                                       
         C964 = '("## ", 77("="),/"## ", A, 25X, I9, " Time steps printed",'                                                                 
         C964 = TRIM (ADJUSTL (C964) ) //'/"## ", 77("="),/"## ",' 
         C964 = TRIM (ADJUSTL (C964) ) //'/"##  Time steps", 27X,' 
         C964 = TRIM (ADJUSTL (C964) ) //'"[Latest time step computed]",' 
         C964 = TRIM (ADJUSTL (C964) ) //'/"## in this file    Time (sec)",'                                                                 
         IF (NSPE.GT.1) THEN 
            C964 = TRIM (ADJUSTL (C964) ) //'9X,A5, 10X,A4,'//TRIM(ADJUSTL(CNSPE) )                                           
            C964 = TRIM (ADJUSTL (C964) ) //'(11X,A4), 11X,A3,' 
            C964 = TRIM (ADJUSTL (C964) ) //'/"## ", 12("-"), 3X, 12("-"), 1X,'
            C964 = TRIM (ADJUSTL (C964) ) //'3(3X,12("-"))'//TRIM(ADJUSTL(CNSPE))
            C964 = TRIM (ADJUSTL (C964) ) //'(3X, 12("-")) )' 
         ELSE 
            C964 = TRIM (ADJUSTL (C964) ) //'9X,A5, 10X,A4, 11X,A3,' 
            C964 = TRIM (ADJUSTL (C964) ) //'/"## ", 12("-"), 3X, 12("-"), 1X,'                                                                 
            C964 = TRIM (ADJUSTL (C964) ) //' 3(3X, 12("-")) )' 
         ENDIF 
!.....C965 FORMAT                                                       
         IF (NSPE.GT.1) THEN 
            C965 = '("## ", 3X, I8, 3X, 1PE13.6, 3(7X, I8),' 
            C965 = TRIM(ADJUSTL(C965))//TRIM(ADJUSTL(CNSPE))//'(7X, I8))'                                                               
         ELSE 
            C965 = '("## ", 3X, I8, 3X, 1PE13.6, 3(7X, I8))' 
         ENDIF 
!.....C966 FORMAT                                                       
         C966 = '("## ",/"## ", 77("="),/"## ",/"## ", 24X, 999(:3X, ' 
         IF ( NSPE.GT.1 ) THEN
              IF ( NSTART.GT.0 ) THEN 
                  C966 = TRIM(ADJUSTL(C966))//TRIM(ADJUSTL(CNS))//'X,' 
              END IF
              C966 = TRIM(ADJUSTL(C966))//'16X, "NODE ", I9, 15X,' 
              C966 = TRIM(ADJUSTL(C966))//TRIM(ADJUSTL(CNE))//'X))'
         ELSE 
              C966 = TRIM(ADJUSTL(C966))//'16X, "NODE ", I9, 15X))' 
         ENDIF 
!.....C967 FORMAT                                                       
         IF (NSPE.GT.1) THEN 
!              C967 = '("## ",25X,999(:2X,'//TRIM(ADJUSTL(CNS))//'X,' 
              C967 = '("## ",25X,999(:2X,'
              IF ( NSTART.GT.0 ) THEN 
                  C967 = TRIM(ADJUSTL(C967))//TRIM(ADJUSTL(CNS))//'X,' 
              END IF
              C967 = TRIM(ADJUSTL(C967))//' "(",2(1PE14.7,","),1PE14.7,")",'                                                               
              C967 = TRIM(ADJUSTL(C967))//TRIM(ADJUSTL(CNE))//'X))' 
         ELSE 
              C967 = '("## ",25X,999(:2X,"(",2(1PE14.7,","),1PE14.7,")"))' 
         ENDIF 
!.....C968 FORMAT                                                       
         IF(NSPE.GT.1)THEN 
            C968 = '("## ",24X,999(:3X,1X,A14,2A15,' 
            C968 = TRIM(ADJUSTL(C968))//TRIM(ADJUSTL(CNSPE))//'A15))'                                                     
         ELSE 
            C968 = '("## ",24X,999(:3X,1X,A14,2A15))' 
         ENDIF 
!.....C969 FORMAT                                                       
         IF(NSPE.GT.1)THEN 
            C969 = '("## ","Time Step",5X,"Time(sec)",999(:10X,A,2X,' 
            C969 = TRIM(ADJUSTL(C969))//'A,'//TRIM(ADJUSTL(CNSPE))                                                          
            C969 = TRIM(ADJUSTL(C969))//'(2X,A),5X,"Saturation"))' 
         ELSE 
            C969 = '("## ","Time Step",5X,"Time(sec)",999(:10X,A,2X,' 
            C969 = TRIM(ADJUSTL(C969))//'A,5X,"Saturation"))' 
         ENDIF 
!.....C980 FORMAT                                                       
         IF(NSPE.GT.1)THEN 
            C980 = '(3X,I9,1PE15.7,999(:3X,3(1PE15.7),' 
            C980 = TRIM(ADJUSTL(C980))//TRIM(ADJUSTL(CNSPE))//'(1PE15.7)))'                                                             
         ELSE 
            C980 = '(3X,I9,1PE15.7,999(:3X,3(1PE15.7)))' 
         ENDIF 
!.....END OF FORMAT BUILD IF (.NOT. ONCEK7)                             
      ENDIF 
!                                                                       
!.....Calculate headers on time step 1                                  
!.....and create output on time steps greater than or equal to 1        
      IF (IT.EQ.0) RETURN 
!                                                                       
      NOBS = NOBSN - 1 
      DKTM2 = DBLE (IABS (KTYPE) - 2) 
!                                                                       
!.....Post-processing information, nodewise and velocity results        
!.....File Header for nodewise and velocity data files                  
      IF (.NOT.ONCEK7) THEN 
!                                                                       
         IF (NOBSN - 1.EQ.0) THEN 
!...........TEST IF OBSERVATION FILE IS OPEN                            
!           IF NOT PRINT NO OBSERVATION NODE TEXT TO OUTPUT FILE        
!           ON UNIT fLST                                                  
            IO = fOBS 
            IF (fOBS.EQ.0) IO = fLST 
            WRITE (IO, 960) TITLE1, TITLE2 
            WRITE (IO, 5) 
    5 FORMAT      (/'  *** NO OBSERVATION NODES SPECIFIED (NOBS=0) ***') 
            ONCEK7 = .TRUE. 
            RETURN 
         ENDIF 
!                                                                       
!.....   Calculate and print Time Step Information                      
         TS = TSTART 
!.....   Time step value                                                
         JT = 0 
!.....   Number of printed time steps                                   
         KT = 0 
!.....   Time step increment                                            
         DELTK = DELT 
!.....   Indicators of when variables are calculated and printed        
         LCHORP = 0 
         LCTORC = 0 
   10    CONTINUE 
         JT = JT + 1 
         IF (MOD (JT, ITCYC) .EQ.0.AND.JT.GT.1) DELTK = DELTK * DTMULT 
         IF (DELTK.GT.DTMAX) DELTK = DTMAX 
         TS = TS + DELTK 
         IF (MOD (JT, NPCYC) .EQ.0.OR.JT.EQ.1) LCHORP = JT 
         IF (MOD (JT, NUCYC) .EQ.0.OR.JT.EQ.1) LCTORC = JT 
         IF (MOD (JT, NOBCYC) .EQ.0.OR.JT.EQ.1) THEN 
            KT = KT + 1 
            TT (KT) = TS 
            ITT (KT) = JT 
            ISHORP (KT) = LCHORP 
            ISTORC (KT) = LCTORC 
            ISSATU (KT) = LCHORP 
         ENDIF 
         IF (JT.LT.ITMAX.AND.TS.LT.TMAX) GOTO 10 
         JTMAX = JT 
!                                                                       
!                                                                       
!.....Print last time step always, unless we already printed it         
         IF (JTMAX.GT.1.AND.MOD (JT, NOBCYC) .NE.0) THEN 
            KT = KT + 1 
            TT (KT) = TS 
            ITT (KT) = JT 
            IF (MOD (JT, NPCYC) .EQ.0) LCHORP = JT 
            IF (MOD (JT, NUCYC) .EQ.0) LCTORC = JT 
            ISHORP (KT) = LCHORP 
            ISTORC (KT) = LCTORC 
            ISSATU (KT) = LCHORP 
         ENDIF 
!.....   Number of printed time steps                                   
         KTMAX = KT 
!                                                                       
!.....If steady-state flow, P and S calculated on t.s. "0" only.        
         IF (ISSFLO.EQ.2) THEN 
            DO 14 KT = 1, KTMAX 
               ISHORP (KT) = 0 
               ISSATU (KT) = 0 
   14       END DO 
         ENDIF 
!                                                                       
!.....Set Temperature or Concentration text string for headers          
!.....TORC - FOR SUMMARY HEADING                                        
!     CORC - SPECIES NAMES DEFINED IN INPUT DATASET 3A                  
!            IN UNIT fNOD                                                 
         DO K = 1, NSPE 
           TORC (K) = "Concentration" 
!..........SET TITLE FOR ENERGY TRANSPORT SPECIES                       
           IF (K.EQ.NESP) TORC(K) = "Temperature" 
           CORC (K)  = '   '//ADJUSTR(SPNAME (K) ) 
         ENDDO 
!                                                                       
!.....Set Pressure text string for header                               
         HORP = "Pressure" 
!                                                                       
!.....Write the header information                                      
         IF (fOBS.GT.0) THEN 
            WRITE (fOBS, 960) TITLE1, TITLE2 
            IF (KTYPE.LT.0) THEN 
               IF (IABS (KTYPE) .EQ.3) THEN 
                  WRITE (fOBS, 961) IABS (KTYPE), NN1, NN2, NN3, NN, "    &
                  Nodes"                                                
               ELSE 
                  WRITE (fOBS, 962) IABS (KTYPE), NN1, NN2, NN, "Nodes" 
               ENDIF 
            ELSE 
               WRITE (fOBS, 963) IABS (KTYPE), NN, "Nodes" 
            ENDIF 
            WRITE (fOBS, C964) "OBSERVATIONNODERESULTS", KTMAX, HORP,     &
            (TORC (K), K = 1, NSPE), "Sat"                              
            DO 20 KT = 1, KTMAX 
               WRITE (fOBS, C965) ITT (KT), TT (KT), ISHORP (KT), (ISTORC &
               (KT), K = 1, NSPE), ISSATU (KT)                          
   20       END DO 
            NOBS = NOBSN - 1 
            WRITE (fOBS, C966) (IOBS (JJ), JJ = 1, NOBS) 
            WRITE (fOBS, C968) ('---------------', JJ = 1, (3 + NSPE-1)   &
            * NOBS)                                                     
            WRITE (fOBS, C967) (X (IOBS (JJ) ), Y (IOBS (JJ) ), Z (IOBS ( &
            JJ) ), JJ = 1, NOBS)                                        
            WRITE (fOBS, C968) ('---------------', JJ = 1, (3 + NSPE-1)   &
            * NOBS)                                                     
            WRITE (fOBS, C969) (HORP, (CORC (K), K = 1, NSPE), JJ = 1,    &
            NOBS)                                                       
         ENDIF 
  960 FORMAT   ("## ", 80A1,                                            &
     &         /"## ", 80A1,                                            &
     &         /"## ")                                                  
  961 FORMAT   ("## ", I1, "-D, REGULAR MESH  ", 5X,                    &
     &                 "(", 2(I9, ")*("), I9, ") = ", I9, A,            &
     &         /"## ")                                                  
  962 FORMAT   ("## ", I1, "-D, REGULAR MESH  ", 17X,                   &
     &                 "(", I9, ")*(", I9, ") = ", I9, A,               &
     &         /"## ")                                                  
  963 FORMAT   ("## ", I1, "-D, IRREGULAR MESH", 43X, I9, A,            &
     &         /"## ")                                                  
                                                                        
         ONCEK7 = .TRUE. 
      ENDIF 
!                                                                       
      IF (NOBSN - 1.EQ.0) RETURN 
!                                                                       
      IF (fOBS.GT.0) THEN 
         WRITE(fOBS, C980)IT, TSEC,(PVEC(IOBS(JJ)),(UVEC(IOBS(JJ),K), K = 1, NSPE), SW(IOBS(JJ)), JJ = 1, NOBS)            
      ENDIF 
!                                                                       
!                                                                       
      RETURN 
!                                                                       
      END SUBROUTINE OUTOBS
