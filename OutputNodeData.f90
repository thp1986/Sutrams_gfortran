!     SUBROUTINE        O  U  T  N  O  D  E    SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO PRINT NODE COORDINATES, PRESSURES OR HEADS, CONCENTRATIONS OR 
! ***  TEMPERATURES, AND SATURATIONS IN A FLEXIBLE, COLUMNWISE FORMAT.  
! ***  OUTPUT IS TO UNIT fNOD.                                            
!                                                                       
      SUBROUTINE OUTNODE(PVEC, UVEC, SW, IN, X, Y, Z, TITLE1, TITLE2) 
      USE ITERAT 
      USE CONTRL 
      USE SOLVI 
      USE JCOLS 
      USE PARAMS 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE GRAVEC
      USE PLT1
      USE SutraMSPrecision
      IMPLICIT NONE
      real (DP) :: &
        PVEC (NN), UVEC (NN, NSPE), SW (NN) 
      real (DP) :: &
        X (NN), Y (NN), Z (NN) 
      real (DP) :: &
        VCOL (NCOLMX), VVAR (7 + NSPE-1) 
      CHARACTER(1) TITLE1 (80), TITLE2 (80) 
      CHARACTER(8) HORP 
      CHARACTER(13) TORC (NSPE) 
      CHARACTER(1) CPHORP, CPTORC, CPSATU 
      integer (I4B) :: &
        IN (NIN), IIN (8)
      real (I4B) :: &
        TT (999)
      integer (I4B) :: &
        ITT (999), ISTORC (999), ISHORP (999), &
        ISSATU (999) 
      LOGICAL PRINTN 
      !LOCAL VARIABLES
      CHARACTER (LEN=15) &
        TCHAR
      INTEGER (I4B) :: &
        LCHAR
      INTEGER (I4B) :: &
        TS, &
        LCHORP, LCTORC, &
        I, &
        JT, JTMAX, &
        K, KK, KT, KTMAX, &
        M
      REAL (DP) :: &
        DKTM2, DELTK


      COLTK5 (1)  = 'Node' 
      COLTK5 (2)  = '              X' 
      COLTK5 (3)  = '              Y' 
      COLTK5 (4)  = '              Z' 
      COLTK5 (5)  = '       Pressure' 
      KK = 0 
      DO K = 6, (5 + NSPE) 
        KK = KK + 1 
        TCHAR = TRIM( SPNAME (KK) )
        LCHAR = LEN_TRIM( TCHAR )
        IF ( LEN_TRIM( TCHAR ) .LT. 15 ) THEN
          TCHAR = ADJUSTR (TCHAR)
          DO I = 1, ( 15 - LCHAR )
            TCHAR (I:I) = ' '
          END DO
        END IF
        COLTK5 (K)  = TCHAR
      ENDDO 
      COLTK5 (6 + NSPE)  = '     Saturation' 
      DO K = 1, (7 + NSPE-1) 
      ENDDO 
!                                                                       
!.....Calculate headers on time step 1                                  
!.....and create output on time steps greater than or equal to 1        
      IF (IT.EQ.0) RETURN 
!                                                                       
      DKTM2 = DBLE (IABS (KTYPE) - 2) 
!                                                                       
!.....Post-processing information, nodewise and velocity results        
!.....File Header for nodewise and velocity data files                  
      IF (.NOT.ONCEK5) THEN 
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
         CPHORP = 'N' 
         CPTORC = 'N' 
         CPSATU = 'N' 
         DO 8 M = 1, NCOLS5 
            IF (J5COL (M) .EQ.5) CPHORP = 'Y' 
            IF (J5COL (M) .EQ.6) CPTORC = 'Y' 
            IF (J5COL (M) .EQ. (7 + NSPE-1) .AND. IUNSAT.EQ.1) CPSATU = 'Y' 
    8    END DO 
   10    CONTINUE 
         JT = JT + 1 
         IF (MOD (JT, ITCYC) .EQ.0.AND.JT.GT.1) DELTK = DELTK * DTMULT 
         IF (DELTK.GT.DTMAX) DELTK = DTMAX 
         TS = TS + DELTK 
         IF (MOD (JT, NPCYC) .EQ.0.OR.JT.EQ.1) LCHORP = JT 
         IF (MOD (JT, NUCYC) .EQ.0.OR.JT.EQ.1) LCTORC = JT 
         IF (MOD (JT, NCOLPR) .EQ.0.OR.JT.EQ.1) THEN 
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
         IF (JTMAX.GT.1.AND.MOD (JT, NCOLPR) .NE.0) THEN 
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
!.....Set Temperature or Concentration text string for header           
         IF (ME.GT.0) THEN 
            TORC (1) = "Temperature" 
         ELSEIF (ME.LT.1) THEN 
            IF (NSPE.EQ.1) THEN 
               TORC (1) = "Concentration" 
            ELSE 
               TORC = "Concentration" 
               IF (NESP.GT.0) TORC (NESP) = "Temperature" 
            ENDIF 
         ENDIF 
!                                                                       
!.....Set Pressure text string for header                               
         HORP = "Pressure" 
!                                                                       
!.....Write the header information                                      
         IF (fNOD.GT.0) THEN 
            WRITE (fNOD, 960) TITLE1, TITLE2 
            IF (KTYPE.LT.0) THEN 
               IF (IABS (KTYPE) .EQ.3) THEN 
                  WRITE (fNOD, 961) IABS (KTYPE), NN1, NN2, NN3, NN, "    &
                  Nodes"                                                
               ELSE 
                  WRITE (fNOD, 962) IABS (KTYPE), NN1, NN2, NN, "Nodes" 
               ENDIF 
            ELSE 
               WRITE (fNOD, 963) IABS (KTYPE), NN, "Nodes" 
            ENDIF 
            WRITE (fNOD, 964) "NODEWISERESULTS", KTMAX, HORP, (TRIM(ADJUSTL(TORC(K))),K=1,NSPE), "Sat"                                         
            DO 20 KT = 1, KTMAX 
               WRITE (fNOD, 965) ITT (KT), TT (KT), CPHORP, ISHORP (KT),  &
               (CPTORC, ISTORC (KT), K = 1, NSPE), CPSATU, ISSATU (KT)  
   20       END DO 
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
  964 FORMAT   ("## ", 77("="),                                         &
     &         /"## ", A, 33X, I9, " Time steps printed",               &
     &         /"## ", 77("="),                                         &
     &         /"## ",                                                  &
     &         /"##  Time steps", 21X,                                  &
     &                 "[Printed? / Latest time step computed]",        &
     &         /"## in this file    Time (sec)", 9X,A5, 10(11X,A4),     &
     &         /"## ", 12("-"), 3X, 12("-"), 1X, 11(3X, 12("-")) )      
  965 FORMAT    ("## ", 3X, I8, 3X, 1PE13.6, 11(5X, A1, 1X, I8)) 
                                                                        
         ONCEK5 = .TRUE. 
      ENDIF 
                                                                        
      IF (fNOD.GT.0) THEN 
!                                                                       
!..... Nodewise Information repeated before each time step              
         WRITE (fNOD, 966) IT, DELT, TSEC 
  966 FORMAT   ('## ',                                                  &
     &         /'## ', 77('='),                                         &
     &         /'## TIME STEP ', I6, 9X, 'Duration: ', 1PD11.4, ' sec', &
     &                               6X, 'Time: ', 1PD11.4, ' sec',     &
     &         /'## ', 77('='))                                         
         PRINTN = (J5COL (1) .EQ.1) 
         IF (PRINTN) THEN 
            WRITE (fNOD, 968) ( COLTK5 (J5COL (M) ), M = 1, NCOLS5) 
  968 FORMAT       ("## ", 2X, A4, 19(A15)) 
         ELSE 
            WRITE (fNOD, 969) COLTK5 (J5COL (1) ) (3:15), &
                   ( COLTK5 (J5COL (M) ), M = 2, NCOLS5)                                        
  969 FORMAT       ("## ", A13, 19(A15)) 
         ENDIF 
!                                                                       
!.....  The nodewise data for this time step                            
         DO 978 I = 1, NN 
            VVAR (1) = DBLE (I) 
            VVAR (2) = X (I) 
            VVAR (3) = Y (I) 
            VVAR (4) = Z (I) 
            VVAR (5) = PVEC (I) 
            KK = 0 
            DO K = 6, (5 + NSPE) 
            KK = KK + 1 
            VVAR (K) = UVEC (I, KK) 
            ENDDO 
            VVAR (6 + NSPE) = SW (I) 
            DO 972 M = 1, NCOLS5 
               VCOL (M) = VVAR (J5COL (M) ) 
  972       END DO 
            IF (PRINTN) THEN 
               WRITE (fNOD, 975) I, (VCOL (M), M = 2, NCOLS5) 
  975 FORMAT          (I9, 19(1PE15.7)) 
            ELSE 
               WRITE (fNOD, 976) (VCOL (M), M = 1, NCOLS5) 
  976 FORMAT          (1X, 20(1PE15.7)) 
            ENDIF 
  978    END DO 
      ENDIF 
!                                                                       
      RETURN 
!                                                                       
      END SUBROUTINE OUTNODE
