!     SUBROUTINE        OUTVELOCITY            SUTRA-MS VERSION 2004.1

! *** PURPOSE :                                                         
! ***  TO PRINT ELEMENT CENTROID COORDINATES AND VELOCITY COMPONENTS    
! ***  IN A FLEXIBLE, COLUMNWISE FORMAT.  OUTPUT IS TO UNIT fELE.         
!                                                                       
      SUBROUTINE OUTVELOCITY(VMAG, VANG1, VANG2, IN, X, Y, Z, TITLE1, TITLE2) 

      USE ITERAT 
      USE CONTRL 
      USE SOLVI 
      USE JCOLS 
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
        VMAG (NE), VANG1 (NE), VANG2 (NEX) 
      real (DP) :: &
        X (NN), Y (NN), Z (NN) 
      real (DP) :: &
        VCOL (NCOLMX), VVAR (7) 
      CHARACTER(1) TITLE1 (80), TITLE2 (80) 
      CHARACTER(15) COLTK6 (7) 
      CHARACTER(1) CPVX, CPVY, CPVZ 
      integer (I4B) :: &
        IN (NIN), IIN (8)
      real (DP) :: &
        TT (999)
      integer (I4B) :: &
        ITT (999) , &
        ISVEL (999) 
      LOGICAL PRINTE 

      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        II, III, &
        JT, JTMAX, &
        KT, KTMAX, &
        L, LL, LCP, LCV, &
        M, MM
      REAL (DP) :: &
        CENTRX, CENTRY, CENTRZ, CVA2, &
        DKTM2, DELTK, &
        RN48, &
        TS, &
        VA1, VA2, VECTRX, VECTRY, VECTRZ


      DATA (COLTK6 (MM) , MM = 1, 7)  / 'Element', &
          '       X origin', '       Y origin', '       Z origin', &
          '     X velocity', '     Y velocity', '     Z velocity' /                                       
      SAVE COLTK6 

!                                                                       
!.....Calculate headers on time step 1                                  
!.....and create output on time steps greater than or equal to 1        
      IF (IT.EQ.0) RETURN 
!                                                                       
      DKTM2 = DBLE (IABS (KTYPE) - 2) 
!                                                                       
!.....Post-processing information, nodewise and velocity results        
!.....File Header for nodewise and velocity data files                  
      IF (.NOT.ONCEK6) THEN 
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
         LCP = 0 
         CPVX = 'N' 
         CPVY = 'N' 
         CPVZ = 'N' 
         DO 8 M = 1, NCOLS6 
            IF (J6COL (M) .EQ.5) CPVX = 'Y' 
            IF (J6COL (M) .EQ.6) CPVY = 'Y' 
            IF (J6COL (M) .EQ.7) CPVZ = 'Y' 
    8    END DO 
   10    CONTINUE 
         JT = JT + 1 
         IF (MOD (JT, ITCYC) .EQ.0.AND.JT.GT.1) DELTK = DELTK * DTMULT 
         IF (DELTK.GT.DTMAX) DELTK = DTMAX 
         TS = TS + DELTK 
         LCV = LCP 
         IF (MOD (JT, NPCYC) .EQ.0.OR.JT.EQ.1) LCP = JT 
         IF (MOD (JT, LCOLPR) .EQ.0.OR.JT.EQ.1) THEN 
            KT = KT + 1 
            TT (KT) = TS 
            ITT (KT) = JT 
            ISVEL (KT) = LCV 
            IF (JT.NE.1.OR.ISSFLO.EQ.2) THEN 
               IF (MOD (JT, NUCYC) .NE.0) ISVEL (KT) = 0 
            ENDIF 
         ENDIF 
         IF (JT.LT.ITMAX.AND.TS.LT.TMAX) GOTO 10 
         JTMAX = JT 
!                                                                       
!                                                                       
!.....Print last time step always, unless we already printed it         
         IF (JTMAX.GT.1.AND.MOD (JT, LCOLPR) .NE.0) THEN 
            KT = KT + 1 
            TT (KT) = TS 
            ITT (KT) = JT 
            ISVEL (KT) = LCP 
         ENDIF 
!.....   Number of printed time steps                                   
         KTMAX = KT 
!                                                                       
!.....If steady-state flow, V based on time step "0" only, and output   
!.....occurs only on time step 1.                                       
         IF (ISSFLO.EQ.2) THEN 
            KTMAX = 1 
            ISVEL (1) = 0 
         ENDIF 
!.....Write the header information                                      
         IF (fELE.GT.0) THEN 
            WRITE (fELE, 960) TITLE1, TITLE2 
            IF (KTYPE.LT.0) THEN 
               IF (IABS (KTYPE) .EQ.3) THEN 
                  WRITE (fELE, 961) IABS (KTYPE), NN1 - 1, NN2 - 1, NN3 - &
                  1, NE, "Elements"                                     
               ELSE 
                  WRITE (fELE, 962) IABS (KTYPE), NN1 - 1, NN2 - 1, NE, " &
                  Elements"                                             
               ENDIF 
            ELSE 
               WRITE (fELE, 963) IABS (KTYPE), NE, "Elements" 
            ENDIF 
            WRITE (fELE, 964) "VELOCITYRESULTS", KTMAX, "Vx", "Vy", "Vz" 
            DO 20 KT = 1, KTMAX 
               WRITE (fELE, 965) ITT (KT), TT (KT), CPVX, ISVEL (KT),     &
               CPVY, ISVEL (KT), CPVZ, ISVEL (KT)                       
   20       END DO 
         ENDIF 
  960 FORMAT   ("## ", 80A1,                                            &
     &         /"## ", 80A1,                                            &
     &         /"## ")                                                  
  961 FORMAT   ("## ", I1, "-D, REGULAR MESH  ", 2X,                    &
     &                 "(", 2(I9, ")*("), I9, ") = ", I9, A,            &
     &         /"## ")                                                  
  962 FORMAT   ("## ", I1, "-D, REGULAR MESH  ", 14X,                   &
     &                 "(", I9, ")*(", I9, ") = ", I9, A,               &
     &         /"## ")                                                  
  963 FORMAT   ("## ", I1, "-D, IRREGULAR MESH", 40X, I9, A,            &
     &         /"## ")                                                  
  964 FORMAT   ("## ", 77("="),                                         &
     &         /"## ", A, 33X, I9, " Time steps printed",               &
     &         /"## ", 77("="),                                         &
     &         /"## ",                                                  &
     &         /"##  Time steps", 20X,                                  &
     &                 "[Printed? / Time step on which V is based]"     &
     &         /"## in this file    Time (sec)",10X,A2, 13X,A2, 13X,A2, &
     &         /"## ", 12("-"), 3X, 12("-"), 1X, 3(3X, 12("-")) )       
  965 FORMAT    ("## ", 3X, I8, 3X, 1PE13.6, 3(5X, A1, 1X, I8)) 
                                                                        
         ONCEK6 = .TRUE. 
      ENDIF 
!                                                                       
!.....   Output velocities for steady flow only on the first time step  
      IF ( (ISSFLO.EQ.2) .AND. (IT.GT.1) ) GOTO 9999 
      IF (fELE.GT.0) THEN 
!                                                                       
!.....  Velocity header information repeated before each time step      
         WRITE (fELE, 966) IT, DELT, TSEC 
  966 FORMAT   ('## ',                                                  &
     &         /'## ', 77('='),                                         &
     &         /'## TIME STEP ', I6, 9X, 'Duration: ', 1PD11.4, ' sec', &
     &                               6X, 'Time: ', 1PD11.4, ' sec',     &
     &         /'## ', 77('='))                                         
         PRINTE = (J6COL (1) .EQ.1) 
         IF (PRINTE) THEN 
            WRITE (fELE, 982) (COLTK6 (J6COL (M) ), M = 1, NCOLS6) 
  982 FORMAT       ("## ", A8, 19(A15)) 
         ELSE 
            WRITE (fELE, 983) COLTK6 (J6COL (1) ) (3:15), (COLTK6 (J6COL (&
            M) ), M = 2, NCOLS6)                                        
  983 FORMAT       ("## ", A13, 19(A15)) 
         ENDIF 
!                                                                       
!.....  The velocity data for this time step                            
         RN48 = 1D0 / DBLE (N48) 
         DO 50 L = 1, NE 
            CENTRX = 0D0 
            CENTRY = 0D0 
            CENTRZ = 0D0 
            DO 40 II = 1, N48 
               III = II + (L - 1) * N48 
               IIN (II) = IN (III) 
               CENTRX = CENTRX + X (IIN (II) ) 
               CENTRY = CENTRY + Y (IIN (II) ) 
               CENTRZ = CENTRZ + Z (IIN (II) ) 
   40       END DO 
            CENTRX = CENTRX * RN48 
            CENTRY = CENTRY * RN48 
            CENTRZ = CENTRZ * RN48 
            VA1 = 0.017453292D0 * VANG1 (L) 
            LL = MIN (L, NEX) 
            VA2 = 0.017453292D0 * VANG2 (LL) * DKTM2 
            CVA2 = DCOS (VA2) 
            VECTRX = VMAG (L) * DCOS (VA1) * CVA2 
            VECTRY = VMAG (L) * DSIN (VA1) * CVA2 
            VECTRZ = VMAG (L) * DSIN (VA2) 
            VVAR (1) = DBLE (L) 
            VVAR (2) = CENTRX 
            VVAR (3) = CENTRY 
            VVAR (4) = CENTRZ 
            VVAR (5) = VECTRX 
            VVAR (6) = VECTRY 
            VVAR (7) = VECTRZ 
            DO 984 M = 1, NCOLS6 
               VCOL (M) = VVAR (J6COL (M) ) 
  984       END DO 
            IF (PRINTE) THEN 
               WRITE (fELE, 985) L, (VCOL (M), M = 2, NCOLS6) 
  985 FORMAT          (I9, 2X, 19(1PE15.7)) 
            ELSE 
               WRITE (fELE, 986) (VCOL (M), M = 1, NCOLS6) 
  986 FORMAT          (1X, 20(1PE15.7)) 
            ENDIF 
   50    END DO 
      ENDIF 
!                                                                       
 9999 CONTINUE 
      RETURN 
!                                                                       
      END SUBROUTINE OUTVELOCITY
