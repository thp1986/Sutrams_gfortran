!     SUBROUTINE        OUTPTCLST_3D           SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO PRINT PRESSURE AND TEMPERATURE OR CONCENTRATION               
! ***  SOLUTIONS AND TO OUTPUT INFORMATION ON TIME STEP, ITERATIONS,    
! ***  SATURATIONS, AND FLUID VELOCITIES FOR 3-D PROBLEMS.              
! ***  OUTPUT IS TO UNIT fLST.                                            
!                                                                       
      SUBROUTINE OUTPTCLST_3D(ML, ISTOP, IGOI, PVEC, UVEC, VMAG, VANG1, VANG2, SW)                                                        
      USE ITERAT 
      USE CONTRL 
      USE PARAMS 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE SutraMSPrecision

      IMPLICIT NONE

      INTEGER (I4B) :: &
        ML, ISTOP, IGOI
      real (DP) :: &
        PVEC (NN), UVEC (NN, NSPE), VMAG (NE), VANG1 (NE), VANG2 (NEX), SW (NN)

      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        I, &
        K, &
        L
!                                                                       
!.....OUTPUT MAJOR HEADINGS FOR CURRENT TIME STEP                       
      IF (IT.GT.0.OR.ISSFLO.EQ.2.OR.ISSTRA.EQ.1) GOTO 100 
      WRITE (fLST, 60) 
   60 FORMAT(1H1////11X,'I N I T I A L   C O N D I T I O N S',          &
     &             /11X,'___________________________________')          
      IF (IREAD.EQ. - 1) WRITE (fLST, 65) 
   65 FORMAT(//11X,'INITIAL CONDITIONS RETRIEVED FROM STORAGE ',        &
     &   'ON UNIT fICS.')                                                 
      GOTO 500 
!                                                                       
  100 WRITE (fLST, 350) IT 
  350 FORMAT(1H1//11X,'RESULTS FOR TIME STEP ',I4/                      &
     &   11X,'_______ ___ ____ ____ ____')                              
!                                                                       
      IF (ITRMAX.EQ.1) GOTO 500 
      IF (ISTOP.GE.0) WRITE (fLST, 355) ITER 
  355 FORMAT(/11X,'CONVERGED AFTER ',I3,' ITERATIONS') 
      IF (ISTOP.EQ. - 1) WRITE (fLST, 356) ITER 
  356 FORMAT(/11X,'N O T  CONVERGED AFTER ',I3,' ITERATIONS') 
      WRITE (fLST, 450) RPM, IPWORS 
  450 FORMAT(/11X,'MAXIMUM P CHANGE FROM PREVIOUS ITERATION ',          &
     &   1PD14.5,' AT NODE ',I9)                                        
      DO 455 K = 1, NSPE 
  455 WRITE (fLST, 460) K, trim(adjustl(SPNAME(K))), RUM (K), IUWORS (K) 
  460 FORMAT(/11X,'[SPECIES ',I3,'] - ',A,' - MAXIMUM CHANGE ',       &
     &   'FROM PREVIOUS ITERATION ',1PD14.5,' AT NODE ',I9)             
!                                                                       
  500 IF (IT.EQ.0.AND.ISSFLO.EQ.2) GOTO 680 
      IF (ISSTRA.EQ.1) GOTO 800 
      WRITE (fLST, 550) DELT, TSEC, TMIN, THOUR, TDAY, TWEEK, TMONTH,     &
      TYEAR                                                             
  550 FORMAT(///11X,'TIME INCREMENT :',T27,1PD15.4,' SECONDS'//11X,     &
     &   'ELAPSED TIME :',T27,1PD15.4,' SECONDS',/T27,1PD15.4,' MINUTES'&
     &   /T27,1PD15.4,' HOURS'/T27,1PD15.4,' DAYS'/T27,1PD15.4,' WEEKS'/&
     &   T27,1PD15.4,' MONTHS'/T27,1PD15.4,' YEARS')                    
!                                                                       
!.....OUTPUT PRESSURES FOR TRANSIENT FLOW SOLUTION (AND POSSIBLY,       
!        SATURATION AND VELOCITY)                                       
      IF (ML.EQ.2.AND.ISTOP.GE.0) GOTO 700 
      IF (ISSFLO.GT.0) GOTO 700 
      WRITE (fLST, 650) (I, PVEC (I), I = 1, NN) 
  650 FORMAT(///11X,'P  R  E  S  S  U  R  E'                            &
     &   //2X,5(6X,'NODE',16X)/(2X,5(1X,I9,1X,1PD15.8)))                
      IF (IUNSAT.NE.0) WRITE (fLST, 651) (I, SW (I), I = 1, NN) 
  651 FORMAT(///11X,'S  A  T  U  R  A  T  I  O  N'                      &
     &   //2X,5(6X,'NODE',16X)/(2X,5(1X,I9,1X,1PD15.8)))                
!.....   WRITE 3-D VELOCITIES.                                          
      IF (KVEL.EQ.1.AND.IT.GT.0) THEN 
         WRITE (fLST, 655) (L, VMAG (L), L = 1, NE) 
         WRITE (fLST, 656) (L, VANG1 (L), L = 1, NE) 
         WRITE (fLST, 657) (L, VANG2 (L), L = 1, NE) 
      ENDIF 
  655 FORMAT(///11X,'F  L  U  I  D     V  E  L  O  C  I  T  Y'//        &
     &   11X,'M A G N I T U D E   AT CENTROID OF ELEMENT'//             &
     &   2X,5(3X,'ELEMENT',16X)/(2X,5(1X,I9,1X,1PD15.8)))               
  656 FORMAT(///11X,'F  L  U  I  D     V  E  L  O  C  I  T  Y'//        &
     &   11X,'A N G L E 1   AT CENTROID OF ELEMENT, IN DEGREES FROM ',  &
     &   '+X-AXIS TO PROJECTION OF FLOW DIRECTION IN XY-PLANE'//        &
     &   2X,5(3X,'ELEMENT',16X)/(2X,5(1X,I9,1X,1PD15.8)))               
  657 FORMAT(///11X,'F  L  U  I  D     V  E  L  O  C  I  T  Y'//        &
     &   11X,'A N G L E 2   AT CENTROID OF ELEMENT, IN DEGREES FROM ',  &
     &   'XY-PLANE TO FLOW DIRECTION'//                                 &
     &   2X,5(3X,'ELEMENT',16X)/(2X,5(1X,I9,1X,1PD15.8)))               
!     END IF                                                            
      GOTO 700 
!                                                                       
!.....OUTPUT PRESSURES FOR STEADY-STATE FLOW SOLUTION                   
  680 WRITE (fLST, 690) (I, PVEC (I), I = 1, NN) 
  690 FORMAT(///11X,'S  T  E  A  D  Y  -  S  T  A  T  E     P  R  E  S',&
     &   '  S  U  R  E'//2X,5(6X,'NODE',16X)/(2X,5(1X,I9,1X,1PD15.8)))  
      IF (IUNSAT.NE.0) WRITE (fLST, 651) (I, SW (I), I = 1, NN) 
      GOTO 1000 
!                                                                       
!.....OUTPUT CONCENTRATIONS OR TEMPERATURES FOR                         
!        TRANSIENT TRANSPORT SOLUTION                                   
  700 IF (ML.EQ.1.AND.ISTOP.GE.0) GOTO 1000 
      DO 740 K = 1, NSPE 
         IF (K.EQ.NESP) GOTO 730 
  720    WRITE (fLST, 725) K, trim(adjustl(SPNAME(K))), (I, UVEC (I, K), I = 1, NN) 
  725 FORMAT(///11X,'C  O  N  C  E  N  T  R  A  T  I  O  N'             &
     &         /20X,'S P E C I E S [',I3,']'/22XA,                      &
     &   //2X,5(6X,'NODE',16X)/(2X,5(1X,I9,1X,1PD15.8)))                
         GOTO 740 
  730    WRITE (fLST, 735) K, trim(adjustl(SPNAME(K))), (I, UVEC (I, K), I = 1, NN) 
  735 FORMAT(///11X,'T  E  M  P  E  R  A  T  U  R  E'                   &
     &         /17X,'S P E C I E S [',I3,']'/22XA,                      &
     &   //2X,5(6X,'NODE',16X)/(2X,5(1X,I9,1X,1PD15.8)))                
  740 END DO 
      GOTO 900 
!                                                                       
!.....OUTPUT CONCENTRATIONS OR TEMPERATURES FOR                         
!        STEADY-STATE TRANSPORT SOLUTION                                
  800 DO 840 K = 1, NSPE 
         IF (K.EQ.NESP) GOTO 830 
  820    WRITE (fLST, 825) K, trim(adjustl(SPNAME(K))), (I, UVEC (I, K), I = 1, NN) 
  825 FORMAT(///11X,'S  T  E  A  D  Y  -  S  T  A  T  E     C  O  N  C',&
     &   '  E  N  T  R  A  T  I  O  N'                                  &
     &   /11X,'S P E C I E S [',I3,']'/11XA,                            &
     &   //2X,5(6X,'NODE',16X)/(2X,5(1X,I9,1X,1PD15.8)))                
         GOTO 840 
  830    WRITE (fLST, 835) K, trim(adjustl(SPNAME(K))), (I, UVEC (I, K), I = 1, NN) 
  835 FORMAT(///11X,'S  T  E  A  D  Y  -  S  T  A  T  E     T  E  M  P',&
     &   '  E  R  A  T  U  R  E'                                        &
     &   /11X,'S P E C I E S [',I3,']'/11XA,                            &
     &   //2X,5(6X,'NODE',16X)/(2X,5(1X,I9,1X,1PD15.8)))                
  840 END DO 
!                                                                       
!.....OUTPUT VELOCITIES FOR STEADY-STATE FLOW SOLUTION                  
  900 IF (ISSFLO.NE.2.OR.IT.NE.1.OR.KVEL.NE.1) GOTO 1000 
!.....   WRITE 3-D VELOCITIES.                                          
      WRITE (fLST, 925) (L, VMAG (L), L = 1, NE) 
      WRITE (fLST, 950) (L, VANG1 (L), L = 1, NE) 
      WRITE (fLST, 951) (L, VANG2 (L), L = 1, NE) 
  925 FORMAT(///11X,'S  T  E  A  D  Y  -  S  T  A  T  E     ',          &
     &   'F  L  U  I  D     V  E  L  O  C  I  T  Y'//                   &
     &   11X,'M A G N I T U D E   AT CENTROID OF ELEMENT'//             &
     &   2X,5(3X,'ELEMENT',16X)/(2X,5(1X,I9,1X,1PD15.8)))               
  950 FORMAT(///11X,'S  T  E  A  D  Y  -  S  T  A  T  E     ',          &
     &   'F  L  U  I  D     V  E  L  O  C  I  T  Y'//                   &
     &   11X,'A N G L E 1   AT CENTROID OF ELEMENT, IN DEGREES FROM ',  &
     &   '+X-AXIS TO PROJECTION OF FLOW DIRECTION IN XY-PLANE'//        &
     &   2X,5(3X,'ELEMENT',16X)/(2X,5(1X,I9,1X,1PD15.8)))               
  951 FORMAT(///11X,'S  T  E  A  D  Y  -  S  T  A  T  E     ',          &
     &   'F  L  U  I  D     V  E  L  O  C  I  T  Y'//                   &
     &   11X,'A N G L E 2   AT CENTROID OF ELEMENT, IN DEGREES FROM ',  &
     &   'XY-PLANE TO FLOW DIRECTION'//                                 &
     &   2X,5(3X,'ELEMENT',16X)/(2X,5(1X,I9,1X,1PD15.8)))               
!                                                                       
 1000 RETURN 
!                                                                       
      END SUBROUTINE OUTPTCLST_3D
