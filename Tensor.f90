!     SUBROUTINE        T  E  N  S  O  R       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO TRANSFORM A SYMMETRIC TENSOR BETWEEN TWO COORDINATE           
! ***  SYSTEMS.  "T" IS THE MATRIX EXPRESSED IN THE FIRST (INPUT)       
! ***  COORDINATE SYSTEM, AND "P" IS THE MATRIX EXPRESSED IN THE SECOND 
! ***  (OUTPUT) COORDINATE SYSTEM.  ANGLE1, ANGLE2, AND ANGLE3 ARE THE  
! ***  "YAW", "PITCH", AND "ROLL" ANGLES THAT ROTATE THE SECOND         
! ***  COORDINATE SYSTEM INTO THE FIRST.                                
!                                                                       
      SUBROUTINE TENSOR (T11, T12, T13, T21, T22, T23, T31, T32, T33,   &
                         ANGLE1, ANGLE2, ANGLE3, &
                         P11, P12, P13, P21, P22, P23, P31, P32, P33, MSTRUC)

      USE SutraMSPrecision

      IMPLICIT NONE

      INTEGER (I4B) :: &
        MSTRUC
      REAL (DP) :: &
        T11, T12, T13, T21, T22, T23, T31, T32, T33,   &
        ANGLE1, ANGLE2, ANGLE3, &
        P11, P12, P13, P21, P22, P23, P31, P32, P33
      
      !LOCAL VARIABLES
      REAL (DP) :: &
        G11, G12, G13, G21, G22, G23, G31, G32, G33, &
        PP11, PP12, PP13, PP21, PP22, PP23, PP31, PP32, PP33


!                                                                       
!.....COMPUTE TRANSFORMATION MATRIX.                                    
!                                                                       
      CALL ROTMAT (ANGLE1, ANGLE2, ANGLE3, &
                   G11, G12, G13, G21, G22, G23,&
                   G31, G32, G33)                                                    
!                                                                       
!.....COMPUTE TRANSFORMED TENSOR.                                       
!                                                                       
      IF (MSTRUC.EQ.0) THEN 
         PP11 = G11 * T11 + G12 * T21 + G13 * T31 
         PP12 = G11 * T12 + G12 * T22 + G13 * T32 
         PP13 = G11 * T13 + G12 * T23 + G13 * T33 
         PP21 = G21 * T11 + G22 * T21 + G23 * T31 
         PP22 = G21 * T12 + G22 * T22 + G23 * T32 
         PP23 = G21 * T13 + G22 * T23 + G23 * T33 
         PP31 = G31 * T11 + G32 * T21 + G33 * T31 
         PP32 = G31 * T12 + G32 * T22 + G33 * T32 
         PP33 = G31 * T13 + G32 * T23 + G33 * T33 
         P11 = PP11 * G11 + PP12 * G12 + PP13 * G13 
         P12 = PP11 * G21 + PP12 * G22 + PP13 * G23 
         P13 = PP11 * G31 + PP12 * G32 + PP13 * G33 
         P22 = PP21 * G21 + PP22 * G22 + PP23 * G23 
         P23 = PP21 * G31 + PP22 * G32 + PP23 * G33 
         P33 = PP31 * G31 + PP32 * G32 + PP33 * G33 
      ELSE 
         P11 = T11 * G11 * G11 + T22 * G12 * G12 + T33 * G13 * G13 
         P12 = T11 * G11 * G21 + T22 * G12 * G22 + T33 * G13 * G23 
         P13 = T11 * G11 * G31 + T22 * G12 * G32 + T33 * G13 * G33 
         P22 = T11 * G21 * G21 + T22 * G22 * G22 + T33 * G23 * G23 
         P23 = T11 * G21 * G31 + T22 * G22 * G32 + T33 * G23 * G33 
         P33 = T11 * G31 * G31 + T22 * G32 * G32 + T33 * G33 * G33 
      ENDIF 
      P21 = P12 
      P31 = P13 
      P32 = P23 
!                                                                       
      RETURN 
      END SUBROUTINE TENSOR                         
