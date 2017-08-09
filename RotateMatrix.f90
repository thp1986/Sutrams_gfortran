!                                                                       
!     SUBROUTINE        R  O  T  M  A  T       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO COMPUTE A TRANSFORMATION MATRIX, [G], THAT CONVERTS           
! ***  COORDINATES OF A VECTOR, {v}, FROM A COORDINATE SYSTEM (X, Y, Z) 
! ***  TO A NEW COORDINATE SYSTEM (X', Y', Z'):  {v'} = [G]{v}.         
! ***  THE OVERALL TRANSFORMATION IS THE RESULT OF THREE ROTATIONS      
! ***  APPLIED CONSECUTIVELY:                                           
! ***  A1 = ROTATION IN THE XY-PLANE, COUNTER-CLOCKWISE FROM THE        
! ***     +X-AXIS (LOOKING DOWN THE +Z-AXIS TOWARD THE ORIGIN),         
! ***  A2 = ROTATION IN THE NEW XZ-PLANE, COUNTER-CLOCKWISE FROM THE    
! ***     NEW +X-AXIS (LOOKING DOWN THE NEW +Y-AXIS TOWARD THE ORIGIN), 
! ***  A3 = ROTATION IN THE NEW YZ-PLANE, COUNTER-CLOCKWISE FROM THE    
! ***     NEW +Y-AXIS (LOOKING DOWN THE NEW +X-AXIS TOWARD THE ORIGIN). 
!                                                                       

      SUBROUTINE ROTMAT (A1, A2, A3, &
                         G11, G12, G13, G21, G22, G23, G31, G32, G33)

      USE SutraMSPrecision

      IMPLICIT NONE

      REAL (DP) :: &
        A1, A2, A3, &
        G11, G12, G13, G21, G22, G23, G31, G32, G33
      
      !LOCAL VARIABLES
      REAL (DP) :: &
        C1, C2, C3, &
        S1, S2, S3

      S1 = DSIN (A1) 
      C1 = DCOS (A1) 
      S2 = DSIN (A2) 
      C2 = DCOS (A2) 
      S3 = DSIN (A3) 
      C3 = DCOS (A3) 
!                                                                       
!.....COMPUTE ROTATION MATRIX.                                          
!                                                                       
      G11 = C1 * C2 
      G12 = - C1 * S2 * S3 - S1 * C3 
      G13 = - C1 * S2 * C3 + S1 * S3 
      G21 = S1 * C2 
      G22 = - S1 * S2 * S3 + C1 * C3 
      G23 = - S1 * S2 * C3 - C1 * S3 
      G31 = S2 
      G32 = C2 * S3 
      G33 = C2 * C3 
      RETURN 
      END SUBROUTINE ROTMAT                         
