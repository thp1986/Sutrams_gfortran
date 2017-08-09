!     SUBROUTINE        R  O  T  A  T  E       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO TRANSFORM THE COORDINATES OF A VECTOR, {x}, BY APPLYING THE   
! ***  ROTATION MATRIX, [G]:  {xp}=[G]{x}.                              
!                                                                       
      SUBROUTINE ROTATE (G11, G12, G13, G21, G22, G23, G31, G32, G33, &
                         X, Y, Z, XP, YP, ZP)                                                 
      
      USE SutraMSPrecision

      IMPLICIT NONE
      !IMPLICIT DOUBLEPRECISION (A - H, O - Z) 

      REAL (DP) :: &
        G11, G12, G13, G21, G22, G23, G31, G32, G33, &
        X, Y, Z, &
        XP, YP, ZP
                                                                      
      XP = G11 * X + G12 * Y + G13 * Z 
      YP = G21 * X + G22 * Y + G23 * Z 
      ZP = G31 * X + G32 * Y + G33 * Z 
                                                                       
      RETURN 
      END SUBROUTINE ROTATE                         
