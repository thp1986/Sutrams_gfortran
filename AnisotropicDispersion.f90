!     SUBROUTINE        A  N  I  D  S  P       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO COMPUTE LONGITUDINAL AND TRANSVERSE DISPERSIVITIES USING AN   
! ***  AD HOC, 3-D ANISOTROPIC DISPERSION MODEL.  THREE DISPERSIVITIES  
! ***  ARE COMPUTED:                                                    
! ***  (1) AL = LONGITUDUNAL DISPERSIVITY                               
! ***  (2) ATH = "HORIZONTAL" TRANSVERSE DISPERSIVITY                   
! ***  (3) ATV = "VERTICAL" TRANSVERSE DISPERSIVITY                     
! ***  ALSO, TO RETURN THE ANGLES OF THE VELOCITY WITH RESPECT TO THE   
! ***  PRINCIPAL PERMEABILITY DIRECTIONS.                               
!                                                                       
      SUBROUTINE ANIDSP (VX, VY, VZ, VMAG, ANG1, ANG2, ANG3, ALMAX, ALMID, ALMIN, &
	                       ATHMAX, ATHMID, ATHMIN, ATVMAX, ATVMID, ATVMIN, &
						             AL, ATH, ATV, THETA, PHI)                                             
        use SutraMSPrecision 
        USE PARAMS 
        IMPLICIT NONE
        REAL (DP) :: &
          VX, VY, VZ, VMAG, &
          ANG1, ANG2, ANG3, &
          ALMAX, ALMID, ALMIN, &
          ATHMAX, ATHMID, ATHMIN, &
          ATVMAX, ATVMID, ATVMIN, &
          AL, ATH, ATV, &
          THETA, PHI
        !LOCALS
          REAL (DP) :: &
            PI, PIH, &
            G11, G12, G13, &
            G21, G22, G23, &
            G31, G32, G33, &
            VXX, VYY, VZZ, VZL, &
            CPHI, SPHI, &
            CTHETA, STHETA, &
            FAC1, FAC2, FAC3, &
            DENOM
          INTEGER (I4B) :: &
            IL, &
            ITH, &
            ITV
                                                                       
        IF (VMAG.EQ.0D0) THEN 
           AL = 0D0 
           ATH = 0D0 
           ATV = 0D0 
           THETA = 0d0 
           PHI = 0d0 
           RETURN 
        ENDIF 
                                                                       
                                                                       
        PI = 3.1415926536D0 
        PIH = 5D-1 * PI 
!                                                                       
!.......COMPUTE ROTATION MATRIX AND ROTATE V TO COORDINATE SYSTEM         
!          ALIGNED WITH THE PERMEABILITY TENSOR.                          
!                                                                       
        CALL ROTMAT (ANG1, ANG2, ANG3, &
	                   G11, G12, G13, G21, G22, G23, G31, G32, G33)                                                         
        CALL ROTATE (G11, G21, G31, G12, G22, G32, G13, G23, G33, &
	                   VX, VY, VZ, VXX, VYY, VZZ)                                                
!                                                                       
!.......COMPUTE VELOCITY DIRECTION COSINES.                               
!                                                                       
        VZL = VZZ / VMAG 
!                                                                       
!.......COMPUTE INCLINATION ANGLE PHI FIRST.                              
!                                                                       
        PHI = ASIN (VZL) 
        CPHI = COS (PHI) 
        SPHI = SIN (PHI) 
!                                                                       
!.......NEXT, COMPUTE AZIMUTH ANGLE THETA.                                
        IF (VXX.NE.0.0D0) THEN 
           THETA = ATAN2 (VYY, VXX) 
           CTHETA = COS (THETA) 
           STHETA = SIN (THETA) 
        ELSE 
           THETA = 0.0D0 
           CTHETA = 0.0D0 
           STHETA = 0.0D0 
        ENDIF 

        FAC1 = CPHI * CPHI * CTHETA * CTHETA 
        FAC2 = CPHI * CPHI * STHETA * STHETA 
        FAC3 = SPHI * SPHI 
!                                                                       
        IL = 1 
        ITH = 1 
        ITV = 1 
        IF (ALMAX.EQ.ALMID.AND.ALMAX.EQ.ALMIN) THEN 
           IL = 0 
           AL = ALMAX * ATSPMULT (KSP) 
        ENDIF 
        IF (ATHMAX.EQ.ATHMID.AND.ATHMAX.EQ.ATHMIN) THEN 
           ITH = 0 
           ATH = ATHMAX * ATSPMULT (KSP) 
        ENDIF 
        IF (ATVMAX.EQ.ATVMID.AND.ATVMAX.EQ.ATVMIN) THEN 
           ITV = 0 
           ATV = ATVMAX * ATSPMULT (KSP) 
        ENDIF 

        IF (IL.EQ.0.AND.ITH.EQ.0.AND.ITV.EQ.0) RETURN 
!                                                                       
!.......COMPUTE VALUE OF AL.                                              
!                                                                       
        IF (IL.EQ.1) THEN 
           DENOM = (ALMID * ALMIN * FAC1 + ALMAX * ALMIN * FAC2 + ALMAX * ALMID * FAC3)                                                  
           AL = ALMAX * ALMID * ALMIN * ATSPMULT (KSP) / DENOM 
        ENDIF 
!                                                                       
!.......COMPUTE VALUE OF ATH.                                             
!                                                                       
        IF (ITH.EQ.1) THEN 
           DENOM = (ATHMID * ATHMIN * FAC1 + ATHMAX * ATHMIN * FAC2 + ATHMAX * ATHMID * FAC3)                                        
           ATH = ATHMAX * ATHMID * ATHMIN * ATSPMULT (KSP) / DENOM 
        ENDIF 
!                                                                       
!.......COMPUTE VALUE OF ATV.                                             
!                                                                       
        IF (ITV.EQ.1) THEN 
           DENOM = (ATVMID * ATVMIN * FAC1 + ATVMAX * ATVMIN * FAC2 + ATVMAX * ATVMID * FAC3)                                        
           ATV = ATVMAX * ATVMID * ATVMIN * ATSPMULT (KSP) / DENOM 
        ENDIF 
!                                                                       
!.......RETURN TO CALLING ROUTINE                                         
        RETURN 
      END SUBROUTINE ANIDSP                         
