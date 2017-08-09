!     SUBROUTINE        B  A  S  I  S  2       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO CALCULATE VALUES OF BASIS AND WEIGHTING FUNCTIONS AND THEIR   
! ***  DERIVATIVES, TRANSFORMATION MATRICES BETWEEN LOCAL AND GLOBAL    
! ***  COORDINATES AND PARAMETER VALUES AT A SPECIFIED POINT IN A       
! ***  QUADRILATERAL FINITE ELEMENT.                                    
! ***  THIS SUBROUTINE HANDLES 2-D CALCULATIONS ONLY.                   
!                                                                       
      SUBROUTINE BASIS2 (ICALL, L, XLOC, YLOC, F, W, DET,     &
                         DFDXG, DFDYG, DWDXG, DWDYG,       &
                         THICKG, VXG, VYG, SWG, RHOG, VISCG, PORG, VGMAG, RELKG, &
                         CJ11, CJ12, CJ21, CJ22, RGXG, RGYG)                                         
      USE PARAMS 
      USE CONTRL 
      USE DIMS
      use SutraZoneModule
      use SutraStorage, ONLY : IN, X, Y, PITER, UITER, PVEL, THICK=>Z, &
                               GXSI, GETA, RCIT, RCITM1
      use SutraMSPrecision
      implicit none
!                                                                       
      integer (I4B) :: &
        ICALL, L
      real (DP) :: &
        XLOC, YLOC 
      real (DP) :: &
        UITERG (NSPE) 
      real (DP) :: &
        F (4), W (4), DFDXG (4), DFDYG (4), DWDXG (4), DWDYG (4) 
      real (DP) :: &
        FX (4), FY (4), AFX (4), AFY (4), DFDXL (4), DFDYL (4), &
        DWDXL (4), DWDYL (4), XDW (4), YDW (4), XIIX (4), YIIY (4)        
      real (DP) :: &
        DET, &
        THICKG, VXG, VYG, SWG, RHOG, VISCG, PORG, VGMAG, RELKG, &
        CJ11, CJ12, CJ21, CJ22, &
        RGXG, RGYG
      !locals
      integer (I4B) :: &
        I, II, IL, &
        K
      real (DP) :: &
        XF1, XF2, YF1, YF2, &
        ODET, &
        CIJ11, CIJ12, CIJ21, CIJ22, &
        RGXL, RGYL, RGXLM1, RGYLM1, &
        ADFDXL, ADFDYL, &
        RGXGM1, RGYGM1, &
        PITERG, DPDXG, DPDYG, &
        TEQN, &
        DSWDPG, &
        DENOM, &
        PGX, PGY, &
        VXG2, VYG2, VXL, VYL, VLMAG, &
        AA, BB, &
        XIXI, YIYI, THAAX, THBBY


      DATA XIIX / - 1.D0, + 1.D0, + 1.D0, - 1.D0 /, &
           YIIY / - 1.D0, - 1.D0, + 1.D0, + 1.D0 /                                          
      SAVE XIIX, YIIY 
!                                                                       
!                                                                       
!.....AT THIS LOCATION IN LOCAL COORDINATES, (XLOC,YLOC),               
!        CALCULATE SYMMETRIC WEIGHTING FUNCTIONS, F(I),                 
!        SPACE DERIVATIVES, DFDXG(I) AND DFDYG(I), AND                  
!        DETERMINANT OF JACOBIAN, DET.                                  
!                                                                       
      XF1 = 1.D0 - XLOC 
      XF2 = 1.D0 + XLOC 
      YF1 = 1.D0 - YLOC 
      YF2 = 1.D0 + YLOC 
!                                                                       
!.....CALCULATE BASIS FUNCTION, F.                                      
      FX (1) = XF1 
      FX (2) = XF2 
      FX (3) = XF2 
      FX (4) = XF1 
      FY (1) = YF1 
      FY (2) = YF1 
      FY (3) = YF2 
      FY (4) = YF2 
      DO I = 1, 4 
         F (I) = 0.250D0 * FX (I) * FY (I) 
  	  END DO
!                                                                       
!.....CALCULATE DERIVATIVES WITH RESPECT TO LOCAL COORDINATES.          
      DO I = 1, 4 
         DFDXL (I) = XIIX (I) * 0.250D0 * FY (I) 
         DFDYL (I) = YIIY (I) * 0.250D0 * FX (I) 
	    END DO
!                                                                       
!.....CALCULATE ELEMENTS OF JACOBIAN MATRIX, CJ.                        
      CJ11 = 0.D0 
      CJ12 = 0.D0 
      CJ21 = 0.D0 
      CJ22 = 0.D0 
      DO IL = 1, 4 
         II = (L - 1) * 4 + IL 
         I = IN (II) 
         CJ11 = CJ11 + DFDXL (IL) * X (I) 
         CJ12 = CJ12 + DFDXL (IL) * Y (I) 
         CJ21 = CJ21 + DFDYL (IL) * X (I) 
         CJ22 = CJ22 + DFDYL (IL) * Y (I) 
	    END DO
!                                                                       
!.....CALCULATE DETERMINANT OF JACOBIAN MATRIX.                         
      DET = CJ11 * CJ22 - CJ21 * CJ12 
!                                                                       
!.....RETURN TO ELEMEN2 WITH JACOBIAN MATRIX ON FIRST TIME STEP.        
      IF (ICALL.EQ.0) RETURN 
!                                                                       
!.....CALCULATE ELEMENTS OF INVERSE JACOBIAN MATRIX, CIJ.               
      ODET = 1.D0 / DET 
      CIJ11 = + ODET * CJ22 
      CIJ12 = - ODET * CJ12 
      CIJ21 = - ODET * CJ21 
      CIJ22 = + ODET * CJ11 
!                                                                       
!.....CALCULATE DERIVATIVES WITH RESPECT TO GLOBAL COORDINATES          
      DO I = 1, 4 
         DFDXG (I) = CIJ11 * DFDXL (I) + CIJ12 * DFDYL (I) 
         DFDYG (I) = CIJ21 * DFDXL (I) + CIJ22 * DFDYL (I) 
      END DO
!                                                                       
!.....CALCULATE CONSISTENT COMPONENTS OF (RHO*GRAV) TERM IN LOCAL       
!        COORDINATES AT THIS LOCATION, (XLOC,YLOC)                      
      RGXL = 0.D0 
      RGYL = 0.D0 
      RGXLM1 = 0.D0 
      RGYLM1 = 0.D0 
      DO IL = 1, 4 
         II = (L - 1) * 4 + IL 
         I = IN (II) 
         ADFDXL = DABS (DFDXL (IL) ) 
         ADFDYL = DABS (DFDYL (IL) ) 
         RGXL = RGXL + RCIT (I) * GXSI (L, IL) * ADFDXL 
         RGYL = RGYL + RCIT (I) * GETA (L, IL) * ADFDYL 
         RGXLM1 = RGXLM1 + RCITM1 (I) * GXSI (L, IL) * ADFDXL 
         RGYLM1 = RGYLM1 + RCITM1 (I) * GETA (L, IL) * ADFDYL 
      END DO 
!                                                                       
!.....TRANSFORM CONSISTENT COMPONENTS OF (RHO*GRAV) TERM TO             
!        GLOBAL COORDINATES                                             
      RGXG = CIJ11 * RGXL + CIJ12 * RGYL 
      RGYG = CIJ21 * RGXL + CIJ22 * RGYL 
      RGXGM1 = CIJ11 * RGXLM1 + CIJ12 * RGYLM1 
      RGYGM1 = CIJ21 * RGXLM1 + CIJ22 * RGYLM1 
!                                                                       
!.....CALCULATE PARAMETER VALUES AT THIS LOCATION, (XLOC,YLOC)          
!                                                                       
      PITERG = 0.D0 
      UITERG = 0.D0 
      DPDXG = 0.D0 
      DPDYG = 0.D0 
      PORG = 0.D0 
      THICKG = 0.0D0 
      DO IL = 1, 4 
         II = (L - 1) * 4 + IL 
         I = IN (II) 
         DPDXG = DPDXG + PVEL (I) * DFDXG (IL) 
         DPDYG = DPDYG + PVEL (I) * DFDYG (IL) 
         PORG = PORG + NodeData(NodeMap(i))%por * F (IL) 
         THICKG = THICKG + THICK (I) * F (IL) 
         PITERG = PITERG + PITER (I) * F (IL) 
         DO K = 1, NSPE 
           UITERG (K) = UITERG (K) + UITER (I, K) * F (IL) 
         END DO
      END DO 
!                                                                       
!.....SET VALUES FOR DENSITY AND VISCOSITY                              
!.....RHOG = FUNCTION(UITER)                                            
      RHOG = RHOW0 
      DO K = 1, NSPE 
        RHOG = RHOG + DRWDU (K) * (UITERG (K) - URHOW0 (K) ) 
	    END DO
!.....VISCG = FUNCTION(UITER)                                           
!        VISCOSITY IN UNITS OF VISC0*(KG/(M*SEC))                       
      VISCG = 0.0D0 
      TEQN = 0.0D0 
      IF (NESP.GT.0) TEQN = VISC0 (NESP) * 239.4D-7 * (10.D0** (248.37D0 / (UITERG (NESP) + 133.15D0) ) )                         
      IF (ME) 1300, 1300, 1200 
 1200 VISCG = TEQN 
      GOTO 1400 
!.....FOR SOLUTE TRANSPORT... AND                                       
!.....FOR SOLUTE AND ENERGY TRANSPORT                                   
!.....VISCG IS A FUNCTION OF CONCENTRATION AND/OR TEMPERATURE           
 1300 VISCG = BVISC0 + TEQN 
      DO K = 1, NSPE 
         IF (K.EQ.NESP) CYCLE
         VISCG = VISCG + VISC0 (K) * (UITERG (K) - URHOW0 (K) ) 
      END DO 
 1400 CONTINUE 
!                                                                       
!.....SET UNSATURATED FLOW PARAMETERS SWG AND RELKG                     
      IF (IUNSAT - 2) 1600, 1500, 1600 
 1500 IF (PITERG) 1550, 1600, 1600 
 1550 CALL UNSAT (SWG, DSWDPG, RELKG, PITERG, ElemMap(L) ) 
      GOTO 1700 
 1600 SWG = 1.0D0 
      RELKG = 1.0D0 
 1700 CONTINUE 
!                                                                       
!.....CALCULATE CONSISTENT FLUID VELOCITIES WITH RESPECT TO GLOBAL      
!        COORDINATES, VXG, VYG, AND VGMAG, AT THIS LOCATION, (XLOC,YLOC)
      DENOM = 1.D0 / (PORG * SWG * VISCG) 
      PGX = DPDXG - RGXGM1 
      PGY = DPDYG - RGYGM1 
!.....ZERO OUT RANDOM BOUYANT DRIVING FORCES DUE TO DIFFERENCING        
!..... NUMBERS PAST PRECISION LIMIT                                     
!..... MINIMUM DRIVING FORCE IS 1.D-10 OF PRESSURE GRADIENT             
!..... (THIS VALUE MAY BE CHANGED DEPENDING ON MACHINE PRECISION)       
      IF (DPDXG) 1720, 1730, 1720 
 1720 IF (DABS (PGX / DPDXG) - 1.0D-10) 1725, 1725, 1730 
 1725 PGX = 0.0D0 
 1730 IF (DPDYG) 1750, 1760, 1750 
 1750 IF (DABS (PGY / DPDYG) - 1.0D-10) 1755, 1755, 1760 
 1755 PGY = 0.0D0 
 1760 VXG = - DENOM * (ElemData(ElemMap(L))%permxx * PGX + ElemData(ElemMap(L))%permxy * PGY) * RELKG 
      VYG = - DENOM * (ElemData(ElemMap(L))%permyx * PGX + ElemData(ElemMap(L))%permyy * PGY) * RELKG 
      VXG2 = VXG * VXG 
      VYG2 = VYG * VYG 
      VGMAG = DSQRT (VXG2 + VYG2) 
!                                                                       
!.....AT THIS POINT IN LOCAL COORDINATES, (XLOC,YLOC),                  
!        CALCULATE ASYMMETRIC WEIGHTING FUNCTIONS, W(I),                
!        AND SPACE DERIVATIVES, DWDXG(I) AND DWDYG(I).                  
!                                                                       
!.....ASYMMETRIC FUNCTIONS SIMPLIFY WHEN  UP=0.0                        
      IF (UP.GT.1.0D-6.AND.NOUMAT.EQ.0) GOTO 1790 
      DO 1780 I = 1, 4 
         W (I) = F (I) 
         DWDXG (I) = DFDXG (I) 
         DWDYG (I) = DFDYG (I) 
 1780 END DO 
!.....RETURN WHEN ONLY SYMMETRIC WEIGHTING FUNCTIONS ARE USED           
      RETURN 
!                                                                       
!.....CALCULATE FLUID VELOCITIES WITH RESPECT TO LOCAL COORDINATES,     
!..... VXL, VYL, AND VLMAG, AT THIS LOCATION, (XLOC,YLOC).              
 1790 VXL = CIJ11 * VXG + CIJ21 * VYG 
      VYL = CIJ12 * VXG + CIJ22 * VYG 
      VLMAG = DSQRT (VXL * VXL + VYL * VYL) 
!                                                                       
      AA = 0.0D0 
      BB = 0.0D0 
      IF (VLMAG) 1900, 1900, 1800 
 1800 AA = UP * VXL / VLMAG 
      BB = UP * VYL / VLMAG 
!                                                                       
 1900 XIXI = .750D0 * AA * XF1 * XF2 
      YIYI = .750D0 * BB * YF1 * YF2 
      DO I = 1, 4 
         AFX (I) = .50D0 * FX (I) + XIIX (I) * XIXI 
         AFY (I) = .50D0 * FY (I) + YIIY (I) * YIYI 
	    END DO
!                                                                       
!.....CALCULATE ASYMMETRIC WEIGHTING FUNCTION, W.                       
      DO I = 1, 4 
         W (I) = AFX (I) * AFY (I) 
      END DO
!                                                                       
      THAAX = 0.50D0 - 1.50D0 * AA * XLOC 
      THBBY = 0.50D0 - 1.50D0 * BB * YLOC 
      DO I = 1, 4 
         XDW (I) = XIIX (I) * THAAX 
		    YDW (I) = YIIY (I) * THBBY 
	    END DO
!                                                                       
!.....CALCULATE DERIVATIVES WITH RESPECT TO LOCAL COORDINATES.          
      DO I = 1, 4 
         DWDXL (I) = XDW (I) * AFY (I) 
         DWDYL (I) = YDW (I) * AFX (I) 
	    END DO
!                                                                       
!.....CALCULATE DERIVATIVES WITH RESPECT TO GLOBAL COORDINATES.         
      DO I = 1, 4 
         DWDXG (I) = CIJ11 * DWDXL (I) + CIJ12 * DWDYL (I) 
         DWDYG (I) = CIJ21 * DWDXL (I) + CIJ22 * DWDYL (I) 
	    END DO
!                                                                       
!                                                                       
      RETURN 
      END SUBROUTINE BASIS2                         
