!     SUBROUTINE        B  A  S  I  S  3       SUTRA-MS VERSION 2004.1 
!                                                                       
! *** PURPOSE :                                                         
! ***  TO CALCULATE VALUES OF BASIS AND WEIGHTING FUNCTIONS AND THEIR   
! ***  DERIVATIVES, TRANSFORMATION MATRICES BETWEEN LOCAL AND GLOBAL    
! ***  COORDINATES AND PARAMETER VALUES AT A SPECIFIED POINT IN A       
! ***  QUADRILATERAL FINITE ELEMENT.                                    
! ***  THIS SUBROUTINE HANDLES 3-D CALCULATIONS ONLY.                   
!                                                                       
      SUBROUTINE BASIS3 (ICALL, L, XLOC, YLOC, ZLOC, F, W, &
                         DET, DFDXG, DFDYG, DFDZG, DWDXG, DWDYG, DWDZG, &
                         VXG, VYG, VZG, SWG, RHOG, VISCG, PORG, VGMAG, RELKG, &
                         CJ11, CJ12, CJ13, CJ21, CJ22, CJ23, CJ31, CJ32, CJ33, &
                         RGXG, RGYG, RGZG)                       
      USE PARAMS 
      USE CONTRL 
      USE DIMS
      use SutraZoneModule
      use SutraStorage, ONLY : IN, X, Y, Z, PITER, UITER, PVEL, &
                               GXSI, GETA, GZET, RCIT, RCITM1
      use SutraMSPrecision
      implicit none
                                                                       
      integer (I4B) :: &
        ICALL, L
      real (DP) :: &
        XLOC, YLOC, ZLOC 
      real (DP) :: &
        F (8), DFDXG (8), DFDYG (8), DFDZG (8) 
      real (DP) :: &
        W (8), DWDXG (8), DWDYG (8), DWDZG (8) 
      real (DP) :: &
        FX (8), FY (8), FZ (8), AFX (8), AFY (8), AFZ (8),      &
        DFDXL (8), DFDYL (8), DFDZL (8), DWDXL (8), DWDYL (8), DWDZL (8), &
        XDW (8), YDW (8), ZDW (8), XIIX (8), YIIY (8), ZIIZ (8)           
      real (DP) :: &
        UITERG (NSPE) 
      real (DP) :: &
        DET, &
        VXG, VYG, VZG, SWG, RHOG, VISCG, PORG, VGMAG, RELKG, &
        CJ11, CJ12, CJ13, CJ21, CJ22, CJ23, CJ31, CJ32, CJ33, &
        RGXG, RGYG, RGZG
      !locals
      integer (I4B) :: &
        I, II, IL, &
        K
      real (DP) :: &
        XF1, XF2, YF1, YF2, ZF1, ZF2, &
        ODET, &
        CIJ11, CIJ12, CIJ13, CIJ21, CIJ22, CIJ23, CIJ31, CIJ32, CIJ33, &
        RGXL, RGYL, RGZL, &
        RGXLM1, RGYLM1, RGZLM1, &
        ADFDXL, ADFDYL, ADFDZL, &
        RGXGM1, RGYGM1, RGZGM1, &
        PITERG, DPDXG, DPDYG, DPDZG, &
        TEQN, &
        DSWDPG, &
        DENOM, &
        PGX, PGY, PGZ, &
        VXG2, VYG2, VZG2, VXL, VYL, VZL, VLMAG, &
        AA, BB, GG, &
        XIXI, YIYI, ZIZI, &
        THAAX, THBBY, THGGZ

      DATA XIIX / - 1.D0, + 1.D0, + 1.D0, - 1.D0, - 1.D0, + 1.D0, + 1.D0, - 1.D0 /                                                  
      DATA YIIY / - 1.D0, - 1.D0, + 1.D0, + 1.D0, - 1.D0, - 1.D0, + 1.D0, + 1.D0 /                                                  
      DATA ZIIZ / - 1.D0, - 1.D0, - 1.D0, - 1.D0, + 1.D0, + 1.D0, + 1.D0, + 1.D0 /                                                  
      SAVE XIIX, YIIY, ZIIZ 
!                                                                       
!                                                                       
!                                                                       
!.....AT THIS LOCATION IN LOCAL COORDINATES, (XLOC,YLOC,ZLOC),          
!        CALCULATE SYMMETRIC WEIGHTING FUNCTIONS, F(I),                 
!        SPACE DERIVATIVES, DFDXG(I), DFDYX(I), AND DFDZG(I),           
!        AND DETERMINANT OF JACOBIAN, DET.                              
!                                                                       
!     CALCULATE BASIS/SYMMETRIC WEIGHTING FUNCTIONS.                    
!                                                                       
      XF1 = 1.D0 - XLOC 
      XF2 = 1.D0 + XLOC 
      YF1 = 1.D0 - YLOC 
      YF2 = 1.D0 + YLOC 
      ZF1 = 1.D0 - ZLOC 
      ZF2 = 1.D0 + ZLOC 
!                                                                       
      FX (1) = XF1 
      FX (2) = XF2 
      FX (3) = XF2 
      FX (4) = XF1 
      FX (5) = XF1 
      FX (6) = XF2 
      FX (7) = XF2 
      FX (8) = XF1 
      FY (1) = YF1 
      FY (2) = YF1 
      FY (3) = YF2 
      FY (4) = YF2 
      FY (5) = YF1 
      FY (6) = YF1 
      FY (7) = YF2 
      FY (8) = YF2 
      FZ (1) = ZF1 
      FZ (2) = ZF1 
      FZ (3) = ZF1 
      FZ (4) = ZF1 
      FZ (5) = ZF2 
      FZ (6) = ZF2 
      FZ (7) = ZF2 
      FZ (8) = ZF2 
      DO I = 1, 8 
        F (I) = 0.125D0 * FX (I) * FY (I) * FZ (I) 
      END DO
!                                                                       
      DO I = 1, 8 
         DFDXL (I) = XIIX (I) * 0.125D0 * FY (I) * FZ (I) 
         DFDYL (I) = YIIY (I) * 0.125D0 * FX (I) * FZ (I) 
         DFDZL (I) = ZIIZ (I) * 0.125D0 * FX (I) * FY (I) 
      END DO
!                                                                       
!     CALCULATE JACOBIAN, CJ(), AT CURRENT LOCAL COORDINATES.           
      CJ11 = 0.D0 
      CJ12 = 0.D0 
      CJ13 = 0.D0 
      CJ21 = 0.D0 
      CJ22 = 0.D0 
      CJ23 = 0.D0 
      CJ31 = 0.D0 
      CJ32 = 0.D0 
      CJ33 = 0.D0 
      DO IL = 1, 8 
         II = (L - 1) * 8 + IL 
         I = IN (II) 
         CJ11 = CJ11 + DFDXL (IL) * X (I) 
         CJ12 = CJ12 + DFDXL (IL) * Y (I) 
         CJ13 = CJ13 + DFDXL (IL) * Z (I) 
         CJ21 = CJ21 + DFDYL (IL) * X (I) 
         CJ22 = CJ22 + DFDYL (IL) * Y (I) 
         CJ23 = CJ23 + DFDYL (IL) * Z (I) 
         CJ31 = CJ31 + DFDZL (IL) * X (I) 
         CJ32 = CJ32 + DFDZL (IL) * Y (I) 
         CJ33 = CJ33 + DFDZL (IL) * Z (I)
      END DO
!                                                                       
!.....CALCULATE DETERMINANT OF JACOBIAN MATRIX.                         
      DET = CJ11 * (CJ22 * CJ33 - CJ32 * CJ23) - CJ21 * (CJ12 * CJ33 -  CJ32 * CJ13) + &
            CJ31 * (CJ12 * CJ23 - CJ22 * CJ13)                 
!                                                                       
!.....RETURN TO ELEMEN3 WITH JACOBIAN MATRIX ON FIRST TIME STEP.        
      IF (ICALL.EQ.0) RETURN 
!                                                                       
!                                                                       
!.....CALCULATE ELEMENTS OF INVERSE JACOBIAN MATRIX, CIJ.               
      ODET = 1.D0 / DET 
      CIJ11 = + ODET * (CJ22 * CJ33 - CJ32 * CJ23) 
      CIJ12 = - ODET * (CJ12 * CJ33 - CJ32 * CJ13) 
      CIJ13 = + ODET * (CJ12 * CJ23 - CJ22 * CJ13) 
      CIJ21 = - ODET * (CJ21 * CJ33 - CJ31 * CJ23) 
      CIJ22 = + ODET * (CJ11 * CJ33 - CJ31 * CJ13) 
      CIJ23 = - ODET * (CJ11 * CJ23 - CJ21 * CJ13) 
      CIJ31 = + ODET * (CJ21 * CJ32 - CJ31 * CJ22) 
      CIJ32 = - ODET * (CJ11 * CJ32 - CJ31 * CJ12) 
      CIJ33 = + ODET * (CJ11 * CJ22 - CJ21 * CJ12) 
!                                                                       
!.....CALCULATE DERIVATIVES WITH RESPECT TO GLOBAL COORDINATES          
      DO I = 1, 8 
         DFDXG (I) = CIJ11 * DFDXL (I) + CIJ12 * DFDYL (I) + CIJ13 * DFDZL (I)                                                      
         DFDYG (I) = CIJ21 * DFDXL (I) + CIJ22 * DFDYL (I) + CIJ23 * DFDZL (I)                                                      
         DFDZG (I) = CIJ31 * DFDXL (I) + CIJ32 * DFDYL (I) + CIJ33 * DFDZL (I)
      END DO
!                                                                       
!.....CALCULATE CONSISTENT COMPONENTS OF (RHO*GRAV) TERM IN LOCAL       
!        COORDINATES AT THIS LOCATION, (XLOC,YLOC,ZLOC)                 
      RGXL = 0.D0 
      RGYL = 0.D0 
      RGZL = 0.D0 
      RGXLM1 = 0.D0 
      RGYLM1 = 0.D0 
      RGZLM1 = 0.D0 
      DO IL = 1, 8 
         II = (L - 1) * 8 + IL 
         I = IN (II) 
         ADFDXL = DABS (DFDXL (IL) ) 
         ADFDYL = DABS (DFDYL (IL) ) 
         ADFDZL = DABS (DFDZL (IL) ) 
         RGXL = RGXL + RCIT (I) * GXSI (L, IL) * ADFDXL 
         RGYL = RGYL + RCIT (I) * GETA (L, IL) * ADFDYL 
         RGZL = RGZL + RCIT (I) * GZET (L, IL) * ADFDZL 
         RGXLM1 = RGXLM1 + RCITM1 (I) * GXSI (L, IL) * ADFDXL 
         RGYLM1 = RGYLM1 + RCITM1 (I) * GETA (L, IL) * ADFDYL 
         RGZLM1 = RGZLM1 + RCITM1 (I) * GZET (L, IL) * ADFDZL 
      END DO 
!                                                                       
!.....TRANSFORM CONSISTENT COMPONENTS OF (RHO*GRAV) TERM TO             
!        GLOBAL COORDINATES                                             
      RGXG = CIJ11 * RGXL + CIJ12 * RGYL + CIJ13 * RGZL 
      RGYG = CIJ21 * RGXL + CIJ22 * RGYL + CIJ23 * RGZL 
      RGZG = CIJ31 * RGXL + CIJ32 * RGYL + CIJ33 * RGZL 
      RGXGM1 = CIJ11 * RGXLM1 + CIJ12 * RGYLM1 + CIJ13 * RGZLM1 
      RGYGM1 = CIJ21 * RGXLM1 + CIJ22 * RGYLM1 + CIJ23 * RGZLM1 
      RGZGM1 = CIJ31 * RGXLM1 + CIJ32 * RGYLM1 + CIJ33 * RGZLM1 
!                                                                       
!.....CALCULATE PARAMETER VALUES AT THIS LOCATION, (XLOC,YLOC,ZLOC)     
!                                                                       
      PITERG = 0.D0 
      UITERG = 0.D0 
      DPDXG = 0.D0 
      DPDYG = 0.D0 
      DPDZG = 0.D0 
      PORG = 0.D0 
      DO IL = 1, 8 
         II = (L - 1) * 8 + IL 
         I = IN (II) 
         DPDXG = DPDXG + PVEL (I) * DFDXG (IL) 
         DPDYG = DPDYG + PVEL (I) * DFDYG (IL) 
         DPDZG = DPDZG + PVEL (I) * DFDZG (IL) 
         PORG = PORG + NodeData(NodeMap(i))%por * F (IL) 
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
      IF (NESP.GT.0) TEQN = VISC0 (NESP) * 239.4D-7 * (10.D0**(248.37D0 / (UITERG (NESP) + 133.15D0) ) )                         
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
!        COORDINATES, VXG, VYG, VZG, AND VGMAG, AT THIS LOCATION,       
!        (XLOC,YLOC,ZLOC)                                               
      DENOM = 1.D0 / (PORG * SWG * VISCG) 
      PGX = DPDXG - RGXGM1 
      PGY = DPDYG - RGYGM1 
      PGZ = DPDZG - RGZGM1 
!.....ZERO OUT RANDOM BOUYANT DRIVING FORCES DUE TO DIFFERENCING        
!..... NUMBERS PAST PRECISION LIMIT                                     
!..... MINIMUM DRIVING FORCE IS 1.D-10 OF PRESSURE GRADIENT             
!..... (THIS VALUE MAY BE CHANGED DEPENDING ON MACHINE PRECISION)       
      IF (DPDXG) 1720, 1727, 1720 
 1720 IF (DABS (PGX / DPDXG) - 1.0D-10) 1725, 1725, 1727 
 1725 PGX = 0.0D0 
 1727 IF (DPDYG) 1730, 1737, 1730 
 1730 IF (DABS (PGY / DPDYG) - 1.0D-10) 1735, 1735, 1737 
 1735 PGY = 0.0D0 
 1737 IF (DPDZG) 1740, 1760, 1740 
 1740 IF (DABS (PGZ / DPDZG) - 1.0D-10) 1745, 1745, 1760 
 1745 PGZ = 0.0D0 
 1760 continue
      VXG = - DENOM * (ElemData(ElemMap(L))%permxx * PGX + ElemData(ElemMap(L))%permxy * PGY + &
ElemData(ElemMap(L))%permxz * PGZ) * RELKG                                                    
      VYG = - DENOM * (ElemData(ElemMap(L))%permyx * PGX + ElemData(ElemMap(L))%permyy * PGY + &
ElemData(ElemMap(L))%permyz * PGZ) * RELKG                                                    
      VZG = - DENOM * (ElemData(ElemMap(L))%permzx * PGX + ElemData(ElemMap(L))%permzy * PGY + &
ElemData(ElemMap(L))%permzz * PGZ) * RELKG                                                    
      VXG2 = VXG * VXG 
      VYG2 = VYG * VYG 
      VZG2 = VZG * VZG 
      VGMAG = DSQRT (VXG2 + VYG2 + VZG2) 
!                                                                       
!.....AT THIS POINT IN LOCAL COORDINATES, (XLOC,YLOC,ZLOC),             
!        CALCULATE ASYMMETRIC WEIGHTING FUNCTIONS, W(I),                
!        AND SPACE DERIVATIVES, DWDXG(I), DWDYG(I), AND DWDZG(I)        
!                                                                       
!.....ASYMMETRIC FUNCTIONS SIMPLIFY WHEN  UP=0.0                        
      IF (UP.GT.1.0D-6.AND.NOUMAT.EQ.0) GOTO 1790 
      DO 1780 I = 1, 8 
         W (I) = F (I) 
         DWDXG (I) = DFDXG (I) 
         DWDYG (I) = DFDYG (I) 
         DWDZG (I) = DFDZG (I) 
 1780 END DO 
!.....RETURN WHEN ONLY SYMMETRIC WEIGHTING FUNCTIONS ARE USED           
      RETURN 
!                                                                       
!.....CALCULATE FLUID VELOCITIES WITH RESPECT TO LOCAL COORDINATES,     
!..... VXL, VYL, VZL, AND VLMAG, AT THIS LOCATION, (XLOC,YLOC,ZLOC).    
 1790 VXL = CIJ11 * VXG + CIJ21 * VYG + CIJ31 * VZG 
      VYL = CIJ12 * VXG + CIJ22 * VYG + CIJ32 * VZG 
      VZL = CIJ13 * VXG + CIJ23 * VYG + CIJ33 * VZG 
      VLMAG = DSQRT (VXL * VXL + VYL * VYL + VZL * VZL) 
!                                                                       
      AA = 0.0D0 
      BB = 0.0D0 
      GG = 0.0D0 
      IF (VLMAG) 1900, 1900, 1800 
 1800 AA = UP * VXL / VLMAG 
      BB = UP * VYL / VLMAG 
      GG = UP * VZL / VLMAG 
!                                                                       
 1900 XIXI = .750D0 * AA * XF1 * XF2 
      YIYI = .750D0 * BB * YF1 * YF2 
      ZIZI = .750D0 * GG * ZF1 * ZF2 
      DO I = 1, 8 
         AFX (I) = .50D0 * FX (I) + XIIX (I) * XIXI 
         AFY (I) = .50D0 * FY (I) + YIIY (I) * YIYI 
         AFZ (I) = .50D0 * FZ (I) + ZIIZ (I) * ZIZI 
      END DO
!                                                                       
!.....CALCULATE ASYMMETRIC WEIGHTING FUNCTION, W.                       
      DO I = 1, 8 
        W (I) = AFX (I) * AFY (I) * AFZ (I) 
      END DO
!                                                                       
      THAAX = 0.50D0 - 1.50D0 * AA * XLOC 
      THBBY = 0.50D0 - 1.50D0 * BB * YLOC 
      THGGZ = 0.50D0 - 1.50D0 * GG * ZLOC 
      DO I = 1, 8 
         XDW (I) = XIIX (I) * THAAX 
         YDW (I) = YIIY (I) * THBBY 
         ZDW (I) = ZIIZ (I) * THGGZ 
      END DO
!                                                                       
!.....CALCULATE DERIVATIVES WITH RESPECT TO LOCAL COORDINATES.          
      DO I = 1, 8 
         DWDXL (I) = XDW (I) * AFY (I) * AFZ (I) 
         DWDYL (I) = YDW (I) * AFX (I) * AFZ (I) 
         DWDZL (I) = ZDW (I) * AFX (I) * AFY (I) 
      END DO
!                                                                       
!.....CALCULATE DERIVATIVES WITH RESPECT TO GLOBAL COORDINATES.         
      DO I = 1, 8 
         DWDXG (I) = CIJ11 * DWDXL (I) + CIJ12 * DWDYL (I) + CIJ13 * DWDZL (I)                                                      
         DWDYG (I) = CIJ21 * DWDXL (I) + CIJ22 * DWDYL (I) + CIJ23 * DWDZL (I)                                                      
         DWDZG (I) = CIJ31 * DWDXL (I) + CIJ32 * DWDYL (I) + CIJ33 * DWDZL (I)
      END DO
!                                                                       
!                                                                       
      RETURN 
      END SUBROUTINE BASIS3                         
