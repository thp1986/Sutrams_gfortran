!     SUBROUTINE        E  L  E  M  N  3       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO CONTROL AND CARRY OUT ALL CALCULATIONS FOR EACH ELEMENT BY    
! ***  OBTAINING ELEMENT INFORMATION FROM THE BASIS FUNCTION ROUTINE,   
! ***  CARRYING OUT GAUSSITRIN INTEGRATION OF FINITE ELEMENT INTEGRALS,   
! ***  AND ASSEMBLING RESULTS OF ELEMENTWISE INTEGRATIONS INTO          
! ***  A GLOBAL MATRIX AND GLOBAL VECTOR FOR BOTH FLOW AND TRANSPORT    
! ***  EQUATIONS.                                                       
! ***  ALSO CALCULATES VELOCITY AT EACH ELEMENT CENTROID FOR            
! ***  PRINTED OUTPUT.                                                  
! ***  THIS SUBROUTINE HANDLES 3-D CALCULATIONS ONLY.                   
!                                                                       
!                                                                       
      SUBROUTINE ELEMN3 (ML)                                                            
      USE ITERAT 
      USE PARAMS 
      USE CONTRL 
      USE SOLVI 
      USE JCOLS 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE GRAVEC
      USE SutraZoneModule
      USE SutraStorage
      USE MSErrorHandler
      USE SutraMSPrecision
      USE TECPLOT
      USE ColumnStorage

      implicit none
      real (DP) :: &
        BFLOWE (8, 8), DFLOWE (8), BTRANE (8, 8), DTRANE (8, 8), VOLE (8)                                                          
      real (DP) :: &
        F (8, 8), W (8, 8), DET (8), DFDXG (8, 8), DFDYG (8, 8),&
        DFDZG (8, 8), DWDXG (8, 8), DWDYG (8, 8), DWDZG (8, 8)            
      real (DP) :: &
        SWG (8), RHOG (8), VISCG (8), PORG (8), VXG (8), VYG (8), VZG (8), &
        RELKG (8), RGXG (8), RGYG (8), RGZG (8),        &
        VGMAG (8)                                                         
      real (DP) :: &
        RXXG (8), RXYG (8), RXZG (8), RYXG (8), RYYG (8),       &
        RYZG (8), RZXG (8), RZYG (8), RZZG (8)                            
      real (DP) :: &
        BXXG (8), BXYG (8), BXZG (8), BYXG (8), BYYG (8),       &
        BYZG (8), BZXG (8), BZYG (8), BZZG (8), EXG (8), EYG (8), EZG (8) 
      real (DP) :: &
        GXLOC (8), GYLOC (8), GZLOC (8) 
      !locals
      integer (I4B) :: &
        INTIM,  ISTOP, &
        IVPRNT, KVPRNT, &
        I, II, IL, I1, IXL, IYL, IZL, &
        J, &
        KG, &
        L, &
        ML
      integer (I4B) :: &
        MSTRUC
      real (DP) :: &
        GLOC, XLOC, YLOC, ZLOC, &
        CJ11, CJ12, CJ13, CJ21, CJ22, CJ23, CJ31, CJ32, CJ33, &
        XIX, YIY, ZIZ, &
        AXSUM, AYSUM, AZSUM, &
        SWTEST, ROMG, &
        DF, VO, BF, &
        ESWG, RHOCWG, ESRCG, &
        DXXG, DXYG, DXZG, &
        DYXG, DYYG, DYZG, &
        DZXG, DZYG, DZZG, &
        ALEFF, AT1EFF, AT2EFF, &
        ANGKV1, ANGKV2, &
        DLG, DTHG, DTVG, &
        ZERO, ANGLE1, ANGLE2, ANGLE3, &
        DDXXG, DDXYG, DDXZG, &
        DDYXG, DDYYG, DDYZG, &
        DDZXG, DDZYG, DDZZG, &
        ESE, BT, DT

      !data
      DATA GLOC / 0.577350269189626D0 / 
      DATA INTIM / 0 /, ISTOP / 0 / 
      DATA GXLOC / - 1.D0, 1.D0, 1.D0, - 1.D0, - 1.D0, 1.D0, 1.D0, - 1.D0 /                                                          
      DATA GYLOC / - 1.D0, - 1.D0, 1.D0, 1.D0, - 1.D0, - 1.D0, 1.D0, 1.D0 /                                                            
      DATA GZLOC / - 1.D0, - 1.D0, - 1.D0, - 1.D0, 1.D0, 1.D0, 1.D0, 1.D0 /                                                            
      SAVE GLOC, INTIM, ISTOP, GXLOC, GYLOC, GZLOC 
!
!
      MSErrorValue%cDataSet='EL3'
!                                                                       
!.....DECIDE WHETHER TO CALCULATE CENTROID VELOCITIES ON THIS CALL      
      IVPRNT = 0 
      IF (MOD (IT, NPRINT) .EQ.0.AND.ML.NE.2.AND.IT.NE.0.AND.ITER.EQ.1) &
      IVPRNT = 1                                                        
      IF (MOD (IT, LCOLPR) .EQ.0.AND.ML.NE.2.AND.IT.NE.0.AND.ITER.EQ.1) &
      IVPRNT = 1                                                        
      IF (IT.EQ.1) IVPRNT = 1 
      IF ( (TSEC.GE.TMAX.OR.IT.EQ.ITMAX) .AND.ML.NE.2.AND.ITER.EQ.1)    &
      IVPRNT = 1                                                        
      KVPRNT = IVPRNT + KVEL 
!                                                                       
!.....ON FIRST TIME STEP, PREPARE GRAVITY VECTOR COMPONENTS,            
!        GXSI, GETA, AND GZET FOR CONSISTENT VELOCITIES,                
!        AND CHECK ELEMENT SHAPES                                       
      IF (INTIM) 100, 100, 2000 
  100 INTIM = 1 
!.....LOOP THROUGH ALL ELEMENTS TO OBTAIN THE JACOBIAN                  
!        AT EACH OF THE EIGHT NODES IN EACH ELEMENT                     
!                                                                       
!.....SCREEN OUTPUT
!      WRITE ( * , * ) 'CALCULATING 3-D JACOBIAN ON FIRST PASS' 
!.....SCREEN OUTPUT
!                                                                       
      DO 1000 L = 1, NE 
         DO 500 IL = 1, 8 
            XLOC = GXLOC (IL) 
            YLOC = GYLOC (IL) 
            ZLOC = GZLOC (IL) 
            CALL BASIS3 (0000, L, XLOC, YLOC, ZLOC, F (1,IL), &
                         W (1, IL), DET (IL), DFDXG (1, IL), DFDYG (1, IL),     &
                         DFDZG (1, IL), DWDXG (1, IL), DWDYG (1, IL), DWDZG (1, IL), &
                         VXG (IL), VYG (IL), VZG (IL),      &
                         SWG (IL), RHOG (IL), VISCG (IL), PORG (IL), VGMAG (IL),     &
                         RELKG (IL), &
                         CJ11, CJ12, CJ13, CJ21, CJ22, CJ23, &
                         CJ31, CJ32, CJ33, RGXG (IL), RGYG (IL), RGZG (IL))
            GXSI (L, IL) = CJ11 * GRAVX + CJ12 * GRAVY + CJ13 * GRAVZ 
            GETA (L, IL) = CJ21 * GRAVX + CJ22 * GRAVY + CJ23 * GRAVZ 
            GZET (L, IL) = CJ31 * GRAVX + CJ32 * GRAVY + CJ33 * GRAVZ 
!.....CHECK FOR NEGATIVE- OR ZERO-AREA ERRORS IN ELEMENT SHAPES         
            IF (DET (IL) ) 200, 200, 500 
  200       ISTOP = ISTOP + 1 
            
            II = (L - 1) * N48 + IL 
            I = IN (II)

            WRITE (fLST, 400) IN ( (L - 1) * 8 + IL), L, DET (IL) 
  400       FORMAT  (11X,'THE DETERMINANT OF THE JACOBIAN AT NODE ',I9,       &
                         ' IN ELEMENT ',I9,' IS NEGATIVE OR ZERO, ',1PE15.7)          
  500    END DO 
 1000 END DO 
!                                                                       
      IF (ISTOP.EQ.0) GOTO 2000 

      WRITE (fLST, 1500) 
 1500 FORMAT(//////11X,'SOME ELEMENTS HAVE INCORRECT GEOMETRY.'         &
     &   //11X,'PLEASE CHECK THE NODE COORDINATES AND ',                &
     &   'INCIDENCE LIST, MAKE CORRECTIONS, AND THEN RERUN.'////////    &
     &   11X,'S I M U L A T I O N   H A L T E D'/                       &
     &   11X,'___________________   ___________')                       
!
!
      call ErrorIO('ELEMN3: Some elements have incorrect geometry')
!                                                                       
!.....LOOP THROUGH ALL ELEMENTS TO CARRY OUT SPATITRIL INTEGRATION        
!        OF FLUX TERMS IN P AND/OR U EQUATIONS                          
 2000 IF (IUNSAT.NE.0) IUNSAT = 2 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SCREEN OUTPUT
      I1 = - 1 
!.....SCREEN OUTPUT
!                                                                       
      DO 9999 L = 1, NE 
         XIX = - 1.D0 
         YIY = - 1.D0 
         ZIZ = - 1.D0 
         KG = 0 
!.....OBTAIN BASIS FUNCTION AND RELATED INFORMATION AT EACH OF          
!        FOUR GAUSS POINTS IN THE ELEMENT                               
!                                                                       
!.....SCREEN OUTPUT
         I1 = I1 + 1 
         IF (I1.EQ.0) THEN 
          WRITE ( * ,  * ) 'CALCULATING 3-D BASIS AND WEIGHTING FUNCTIONS' 
         ENDIF 
!.....SCREEN OUTPUT
!                                                                       
         DO 2250 IZL = 1, 2 
            DO 2200 IYL = 1, 2 
               DO 2100 IXL = 1, 2 
                  KG = KG + 1 
                  XLOC = XIX * GLOC 
                  YLOC = YIY * GLOC 
                  ZLOC = ZIZ * GLOC 
                  CALL BASIS3 (0001, L, XLOC, YLOC, ZLOC,  &
                               F (1, KG), W (1, KG), DET (KG), DFDXG (1, KG),        &
                               DFDYG (1, KG), DFDZG (1, KG), DWDXG (1, KG), DWDYG (1,KG), &
                               DWDZG (1, KG), VXG (KG),&
                               VYG (KG), VZG (KG), SWG (KG), RHOG (KG), VISCG (KG),  &
                               PORG (KG), VGMAG (KG), RELKG (KG), &
                               CJ11, CJ12, CJ13, CJ21, CJ22, CJ23, CJ31,     &
                               CJ32, CJ33, RGXG (KG), RGYG (KG), RGZG (KG))
 2100          XIX = - XIX 
 2200       YIY = - YIY 
 2250    ZIZ = - ZIZ 
!.....VERSION 1.1
!.....SAVE NODAL VELOCITIES IF REQUIRED
      IF (LTpPrintVelocityAtNodes) THEN
        DO KG = 1, N48
         II = (L - 1) * N48 + KG 
         I = IN (II) 
         CALL MakeVelocityVector(I, VXG(KG), VYG(KG), VZG(KG))
        END DO
      END IF
!                                                                       
!.....CALCULATE VELOCITY AT ELEMENT CENTROID WHEN REQUIRED              
         IF (KVPRNT - 2) 3000, 2300, 3000 
 2300    AXSUM = 0.0D0 
         AYSUM = 0.0D0 
         AZSUM = 0.0D0 
         DO 2400 KG = 1, 8 
            AXSUM = AXSUM + VXG (KG) 
            AZSUM = AZSUM + VZG (KG) 
 2400    AYSUM = AYSUM + VYG (KG) 
         VMAG (L) = DSQRT (AXSUM * AXSUM + AYSUM * AYSUM + AZSUM *      &
         AZSUM)                                                         
         IF (VMAG (L) .NE.0D0) THEN 
            VANG2 (L) = DASIN (AZSUM / VMAG (L) ) *                     &
            5.729577951308232D+1                                        
            VMAG (L) = VMAG (L) * 1.25D-1 
            VANG1 (L) = DATAN2 (AYSUM, AXSUM) * 5.729577951308232D+1 
         ELSE 
            VANG1 (L) = 0D0 
            VANG2 (L) = 0D0 
         ENDIF 
!                                                                       
 3000    CONTINUE 
!                                                                       
!.....CALCULATE PARAMETERS FOR FLUID MASS BALANCE AT GAUSS POINTS       
         IF (ML - 1) 3400, 3400, 6100 
 3400    SWTEST = 0.D0 
!                                                                       
!.....SCREEN OUTPUT
!         IF (I1.EQ.0) THEN 
!            WRITE ( * , * ) 'CALCULATING 3-D PARAMETERS FOR FLUID MASS' 
!            WRITE ( * , * ) 'BALANCE AT GAUSS POINTS' 
!         ENDIF 
!.....SCREEN OUTPUT
!                                                                       
         DO 4000 KG = 1, 8 
            SWTEST = SWTEST + SWG (KG) 
            ROMG = RHOG (KG) * RELKG (KG) / VISCG (KG) 

            RXXG (KG) = ElemData(ElemMap(L))%permxx * ROMG 
            RXYG (KG) = ElemData(ElemMap(L))%permxy * ROMG 
            RXZG (KG) = ElemData(ElemMap(L))%permxz * ROMG 
            RYXG (KG) = ElemData(ElemMap(L))%permyx * ROMG 
            RYYG (KG) = ElemData(ElemMap(L))%permyy * ROMG 
            RYZG (KG) = ElemData(ElemMap(L))%permyz * ROMG 
            RZXG (KG) = ElemData(ElemMap(L))%permzx * ROMG 
            RZYG (KG) = ElemData(ElemMap(L))%permzy * ROMG 
            RZZG (KG) = ElemData(ElemMap(L))%permzz * ROMG 
 4000    END DO 
!                                                                       
!.....INTEGRATE FLUID MASS BALANCE IN AN UNSATURATED ELEMENT            
!        USING ASYMMETRIC WEIGHTING FUNCTIONS                           
         IF (UP.LE.1.0D-6) GOTO 5200 
         IF (SWTEST - 3.999D0) 4200, 5200, 5200 
 4200    DO 5000 I = 1, 8 
            DF = 0.D0 
            VO = 0.D0 
            DO 4400 KG = 1, 8 
               VO = VO + F (I, KG) * DET (KG) 
 4400       DF = DF + ( (RXXG (KG) * RGXG (KG) + RXYG (KG) * RGYG (KG) + &
                         RXZG (KG) * RGZG (KG) ) * DWDXG (I, KG) + &
                        (RYXG (KG) * RGXG (KG) + RYYG (KG) * RGYG (KG) + &
                         RYZG (KG) * RGZG (KG) ) * DWDYG (I, KG) + &
                        (RZXG (KG) * RGXG (KG) + RZYG (KG) * RGYG (KG) + &
                         RZZG (KG) * RGZG (KG) ) * DWDZG (I, KG) ) * DET (KG)                                                  
            DO 4800 J = 1, 8 
               BF = 0.D0 
               DO 4600 KG = 1, 8 
 4600          BF = BF + ( (RXXG (KG) * DFDXG (J, KG) + RXYG (KG) * DFDYG (J, KG) + &
                            RXZG (KG) * DFDZG (J, KG) ) * DWDXG (I,KG) + &
                           (RYXG (KG) * DFDXG (J, KG) + RYYG (KG) * DFDYG (J, KG) + &
                            RYZG (KG) * DFDZG (J, KG) ) * DWDYG (I, KG) + &
                           (RZXG (KG) * DFDXG (J, KG) + RZYG (KG) * DFDYG (J, KG) + &
                            RZZG (KG) * DFDZG (J, KG) ) * DWDZG (I, KG) ) * DET (KG)                                                      
 4800       BFLOWE (I, J) = BF 
            VOLE (I) = VO 
 5000    DFLOWE (I) = DF 
         GOTO 6200 
!                                                                       
!.....INTEGRATE FLUID MASS BALANCE IN A SATURATED OR UNSATURATED        
!        ELEMENT USING SYMMETRIC WEIGHTING FUNCTIONS                    
 5200    DO 6000 I = 1, 8 
            DF = 0.D0 
            VO = 0.D0 
            DO 5400 KG = 1, 8 
               VO = VO + F (I, KG) * DET (KG) 
 5400       DF = DF + ( (RXXG (KG) * RGXG (KG) + RXYG (KG) * RGYG (KG)    &
            + RXZG (KG) * RGZG (KG) ) * DFDXG (I, KG) + (RYXG (KG)        &
            * RGXG (KG) + RYYG (KG) * RGYG (KG) + RYZG (KG) * RGZG (KG) ) &
            * DFDYG (I, KG) + (RZXG (KG) * RGXG (KG) + RZYG (KG)          &
            * RGYG (KG) + RZZG (KG) * RGZG (KG) ) * DFDZG (I, KG) )       &
            * DET (KG)                                                  
            DO 5800 J = 1, 8 
               BF = 0.D0 
               DO 5600 KG = 1, 8 
 5600          BF = BF + ( (RXXG (KG) * DFDXG (J, KG) + RXYG (KG)           &
               * DFDYG (J, KG) + RXZG (KG) * DFDZG (J, KG) ) * DFDXG (I,KG) &
               + (RYXG (KG) * DFDXG (J, KG) + RYYG (KG) * DFDYG (J, KG)     &
               + RYZG (KG) * DFDZG (J, KG) ) * DFDYG (I, KG) + (RZXG (KG)   &
               * DFDXG (J, KG) + RZYG (KG) * DFDYG (J, KG) + RZZG (KG)      &
               * DFDZG (J, KG) ) * DFDZG (I, KG) ) * DET (KG)                                                      
 5800       BFLOWE (I, J) = BF 
            VOLE (I) = VO 
 6000    DFLOWE (I) = DF 
 6200    CONTINUE 
         IF (ML - 1) 6100, 9000, 6100 
 6100    IF (NOUMAT.EQ.1) GOTO 9000 
!                                                                       
!                                                                       
!.....CALCULATE PARAMETERS FOR ENERGY BALANCE OR SOLUTE MASS BALANCE    
!        AT GAUSS POINTS                                                
!                                                                       
!.....SCREEN OUTPUT
!         IF (I1.EQ.0) THEN 
!            WRITE ( * , * ) 'CALCULATING 3-D PARAMETERS FOR ENERGY OR' 
!            WRITE ( * , * ) 'SOLUTE MASS BALANCE AT GAUSS POINTS' 
!            WRITE ( * , * ) 'SPECIES [', KSP, ']' 
!         ENDIF 
!.....SCREEN OUTPUT
!                                                                       
         DO 7000 KG = 1, 8 
            ESWG = PORG (KG) * SWG (KG) 
            RHOCWG = RHOG (KG) * CW 
            ESRCG = ESWG * RHOCWG 
            IF (VGMAG (KG) ) 6300, 6300, 6600 
 6300       EXG (KG) = 0.0D0 
            EYG (KG) = 0.0D0 
            EZG (KG) = 0.0D0 
            DXXG = 0.0D0 
            DXYG = 0.0D0 
            DXZG = 0.0D0 
            DYXG = 0.0D0 
            DYYG = 0.0D0 
            DYZG = 0.0D0 
            DZXG = 0.0D0 
            DZYG = 0.0D0 
            DZZG = 0.0D0 
            GOTO 6900 
 6600       EXG (KG) = ESRCG * VXG (KG) 
            EYG (KG) = ESRCG * VYG (KG) 
            EZG (KG) = ESRCG * VZG (KG) 
!                                                                       
!                                                                       
!.....DISPERSIVITY MODEL FOR ANISOTROPIC MEDITRI                          
!        WITH PRINCIPAL DISPERSIVITIES:  ALMAX,ALMID,ALMIN,             
!        AT1MAX,AT1MID,AT1MIN,AT2MAX,AT2MID,AT2MIN.                     
            CALL ANIDSP (VXG (KG), VYG (KG), VZG (KG), VGMAG (KG),      &
              ElemData(ElemMap(L))%pangl1, ElemData(ElemMap(L))%pangl2, ElemData(ElemMap(L))%pangl3, &
              ElemData(ElemMap(L))%almax, ElemData(ElemMap(L))%almid, ElemData(ElemMap(L))%almin, &
              ElemData(ElemMap(L))%at1max, ElemData(ElemMap(L))%at1mid, ElemData(ElemMap(L))%at1min, &
              ElemData(ElemMap(L))%at2max, ElemData(ElemMap(L))%at2mid, ElemData(ElemMap(L))%at2min, &
              ALEFF, AT1EFF, AT2EFF, ANGKV1, ANGKV2)                                                     
            DLG = ALEFF * VGMAG (KG) 
            DTHG = AT1EFF * VGMAG (KG) 
            DTVG = AT2EFF * VGMAG (KG) 
!.....DISPERSION TENSOR                                                 
            ZERO = 0D0 
            ANGLE1 = ANGKV1 
            ANGLE2 = ANGKV2 
            ANGLE3 = 0D0 
            MSTRUC = 1 
            CALL TENSOR (DLG, ZERO, ZERO, ZERO, DTHG, ZERO, ZERO, ZERO, &
            DTVG, ANGLE1, ANGLE2, ANGLE3, DDXXG, DDXYG, DDXZG, DDYXG,   &
            DDYYG, DDYZG, DDZXG, DDZYG, DDZZG, MSTRUC)                  
            ANGLE1 = ElemData(ElemMap(L))%pangl1
            ANGLE2 = ElemData(ElemMap(L))%pangl2
            ANGLE3 = ElemData(ElemMap(L))%pangl3
            MSTRUC = 0 
            CALL TENSOR (DDXXG, DDXYG, DDXZG, DDYXG, DDYYG, DDYZG,      &
                         DDZXG, DDZYG, DDZZG, ANGLE1, ANGLE2, ANGLE3, DXXG, DXYG,    &
                         DXZG, DYXG, DYYG, DYZG, DZXG, DZYG, DZZG, MSTRUC)           
!                                                                       
!.....IN-PARALLEL CONDUCTIVITIES (DIFFUSIVITIES) FORMULA                
 6900       IF (KSP.EQ.NESP) THEN 
!..........FOR ENERGY TRANSPORT:                                        
               if(LVolAvgLambda) then
                 ESE = ESWG * SIGMAW (KSP) + (1D0 - PORG (KG) ) * ElemData(ElemMap(L))%lambdas  !Volumetric Average
               else
                 !lambda(bulk) = (la^(1-SW) * lw^Sw)^por * ls(1-por)
                 !lambda(air)  = 0.026 W/(m degC) - parameter in PARAMS
                 ESE = (SigmaAir**(1-SWG(KG)) * SIGMAW(KSP)**SWG(KG))**PORG(KG) * (ElemData(ElemMap(L))%lambdas**(1D0-PORG(KG)))       !Geometric Mean
               end if
            ELSE 
!..........FOR SOLUTE TRANSPORT:                                        
               ESE = ESRCG * SIGMAW (KSP) + (1D0 - PORG (KG) ) * RHOCWG * SIGMAS                                                 
            ENDIF 
!.....ADD DIFFUSION AND DISPERSION TERMS TO TOTAL DISPERSION TENSOR     
            BXXG (KG) = ESRCG * DXXG + ESE 
            BXYG (KG) = ESRCG * DXYG 
            BXZG (KG) = ESRCG * DXZG 
            BYXG (KG) = ESRCG * DYXG 
            BYYG (KG) = ESRCG * DYYG + ESE 
            BYZG (KG) = ESRCG * DYZG 
            BZXG (KG) = ESRCG * DZXG 
            BZYG (KG) = ESRCG * DZYG 
 7000       BZZG (KG) = ESRCG * DZZG + ESE 
!                                                                       
!.....INTEGRATE SOLUTE MASS BALANCE OR ENERGY BALANCE                   
!        USING SYMMETRIC WEIGHTING FUNCTIONS FOR DISPERSION TERM AND    
!        USING EITHER SYMMETRIC OR ASYMMETRIC WEIGHTING FUNCTIONS       
!        FOR ADVECTION TERM                                             
         DO 8000 I = 1, 8 
            DO 8000 J = 1, 8 
               BT = 0.D0 
               DT = 0.D0 
               DO 7500 KG = 1, 8 
                  BT = BT + ( (BXXG (KG) * DFDXG (J, KG) + BXYG (KG)            &
                  * DFDYG (J, KG) + BXZG (KG) * DFDZG (J, KG) ) * DFDXG (I, KG) &
                  + (BYXG (KG) * DFDXG (J, KG) + BYYG (KG) * DFDYG (J, KG)      &
                  + BYZG (KG) * DFDZG (J, KG) ) * DFDYG (I, KG) + (BZXG (KG)    &
                  * DFDXG (J, KG) + BZYG (KG) * DFDYG (J, KG) + BZZG (KG)       &
                  * DFDZG (J, KG) ) * DFDZG (I, KG) ) * DET (KG)                                  
 7500          DT = DT + (EXG (KG) * DFDXG (J, KG) + EYG (KG) * DFDYG (J, KG)   &
                    + EZG (KG) * DFDZG (J, KG) ) * W (I, KG) * DET (KG)                                                      
               BTRANE (I, J) = BT 
 8000    DTRANE (I, J) = DT 
 9000    CONTINUE 
!                                                                       
!                                                                       
!.....SEND RESULTS OF INTEGRATIONS FOR THIS ELEMENT                     
!        TO GLOBAL ASSEMBLY ROUTINE:                                    
!        GLOBAN -- SUTRA'S ORIGINAL BANDED FORMAT                       
!        GLOTRI -- SLAP TRIAD FORMAT                                    
!        GLO27  -- COMPRESSED FORMAT BASED ON 27-NODE "MOLECULE"        
         IF (KSOLVP.EQ.0) THEN 
            CALL GLOBAN (L, ML, VOLE, BFLOWE, DFLOWE, BTRANE, DTRANE)
         ELSE 
            IF ( .not.lColumnStorage ) THEN 
              CALL GLOTRI (L, ML, VOLE, BFLOWE, DFLOWE, BTRANE, DTRANE)              
            ELSE
              CALL GLOCOL (L,ML,VOLE,BFLOWE,DFLOWE,BTRANE,DTRANE)
            END IF
         ENDIF 
 9999 END DO 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
      RETURN 
      END SUBROUTINE ELEMN3                         
