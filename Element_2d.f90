!     SUBROUTINE        E  L  E  M  N  2       SUTRA-MS VERSION 2004.1
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
! ***  THIS SUBROUTINE HANDLES 2-D CALCULATIONS ONLY.                   
!                                                                       
!                                                                       
      SUBROUTINE ELEMN2 (ML)
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
      USE MSErrorHandler
      USE SutraMSPrecision
      USE SutraZoneModule
      USE SutraStorage, ONLY: IN, X, Y, THICK=>Z, PITER, UITER, RCIT, RCITM1, &
                              VMAG, VANG=>VANG1, VOL, PMAT, PVEC, UMAT, UVEC, &
                              GXSI, GETA, PVEL, NBI27, MIOFF, ITRI, JTRI
      USE TECPLOT
      USE ColumnStorage

      implicit none
      integer (I4B) :: &
        ML
      real (DP) :: &
        BFLOWE (8, 8), DFLOWE (8), BTRANE (8, 8), DTRANE (8, 8), VOLE (8)                                                          
      real (DP) :: &
        F (4, 4), W (4, 4), DET (4), DFDXG (4, 4), DFDYG (4, 4),&
        DWDXG (4, 4), DWDYG (4, 4)                                        
      real (DP) :: &
        SWG (4), RHOG (4), VISCG (4), PORG (4), VXG (4),        &
        VYG (4), RELKG (4), RGXG (4), RGYG (4), VGMAG (4), THICKG (4)     
      real (DP) :: &
        RXXG (4), RXYG (4), RYXG (4), RYYG (4) 
      real (DP) :: &
        BXXG (4), BXYG (4), BYXG (4), BYYG (4), EXG (4), EYG (4)                                                           
      real (DP) :: &
        GXLOC (4), GYLOC (4) 
      REAL (DP) :: &
        ESE
      !locals
      integer (I4B) :: &
        INTIM,  ISTOP, &
        IVPRNT, KVPRNT
      integer (I4B) :: &
        I,    I1,   II,   IL, &
        IXL, IYL, &
        J, &
        KG, &
        L
      real (DP) :: &
        GLOC, &
        XLOC, YLOC
      real (DP) :: &
        CJ11, CJ12, CJ21, CJ22
      real (DP) :: &
        XIX,  YIY
      real (DP) :: &
        AXSUM, AYSUM
      real (DP) :: &
        SWTEST, &
        ROMG, &
        DF, &
        VO, &
        BF, &
        ESWG, &
        RHOCWG, &
        ESRCG
      real (DP) :: &
        DXXG,  DXYG, DYXG, DYYG, &
        ALEFF, ATEFF, &
        DLG,   DTG, &
        BT,    DT
      real (DP) :: &
        VANGG, VKANGG, &
        VXVG,  VYVG, &
        VXVG2, VYVG2, &
        DCO,   DSI

      !data
      DATA GLOC / 0.577350269189626D0 / 
      DATA INTIM / 0 /, ISTOP / 0 /, &
            GXLOC / - 1.D0, 1.D0, 1.D0, - 1.D0 /, &
            GYLOC / - 1.D0, - 1.D0, 1.D0, 1.D0 /                    
      SAVE GLOC, INTIM, ISTOP, GXLOC, GYLOC 
!
!
      MSErrorValue%cDataSet='EL2'
!
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
!        GXSI AND GETA, FOR CONSISTENT VELOCITIES,                      
!        AND CHECK ELEMENT SHAPES                                       
      IF (INTIM) 100, 100, 2000 
  100 INTIM = 1 
!.....LOOP THROUGH ALL ELEMENTS TO OBTAIN THE JACOBIAN                  
!        AT EACH OF THE FOUR NODES IN EACH ELEMENT                      
!                                                                       
!.......SCREEN OUTPUT                                               
      WRITE ( * , * ) 'CALCULATING 2-D JACOBIAN ON FIRST PASS' 
!.......SCREEN OUTPUT                                               
!                                                                       
      DO 1000 L = 1, NE 
         DO 500 IL = 1, 4 
            XLOC = GXLOC (IL) 
            YLOC = GYLOC (IL) 
            CALL BASIS2 (0000, L, XLOC, YLOC, F (1, IL),      &
                         W (1, IL), DET (IL), DFDXG (1, IL), DFDYG (1, IL), DWDXG (1,IL), &
                         DWDYG (1, IL), THICKG (IL), &
                         VXG (IL), VYG (IL), SWG (IL), RHOG (IL), VISCG (IL),   &
                         PORG (IL), VGMAG (IL), RELKG (IL),  &
                         CJ11, CJ12, CJ21, CJ22, &
                         RGXG (IL), RGYG (IL))                                 
            GXSI (L, IL) = CJ11 * GRAVX + CJ12 * GRAVY 
            GETA (L, IL) = CJ21 * GRAVX + CJ22 * GRAVY 
!.....CHECK FOR NEGATIVE- OR ZERO-AREA ERRORS IN ELEMENT SHAPES         
            IF (DET (IL) ) 200, 200, 500 
  200       ISTOP = ISTOP + 1 

            II = (L - 1) * N48 + IL 
            I = IN (II)

            WRITE (fLST, 400) IN ( (L - 1) * 4 + IL), L, DET (IL) 
  400 FORMAT  (11X,'THE DETERMINANT OF THE JACOBIAN AT NODE ',I9,       &
     &     ' IN ELEMENT ',I9,' IS NEGATIVE OR ZERO, ',1PE15.7)          
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
      call ErrorIO('ELEMN2: Some elements have incorrect geometry')
!                                                                       
!.....LOOP THROUGH ALL ELEMENTS TO CARRY OUT SPATITRIL INTEGRATION        
!        OF FLUX TERMS IN P AND/OR U EQUATIONS                          
 2000 IF (IUNSAT.NE.0) IUNSAT = 2 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SCREEN OUTPUT                                                                    
      WRITE ( * , * ) 'PREPARING TO CALCULATE WEIGHTING FUNCTIONS' 
      I1 = - 1 
!.....SCREEN OUTPUT
!                                                                       
      DO 9999 L = 1, NE 
         XIX = - 1.D0 
         YIY = - 1.D0 
         KG = 0 
!.....OBTAIN BASIS FUNCTION AND RELATED INFORMATION AT EACH OF          
!        FOUR GAUSS POINTS IN THE ELEMENT                               
!                                                                       
!.....SCREEN OUTPUT
         I1 = I1 + 1 
         IF (I1.EQ.0) THEN 
           WRITE ( * ,  * ) 'CALCULATING 2-D BASIS AND WEIGHTING FUNCTIONS' 
         ENDIF 
!.....SCREEN OUTPUT
!                                                                       
         DO 2200 IYL = 1, 2 
            DO 2100 IXL = 1, 2 
               KG = KG + 1 
               XLOC = XIX * GLOC 
               YLOC = YIY * GLOC 
               CALL BASIS2 (0001, L, XLOC, YLOC, F (1, KG),   &
                            W (1, KG), DET (KG), DFDXG (1, KG), DFDYG (1, KG),       &
                            DWDXG (1, KG), DWDYG (1, KG), &
                            THICKG (KG), VXG (KG), VYG (KG), SWG (KG),        &
                            RHOG (KG), VISCG (KG), PORG (KG), VGMAG (KG), RELKG (KG),&
                            CJ11, CJ12, CJ21, CJ22,  &
                            RGXG (KG), RGYG (KG))    
 2100       XIX = - XIX 
 2200    YIY = - YIY 
!.....VERSION 1.1
!.....SAVE NODAL VELOCITIES IF REQUIRED
      IF (LTpPrintVelocityAtNodes) THEN
        DO KG = 1, N48
         II = (L - 1) * N48 + KG
         I = IN (II) 
         CALL MakeVelocityVector(I, VXG(KG), VYG(KG))
        END DO
      END IF
!                                                                       
!.....CALCULATE VELOCITY AT ELEMENT CENTROID WHEN REQUIRED              
         IF (KVPRNT - 2) 3000, 2300, 3000 
 2300    AXSUM = 0.0D0 
         AYSUM = 0.0D0 
         DO 2400 KG = 1, 4 
            AXSUM = AXSUM + VXG (KG) 
 2400    AYSUM = AYSUM + VYG (KG) 
         VMAG (L) = DSQRT (AXSUM * AXSUM + AYSUM * AYSUM) 
         IF (VMAG (L) .NE.0D0) THEN 
            VMAG (L) = VMAG (L) * 2.5D-1 
            VANG (L) = DATAN2 (AYSUM, AXSUM) * 5.729577951308232D+1 
         ELSE 
            VANG (L) = 0D0 
         ENDIF 
!                                                                       
!.....INCLUDE MESH THICKNESS IN NUMERICAL INTEGRATION                   
 3000    DO 3300 KG = 1, 4 
 3300    DET (KG) = THICKG (KG) * DET (KG) 
!                                                                       
!.....CALCULATE PARAMETERS FOR FLUID MASS BALANCE AT GAUSS POINTS       
         IF (ML - 1) 3400, 3400, 6100 
 3400    SWTEST = 0.D0 
!                                                                       
!.....SCREEN OUTPUT
!         IF (I1.EQ.0) THEN 
!            WRITE ( * , * ) 'CALCULATING 2-D PARAMETERS FOR FLUID MASS' 
!            WRITE ( * , * ) 'BALANCE AT GAUSS POINTS' 
!         ENDIF 
!.....SCREEN OUTPUT
!                                                                       
         DO 4000 KG = 1, 4 
            SWTEST = SWTEST + SWG (KG) 
            ROMG = RHOG (KG) * RELKG (KG) / VISCG (KG) 
            RXXG (KG) = ElemData(ElemMap(L))%permxx * ROMG 
            RXYG (KG) = ElemData(ElemMap(L))%permxy * ROMG 
            RYXG (KG) = ElemData(ElemMap(L))%permyx * ROMG 
            RYYG (KG) = ElemData(ElemMap(L))%permyy * ROMG 
 4000    END DO 
!                                                                       
!.....INTEGRATE FLUID MASS BALANCE IN AN UNSATURATED ELEMENT            
!        USING ASYMMETRIC WEIGHTING FUNCTIONS                           
         IF (UP.LE.1.0D-6) GOTO 5200 
         IF (SWTEST - 3.999D0) 4200, 5200, 5200 
 4200    DO 5000 I = 1, 4 
            DF = 0.D0 
            VO = 0.D0 
            DO 4400 KG = 1, 4 
               VO = VO + F (I, KG) * DET (KG) 
 4400       DF = DF + ( (RXXG (KG) * RGXG (KG) + RXYG (KG) * RGYG (KG) )&
            * DWDXG (I, KG) + (RYXG (KG) * RGXG (KG) + RYYG (KG)        &
            * RGYG (KG) ) * DWDYG (I, KG) ) * DET (KG)                  
            DO 4800 J = 1, 4 
               BF = 0.D0 
               DO 4600 KG = 1, 4 
 4600          BF = BF + ( (RXXG (KG) * DFDXG (J, KG) + RXYG (KG)       &
               * DFDYG (J, KG) ) * DWDXG (I, KG) + (RYXG (KG) * DFDXG ( &
               J, KG) + RYYG (KG) * DFDYG (J, KG) ) * DWDYG (I, KG) )   &
               * DET (KG)                                               
 4800       BFLOWE (I, J) = BF 
            VOLE (I) = VO 
 5000    DFLOWE (I) = DF 
         GOTO 6200 
!                                                                       
!.....INTEGRATE FLUID MASS BALANCE IN A SATURATED OR UNSATURATED        
!        ELEMENT USING SYMMETRIC WEIGHTING FUNCTIONS                    
 5200    DO 6000 I = 1, 4 
            DF = 0.D0 
            VO = 0.D0 
            DO 5400 KG = 1, 4 
               VO = VO + F (I, KG) * DET (KG) 
 5400       DF = DF + ( (RXXG (KG) * RGXG (KG) + RXYG (KG) * RGYG (KG) )&
            * DFDXG (I, KG) + (RYXG (KG) * RGXG (KG) + RYYG (KG)        &
            * RGYG (KG) ) * DFDYG (I, KG) ) * DET (KG)                  
            DO 5800 J = 1, 4 
               BF = 0.D0 
               DO 5600 KG = 1, 4 
 5600          BF = BF + ( (RXXG (KG) * DFDXG (J, KG) + RXYG (KG)       &
               * DFDYG (J, KG) ) * DFDXG (I, KG) + (RYXG (KG) * DFDXG ( &
               J, KG) + RYYG (KG) * DFDYG (J, KG) ) * DFDYG (I, KG) )   &
               * DET (KG)                                               
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
!            WRITE ( * , * ) 'CALCULATING 2-D PARAMETERS FOR ENERGY OR' 
!            WRITE ( * , * ) 'SOLUTE MASS BALANCE AT GAUSS POINTS' 
!            WRITE ( * , * ) 'SPECIES [', KSP, ']' 
!         ENDIF 
!.....SCREEN OUTPUT
!                                                                       
         DO 7000 KG = 1, 4 
            ESWG = PORG (KG) * SWG (KG) 
            RHOCWG = RHOG (KG) * CW 
            ESRCG = ESWG * RHOCWG 
            IF (VGMAG (KG) ) 6300, 6300, 6600 
 6300       EXG (KG) = 0.0D0 
            EYG (KG) = 0.0D0 
            DXXG = 0.0D0 
            DXYG = 0.0D0 
            DYXG = 0.0D0 
            DYYG = 0.0D0 
            GOTO 6900 
 6600       EXG (KG) = ESRCG * VXG (KG) 
            EYG (KG) = ESRCG * VYG (KG) 
!                                                                       
!.....DISPERSIVITY MODEL FOR ANISOTROPIC MEDITRI                          
!        WITH PRINCIPAL DISPERSIVITIES:  ALMAX,ALMIN, AND ATAVG         
            VANGG = 1.570796327D0 
            IF (VXG (KG) * VXG (KG) .GT.0.D0) VANGG = DATAN (VYG (KG) / VXG (KG) )                                                
            VKANGG = VANGG - ElemData(ElemMap(L))%pangl1
            DCO = DCOS (VKANGG) 
            DSI = DSIN (VKANGG) 
!.....EFFECTIVE LONGITUDINAL DISPERSIVITY IN FLOW DIRECTION, ALEFF      
            ALEFF = 0.0D0 
            IF (ElemData(ElemMap(L))%almax + ElemData(ElemMap(L))%almin) 6800, 6800, 6700 
 6700       ALEFF = ElemData(ElemMap(L))%almax * ElemData(ElemMap(L))%almin / (ElemData(ElemMap(L))%almin * DCO * DCO +    &
            ElemData(ElemMap(L))%almax * DSI * DSI)                                      
!.......SCALE LONGITUDINAL DISPERSIVITY FOR SPECIES KSP                 
            ALEFF = ALEFF * ATSPMULT (KSP) 
 6800       DLG = ALEFF * VGMAG (KG) 
!.....EFFECTIVE TRANSVERSE DISPERSIVITY IN FLOW DIRECTION, ATEFF        
            ATEFF = 0.0D0 
            IF (ElemData(ElemMap(L))%at1max + ElemData(ElemMap(L))%at1min ) 6860, 6860, 6840 
 6840       ATEFF = ElemData(ElemMap(L))%at1max * ElemData(ElemMap(L))%at1min / (ElemData(ElemMap(L))%at1min * DCO * DCO +    &
            ElemData(ElemMap(L))%at1max * DSI * DSI)                                      
!.......SCALE TRANSVERSE DISPERSIVITY FOR SPECIES KSP                   
            ATEFF = ATEFF * ATSPMULT (KSP) 
 6860       DTG = ATEFF * VGMAG (KG) 
!                                                                       
            VXVG = VXG (KG) / VGMAG (KG) 
            VYVG = VYG (KG) / VGMAG (KG) 
            VXVG2 = VXVG * VXVG 
            VYVG2 = VYVG * VYVG 
!.....DISPERSION TENSOR                                                 
            DXXG = DLG * VXVG2 + DTG * VYVG2 
            DYYG = DTG * VXVG2 + DLG * VYVG2 
            DXYG = (DLG - DTG) * VXVG * VYVG 
            DYXG = DXYG 
!                                                                       
!.....IN-PARALLEL CONDUCTIVITIES (DIFFUSIVITIES) FORMULA                
 6900       IF (KSP.EQ.NESP) THEN 
!..........FOR ENERGY TRANSPORT:                                        
               if(LVolAvgLambda) then
                 ESE = ESWG * SIGMAW(KSP) + (1D0 - PORG (KG) ) * ElemData(ElemMap(L))%lambdas  !Volumetric Average
               else
                 !lambda(bulk) = (la^(1-SW) * lw^Sw)^por * ls(1-por)
                 !lambda(air)  = 0.026 W/(m degC) - parameter in PARAMS
                 ESE = (SigmaAir**(1-SWG(KG)) * SIGMAW(KSP)**SWG(KG))**PORG(KG) * (ElemData(ElemMap(L))%lambdas**(1D0-PORG(KG)))       !Geometric Mean
               end if
            ELSE 
!..........FOR SOLUTE TRANSPORT:                                        
               ESE = ESRCG * SIGMAW(KSP) 
            ENDIF 
!.....ADD DIFFUSION AND DISPERSION TERMS TO TOTAL DISPERSION TENSOR     
            BXXG (KG) = ESRCG * DXXG + ESE 
            BXYG (KG) = ESRCG * DXYG 
            BYXG (KG) = ESRCG * DYXG 
 7000    BYYG (KG) = ESRCG * DYYG + ESE 
!                                                                       
!.....INTEGRATE SOLUTE MASS BALANCE OR ENERGY BALANCE                   
!        USING SYMMETRIC WEIGHTING FUNCTIONS FOR DISPERSION TERM AND    
!        USING EITHER SYMMETRIC OR ASYMMETRIC WEIGHTING FUNCTIONS       
!        FOR ADVECTION TERM                                             
         DO 8000 I = 1, 4 
            DO 8000 J = 1, 4 
               BT = 0.D0 
               DT = 0.D0 
               DO 7500 KG = 1, 4 
                  BT = BT + ( (BXXG (KG) * DFDXG (J, KG) + BXYG (KG)    &
                  * DFDYG (J, KG) ) * DFDXG (I, KG) + (BYXG (KG)        &
                  * DFDXG (J, KG) + BYYG (KG) * DFDYG (J, KG) ) * DFDYG &
                  (I, KG) ) * DET (KG)                                  
 7500          DT = DT + (EXG (KG) * DFDXG (J, KG) + EYG (KG) * DFDYG ( &
               J, KG) ) * W (I, KG) * DET (KG)                          
               BTRANE (I, J) = BT 
 8000    DTRANE (I, J) = DT 
 9000    CONTINUE 
!                                                                       
!.....SCREEN OUTPUT
!         IF (I1.EQ.0) THEN 
!            WRITE ( * , * ) 'PLACING INTEGRATION RESULTS FOR ELEMENT' 
!            WRITE ( * , * ) 'IN GLOBAL MATRICES' 
!         ENDIF 
!.....SCREEN OUTPUT
!                                                                       
!                                                                       
!                                                                       
!.....PLACE RESULTS OF INTEGRATIONS FOR THIS ELEMENT IN                 
!        GLOBAL MATRICES (FORMERLY SUBSROUTINE GLOBAN)                  
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
!                                                                       
 9999 END DO 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!.....SCREEN OUTPUT
      WRITE ( * , * ) 'FINISHED CALCULATING WEIGHTING FUNCTIONS' 
!.....SCREEN OUTPUT
      RETURN 
      END SUBROUTINE ELEMN2                         
