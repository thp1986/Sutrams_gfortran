!     SUBROUTINE        N  O  D  A  L          SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  (1) TO CARRY OUT ALL CELLWISE CALCULATIONS AND TO ADD CELLWISE   
! ***      TERMS TO THE GLOBAL MATRIX AND GLOBAL VECTOR FOR BOTH FLOW   
! ***      AND TRANSPORT EQUATIONS.                                     
! ***  (2) TO ADD FLUID SOURCE AND SOLUTE MASS OR ENERGY SOURCE TERMS   
! ***      TO THE MATRIX EQUATIONS.                                     
!                                                                       
      SUBROUTINE NODAL (ML)
      USE PARAMS 
      USE CONTRL 
      USE SOLVI 
      USE DIMS
      USE DIMX
      USE TIMES
      USE SutraZoneModule
      USE SutraStorage, ONLY: VOL, PMAT, PVEC, UMAT, UVEC, PITER, UITER,  &
                              PM1, UM1, UM2, QIN, UIN, QUIN, &
                              CS1, CS2, CS3, SL, SR, SW,    &
                              DSWDP, RHO, MIOFF, JTRI
      use SutraMSPrecision
      USE ColumnStorage
      implicit none
      integer (I4B) :: &
        ML
      !locals
      integer (I4B) :: &
        I, IMID, IMID0, &
        JMID, &
        K
      integer (I4B) :: &
        imap
      real (DP) :: &
        RELK, SWRHON, &
        AFLN, CFLNT, CFLN, &
        DUDT, EPRS, &
        ATRN, GTRN, GSV, GSLTRN, GSRTRN, ETRN, QUR, QUL
!                                                                       
      IF (IUNSAT.NE.0) IUNSAT = 1 
!                                                                       
      IMID0 = 0 
      IF (KSOLVP.EQ.0) THEN 
         IMID0 = 0 
         JMID = NBHALF 
      ELSE 
         IF( .not.lColumnStorage ) THEN
           IF ( IABS(KTYPE).EQ.3 ) THEN 
              IMID0 = MIOFF(14) 
           ELSE 
              IMID0 = MIOFF(5) 
           ENDIF 
         END IF
         JMID = 1 
      ENDIF 
!                                                                       
!.....DO NOT UPDATE NODAL PARAMETERS ON A TIME STEP WHEN ONLY U IS      
!     SOLVED FOR BY BACK SUBSTITUTION (IE: WHEN NOUMAT=1)               
      IF (NOUMAT) 50, 50, 200 
!.....SET UNSATURATED FLOW PARAMETERS AT NODES, SW(I) AND DSWDP(I)      
   50 DO 120 I = 1, NN 
         IF (IUNSAT - 1) 120, 100, 120 
  100    IF (PITER (I) ) 110, 115, 115 
  110    CALL UNSAT (SW (I), DSWDP (I), RELK, PITER (I), NodeMap(I) ) 
         GOTO 120 
  115    SW(I)    = 1.0D0 
         DSWDP(I) = 0.0D0 
  120 END DO 
!.....SET FLUID DENSITY AT NODES, RHO(I)                                
!     RHO = F (UITER(I))                                                
      DO 150 I = 1, NN 
!.......CALCULATE FLUID DENSITY BASED ON PARAMETERS FOR EACH SPECIES    
         RHO (I) = RHOW0 
         DO 140 K = 1, NSPE 
  140    RHO(I) = RHO(I) + DRWDU(K) * ( UITER(I, K) - URHOW0(K) ) 
  150 END DO 
  200 CONTINUE 
!                                                                       
      DO 1000 I = 1, NN 
         imap = NodeMap(i)
         if( lColumnStorage .AND. KSOLVP.NE.0 ) then
           IMID = JTRI(I)
         else
           IMID = IMID0 + I 
         end if
!         IF(KSOLVP.NE.0) then
!           IMID = JTRI(I)
!         ELSE
!           IMID = IMID0 + I 
!         END IF
!                                                                       
         SWRHON = SW(I) * RHO(I) 
!                                                                       
         IF (ML - 1) 220, 220, 230 
!                                                                       
!.....CALCULATE CELLWISE TERMS FOR P EQUATION                           
!.....FOR STEADY-STATE FLOW, ISSFLO=2; FOR TRANSIENT FLOW, ISSFLO=0     
  220    AFLN = (1 - ISSFLO / 2) * (SWRHON * NodeData(imap)%sop + NodeData(imap)%por * &
                RHO(I) * DSWDP(I) ) * VOL(I) / DELTP                                
         CFLNT = 0.0D0 
         DO 225 K = 1, NSPE 
            CFLN = NodeData(imap)%por * SW(I) * DRWDU(K) * VOL(I) 
            DUDT = (1 - ISSFLO / 2) * ( UM1(I, K) - UM2(I, K) ) / DLTUM1                                                    
            CFLNT = CFLNT + CFLN * DUDT 
  225    END DO 
!.....ADD CELLWISE TERMS AND FLUID SOURCES OR FLUXES TO P EQUATION      
         PMAT (IMID, JMID) = PMAT (IMID, JMID) + AFLN 
         PVEC (I) = PVEC (I) - CFLNT + AFLN * PM1(I) + QIN(I) 
!                                                                       
         IF (ML - 1) 230, 1000, 230 
!                                                                       
!.....CALCULATE CELLWISE TERMS FOR U-EQUATION                           
  230    EPRS = (1.D0 - NodeData(imap)%por ) * NodeData(imap)%rhos 
         ATRN = (1 - ISSTRA) * (NodeData(imap)%por * SWRHON * CW + EPRS * CS1(I, KSP) ) * &
                VOL(I) / DELTU                                       
!         GTRN = NodeData(imap)%por * SWRHON * PRODF1 (KSP) * VOL (I) 
!         GSV = EPRS * PRODS1 (KSP) * VOL (I) 
         GTRN = NodeData(imap)%por * SWRHON * ProdSorp(imap)%prodf1(KSP) * VOL (I) 
         GSV = EPRS * ProdSorp(imap)%prods1(KSP) * VOL (I) 
         GSLTRN = GSV * SL (I) 
         GSRTRN = GSV * SR (I) 
!         ETRN = (NodeData(imap)%por * SWRHON * PRODF0 (KSP) + EPRS * PRODS0 (KSP) ) * VOL (I)                                                      
         ETRN = (NodeData(imap)%por * SWRHON * ProdSorp(imap)%prodf0(KSP) + EPRS * ProdSorp(imap)%prods0(KSP) ) * VOL (I)                                                      
!.....CALCULATE SOURCES OF SOLUTE OR ENERGY CONTAINED IN                
!     SOURCES OF FLUID (ZERO CONTRIBUTION FOR OUTFLOWING FLUID)         
         QUR = 0.0D0 
         QUL = 0.0D0 
         IF ( QIN(I) ) 360, 360, 340 
  340    QUL = - CW  * QIN(I) 
         QUR = - QUL * UIN(I, KSP) 
!.....ADD CELLWISE TERMS, SOURCES OF SOLUTE OR ENERGY IN FLUID INFLOWS, 
!        AND PURE SOURCES OR FLUXES OF SOLUTE OR ENERGY TO U-EQUATION   
  360    IF (NOUMAT) 370, 370, 380 
  370    UMAT (IMID, JMID) = UMAT (IMID, JMID) + ATRN - GTRN - GSLTRN - QUL                                                            
  380    UVEC (I, KSP) = UVEC (I, KSP) + ATRN * UM1 (I, KSP) + ETRN + &
                         GSRTRN + QUR + QUIN (I, KSP)                                   
!                                                                       
 1000 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE NODAL                          
