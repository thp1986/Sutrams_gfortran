!     SUBROUTINE        B  C                   SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO IMPLEMENT SPECIFIED PRESSURE AND SPECIFIED TEMPERATURE OR     
! ***  CONCENTRATION CONDITIONS BY MODIFYING THE GLOBAL FLOW AND        
! ***  TRANSPORT MATRIX EQUATIONS.                                      
!                                                                       
      SUBROUTINE BC(ML)
      USE PARAMS 
      USE CONTRL 
      USE SOLVI 
      USE DIMS
      USE DIMX
      USE TIMES
      USE SutraStorage, ONLY: PMAT, PVEC, UMAT, UVEC, &
                              IPBC, PBC, IUBC, UBC,  &
                              QPLITR, MIOFF, JTRI, &
                              SpecifiedPBC, &
                              MultiSpeciesBC
      use SutraMSPrecision
      USE ColumnStorage
      implicit none
      !locals
      integer (I4B) :: &
        IMID0, &
        IMID, &
        JMID
      integer (I4B) :: &
        I, &
        IP, &
        IU, &
        IUP, &
        ML
      real (DP) :: &
        GPINL, &
        GPINR, &
        GUINL, &
        GUINR, &
        GUR, &
        GUL
             
      IF ( KSOLVP.EQ.0 ) THEN 
         IMID0 = 0 
         JMID = NBHALF 
      ELSE 
         IF ( .not.lColumnStorage ) THEN
             IF ( IABS(KTYPE).EQ.3 ) THEN 
                IMID0 = MIOFF(14) 
             ELSE 
                IMID0 = MIOFF(5) 
             ENDIF 
         END IF
         JMID = 1 
      ENDIF 
!                                                                       
      IF (NPBC.EQ.0) GOTO 1050 
!.....SPECIFIED P BOUNDARY CONDITIONS                                   
      DO 1000 IP = 1, NPBC 
         I = IABS( SpecifiedPBC(IP)%node ) 
         IMID = IMID0 + I 
         if ( KSOLVP.NE.0 .and. lColumnStorage ) IMID = JTRI(I)
!                                                                       
         IF (ML - 1) 100, 100, 200 
!.....MODIFY EQUATION FOR P BY ADDING FLUID SOURCE AT SPECIFIED         
!        PRESSURE NODE                                                  
  100    GPINL = - GNUP 
         GPINR = GNUP * SpecifiedPBC(IP)%P
         PMAT (IMID, JMID) = PMAT (IMID, JMID) - GPINL 
         PVEC (I) = PVEC (I) + GPINR 
!                                                                       
         IF (ML - 1) 200, 1000, 200 
!.....MODIFY EQUATION FOR U BY ADDING U SOURCE WHEN FLUID FLOWS IN      
!        AT SPECIFIED PRESSURE NODE                                     
  200    GUR = 0.0D0 
         GUL = 0.0D0 
         IF ( QPLITR(IP) ) 360, 360, 340 
  340    GUL = - CW  * QPLITR(IP) 
         GUR = - GUL * SpecifiedPBC(IP)%U(KSP) 
  360    IF (NOUMAT) 370, 370, 380 
  370    UMAT (IMID, JMID) = UMAT (IMID, JMID) - GUL 
  380    UVEC (I, KSP) = UVEC (I, KSP) + GUR 
! 1000 END DO 
 1000 CONTINUE 
!                                                                       
!                                                                       
 1050 IF (ML - 1) 1100, 3000, 1100 
 1100 IF ( NUBC(KSP).EQ.0 ) GOTO 3000 
!.....SPECIFIED U BOUNDARY CONDITIONS                                   
!        MODIFY EQUATION FOR U BY ADDING ENERGY/SOLUTE MASS SOURCE      
!        AT SPECIFIED U NODE                                            
      DO 2500 IU = 1, NUBC(KSP) 
!         IUP = IU + NPBC 
         I = IABS ( MultiSpeciesBC(KSP)%SpecifiedU(IU)%node )
         IMID = IMID0 + I 
         if ( KSOLVP.NE.0 .and. lColumnStorage ) IMID = JTRI(I)
         IF (NOUMAT) 1200, 1200, 2000 
 1200    GUINL = - GNUU (KSP) 
         UMAT(IMID, JMID) = UMAT(IMID, JMID) - GUINL 
 2000    GUINR = GNUU(KSP) * MultiSpeciesBC(KSP)%SpecifiedU(IU)%U
 2500    UVEC(I, KSP) = UVEC(I, KSP) + GUINR 
 
!                                                                       
 3000 CONTINUE 
!                                                                       
!                                                                       
      RETURN 
      END SUBROUTINE BC                             
