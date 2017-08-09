!     SUBROUTINE        B  C  T  I  M  E       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  USER-PROGRAMMED SUBROUTINE WHICH ALLOWS THE USER TO SPECIFY:     
! ***   (1) TIME-DEPENDENT SPECIFIED PRESSURES AND TIME-DEPENDENT       
! ***       CONCENTRATIONS OR TEMPERATURES OF INFLOWS AT THESE POINTS   
! ***   (2) TIME-DEPENDENT SPECIFIED CONCENTRATIONS OR TEMPERATURES     
! ***   (3) TIME-DEPENDENT FLUID SOURCES AND CONCENTRATIONS             
! ***       OR TEMPERATURES OF INFLOWS AT THESE POINTS                  
! ***   (4) TIME-DEPENDENT ENERGY OR SOLUTE MASS SOURCES                
!                                                                       
      SUBROUTINE BCTIME (IPBCT, IUBCT, IQSOPT, IQSOUT)
      USE PARAMS 
      USE FUNITS 
      USE DIMS
      USE TIMES
      USE GRAVEC
      USE SutraStorage, ONLY : IPBC, PBC, IUBC, UBC, &
                               QIN, UIN, QUIN, IQSOP, IQSOU, &
                               X, Y, Z, &
                               SpecifiedPBC, &
                               MultiSpeciesBC
      USE SutraMSPrecision
      USE M_TIDE


      IMPLICIT NONE

      integer (I4B) :: &
        IPBCT, IUBCT, IQSOPT, IQSOUT

!     LOCAL VARIABLES
      INTEGER (I4B) :: &
        I, IP, IU, IUP, IQP, IQU, &
        K, &
        NSOPI, NSOUI
      REAL (DP) :: &
	TDLEVEL

!.....DEFINITION OF REQUIRED VARIABLES                                  
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     NN      = EXACT NUMBER OF NODES IN MESH                           
!     NPBC    = EXACT NUMBER OF SPECIFIED PRESSURE NODES                
!     NUBC(K) = EXACT NUMBER OF SPECIFIED CONCENTRATION                 
!               OR TEMPERATURE NODES FOR EACH SPECIES                   
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     IT = NUMBER OF CURRENT TIME STEP                                  
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     TSEC = TIME AT END OF CURRENT TIME STEP IN SECONDS                
!     TMIN = TIME AT END OF CURRENT TIME STEP IN MINUTES                
!     THOUR = TIME AT END OF CURRENT TIME STEP IN HOURS                 
!     TDAY = TIME AT END OF CURRENT TIME STEP IN DAYS                   
!     TWEEK = TIME AT END OF CURRENT TIME STEP IN WEEKS                 
!     TMONTH = TIME AT END OF CURRENT TIME STEP IN MONTHS               
!     TYEAR = TIME AT END OF CURRENT TIME STEP IN YEARS                 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     PBC(IP)   = SPECIFIED PRESSURE VALUE AT IP(TH) SPECIFIED          
!                 PRESSURE NODE                                         
!     UBC(IP,K) = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE OF ANY   
!                 INFLOW OCCURRING AT IP(TH) SPECIFIED PRESSURE NODE    
!                 FOR EACH SPECIES                                      
!     IPBC(IP)  = ACTUAL NODE NUMBER OF IP(TH) SPECIFIED PRESSURE NODE  
!                 {WHEN NODE NUMBER I=IPBC(IP) IS NEGATIVE (I<0),       
!                 VALUES MUST BE SPECIFIED FOR PBC AND UBC.}            
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     UBC(IUP,K)  = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE AT     
!                   IU(TH) SPECIFIED CONCENTRATION OR TEMPERATURE NODE  
!                   (WHERE IUP=IU+NPBC)                                 
!                   FOR EACH SPECIES                                    
!     IUBC(IUP,K) = ACTUAL NODE NUMBER OF IU(TH) SPECIFIED              
!                   CONCENTRATION OR TEMPERATURE NODE                   
!                   (WHERE IUP=IU+NPBC)                                 
!                   {WHEN NODE NUMBER I=IUBC(IU) IS NEGATIVE (I<0),     
!                   A VALUE MUST BE SPECIFIED FOR UBC.}                 
!                   FOR EACH SPECIES                                    
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     IQSOP(IQP) = NODE NUMBER OF IQP(TH) FLUID SOURCE NODE.            
!                  {WHEN NODE NUMBER I=IQSOP(IQP) IS NEGATIVE (I<0),    
!                  VALUES MUST BE SPECIFIED FOR QIN AND UIN.}           
!     QIN(-I)    = SPECIFIED FLUID SOURCE VALUE AT NODE (-I)            
!     UIN(-I,K)  = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE OF ANY  
!                  INFLOW OCCURRING AT FLUID SOURCE NODE (-I)           
!                  FOR EACH SPECIES                                     
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     IQSOU(IQU,K) = NODE NUMBER OF IQU(TH) ENERGY OR                   
!                    SOLUTE MASS SOURCE NODE                            
!                    {WHEN NODE NUMBER I=IQSOU(IQU) IS NEGATIVE (I<0),  
!                    A VALUE MUST BE SPECIFIED FOR QUIN.}               
!                    FOR EACH SPECIES                                   
!     QUIN(-I,K)  = SPECIFIED ENERGY OR SOLUTE MASS SOURCE VALUE        
!                   AT NODE (-I)                                        
!                   FOR EACH SPECIES                                    
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!                                                                       
!.....ADDITIONAL USEFUL VARIABLES                                       
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     "FUNITS" ARE UNIT NUMBERS FOR INPUT AND OUTPUT FILES              
!         AS ASSIGNED IN THE INPUT FILE, "SUTRA.FIL"                    
!                                                                       
!     X(I), Y(I), AND Z(I) ARE THE X-, Y-, AND Z-COORDINATES OF NODE I  
!     (FOR 2-D PROBLEMS, Z(I) IS THE THICKNESS AT NODE I)               
!                                                                       
!     GRAVX, GRAVY AND GRAVZ ARE THE X-, Y-, AND Z-COMPONENTS OF THE    
!     GRAVITY VECTOR                                                    
!     (FOR 2-D PROBLEMS, GRAVZ = 0)                                     
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!                                                                       
!                                                                       
!.....NSOPI IS ACTUAL NUMBER OF FLUID SOURCE NODES                      
      NSOPI = NSOP - 1 
!.....NSOUI IS ACTUAL NUMBER OF ENERGY OR SOLUTE MASS SOURCE NODES      
      NSOUI = MNSOU - 1 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      IF (IPBCT) 50, 240, 240 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (1):  SET TIME-DEPENDENT SPECIFIED PRESSURES OR           
!     CONCENTRATIONS (TEMPERATURES) OF INFLOWS AT SPECIFIED             
!     PRESSURE NODES                                                    
!                                                                       
   50 CONTINUE 
      DO 200 IP = 1, NPBC 
         I = SpecifiedPBC(IP)%node
         IF (I) 100, 200, 200 
  100    CONTINUE 
!     NOTE : A FLOW AND TRANSPORT SOLUTION MUST OCCUR FOR ANY           
!            TIME STEP IN WHICH PBC( ) CHANGES.                         
!     SpecifiedPBC(IP)%P =  ((          ))                                         
!     DO 150 K=1,NSPE                                                   
! 150 SpecifiedPBC(IP)%U(K) =  ((          ))        
!.....TIDAL ACTION
	TDLEVEL=MSL+TAMP*SIN(2*3.1415926*TSEC/TPESEC)     
	IF ((Y(IABS(I))).LE.TDLEVEL)THEN
        	SpecifiedPBC(IP)%P = PBCRHO*9.81*(TDLEVEL-Y(IABS(I)))
	ELSE 
        	SpecifiedPBC(IP)%P = 0
	ENDIF
!.....TEMPERATURE
      DO 150 K=1, NSPE
      IF(K.EQ.1) THEN
         SpecifiedPBC(IP)%U(K) = PBCTEMP
!.....CONCENTRATION
      ELSEIF(K.EQ.2) THEN
            SpecifiedPBC(IP)%U(K) = PBCSAL
      ENDIF
  150 END DO   
  200 END DO 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
  240 IF (IUBCT) 250, 440, 440 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (2):  SET TIME-DEPENDENT SPECIFIED                        
!     CONCENTRATIONS (TEMPERATURES)                                     
!                                                                       
  250 CONTINUE 
      DO 400 K = 1, NSPE 
         DO 350 IU = 1, NUBC (K) 
            IUP = NPBC + IU 
            I = MultiSpeciesBC(K)%SpecifiedU(IU)%node
            IF (I) 300, 400, 400 
  300       CONTINUE 
!       NOTE : A TRANSPORT SOLUTION MUST OCCUR FOR ANY TIME STEP        
!              IN WHICH UBC( ) CHANGES.  IN ADDITION, IF FLUID          
!              PROPERTIES ARE SENSITIVE TO 'U' THEN A FLOW SOLUTION     
!              MUST OCCUR AS WELL                                       
!           MultiSpeciesBC(K)%SpecifiedU(IU)%U =   ((          ))                                   
  350    END DO 
  400 END DO 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
  440 IF (IQSOPT) 450, 640, 640 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (3):  SET TIME-DEPENDENT FLUID SOURCES/SINKS,             
!      OR CONCENTRATIONS (TEMPERATURES) OF SOURCE FLUID                 
!                                                                       
  450 CONTINUE 
      DO 600 IQP = 1, NSOPI 
         I = IQSOP (IQP) 
         IF (I) 500, 600, 600 
  500    CONTINUE 
!     NOTE : A FLOW AND TRANSPORT SOLUTION MUST OCCUR FOR ANY           
!            TIME STEP IN WHICH QIN( ) CHANGES.                         
!     QIN(-I) =   ((           ))                                       
!     NOTE : A TRANSPORT SOLUTION MUST OCCUR FOR ANY                    
!            TIME STEP IN WHICH UIN( ) CHANGES.                         
!     DO 550 K=1,NSPE                                                   
!       UIN(-I,K) =   ((           ))                                   
! 550 CONTINUE                                                          
  600 END DO 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
  640 IF (IQSOUT) 650, 840, 840 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (4):  SET TIME-DEPENDENT SOURCES/SINKS                    
!     OF SOLUTE MASS OR ENERGY                                          
!                                                                       
  650 CONTINUE 
      DO 800 K = 1, NSPE 
         DO 750 IQU = 1, NSOU (K) 
            I = IQSOU (IQU, K) 
            IF (I) 700, 800, 800 
  700       CONTINUE 
!       NOTE : A TRANSPORT SOLUTION MUST OCCUR FOR ANY                  
!              TIME STEP IN WHICH QUIN( ) CHANGES.                      
!       QUIN(-I,K) =   ((           ))                                  
  750    END DO 
  800 END DO 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
  840 CONTINUE 
!                                                                       
      RETURN 
      END SUBROUTINE BCTIME 

!                                                                       
!     SUBROUTINE        U  N  S  A  T              SUTRA VERSION 2D3D.1 
!                                                                       
! *** PURPOSE :                                                         
! ***  USER-PROGRAMMED SUBROUTINE GIVING:                               
! ***  (1)  SATURATION AS A FUNCTION OF PRESSURE ( SW(PRES) )           
! ***  (2)  DERIVATIVE OF SATURATION WITH RESPECT TO PRESSURE           
! ***       AS A FUNCTION OF EITHER PRESSURE OR SATURATION              
! ***       ( DSWDP(PRES), OR DSWDP(SW) )                               
! ***  (3)  RELATIVE PERMEABILITY AS A FUNCTION OF EITHER               
! ***       PRESSURE OR SATURATION ( REL(PRES) OR RELK(SW) )            
! ***                                                                   
! ***  CODE BETWEEN DASHED LINES MUST BE REPLACED TO GIVE THE           
! ***  PARTICULAR UNSATURATED RELATIONSHIPS DESIRED.                    
! ***                                                                   
! ***  DIFFERENT FUNCTIONS MAY BE GIVEN FOR EACH REGION OF THE MESH.    
! ***  REGIONS ARE SPECIFIED BY BOTH NODE NUMBER AND ELEMENT NUMBER     
! ***  IN INPUT DATA FILE FOR UNIT fINP.                                  
!                                                                       
      SUBROUTINE UNSAT (SW, DSWDP, RELK, PRES, KREG) 
      
      USE CONTRL
      USE SutraMSPrecision

      IMPLICIT NONE
!                                                                       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
!     E X A M P L E   C O D I N G   FOR                                 
!     MESH WITH TWO REGIONS OF UNSATURATED PROPERTIES USING             
!     THREE PARAMETER-UNSATURATED FLOW RELATIONSHIPS OF                 
!     VAN GENUCHTEN(1980)                                               
!        RESIDUAL SATURATION, SWRES, GIVEN IN UNITS {L**0}              
!        PARAMETER, AA, GIVEN IN INVERSE PRESSURE UNITS {m*(s**2)/kg}   
!        PARAMETER, VN, GIVEN IN UNITS {L**0}                           
!                                                                       
      INTEGER (I4B) :: &
        KREG
      REAL (DP) :: &
        SW, DSWDP, RELK, PRES

      !LOCAL VARIABLES
      REAL (DP) :: &
        SWRES, AA, VN, SWRM1, AAPVN, VNF, AAPVNN, DNUM, DNOM, SWSTAR 
      REAL (DP) :: &
        SWRES1, SWRES2, AA1, AA2, VN1, VN2 
!                                                                       
!     DATA FOR REGION 1:                                                
      DATA SWRES1 / 0.30E0 /, AA1 / 5.0E-5 /, VN1 / 2.0E0 / 
      SAVE SWRES1, AA1, VN1 
!     DATA FOR REGION 2:                                                
      DATA SWRES2 / 0.30E0 /, AA2 / 5.0E-5 /, VN2 / 2.0E0 / 
      SAVE SWRES2, AA2, VN2 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
!                                                                       
! *** BECAUSE THIS ROUTINE IS CALLED OFTEN FOR UNSATURATED FLOW RUNS,   
! *** EXECUTION TIME MAY BE SAVED BY CAREFUL CODING OF DESIRED          
! *** RELATIONSHIPS USING ONLY INTEGER AND SINGLE PRECISION VARIABLES!  
! *** RESULTS OF THE CALCULATIONS MUST THEN BE PLACED INTO DOUBLE       
! *** PRECISION VARIABLES SW, DSWDP AND RELK BEFORE LEAVING             
! *** THIS SUBROUTINE.                                                  
!                                                                       
!                                                                       
!***********************************************************************
!***********************************************************************
!                                                                       
!     SET PARAMETERS FOR CURRENT REGION, KREG                           
      GOTO (10, 20), KREG 
   10 SWRES = SWRES1 
      AA = AA1 
      VN = VN1 
      GOTO 100 
   20 SWRES = SWRES2 
      AA = AA2 
      VN = VN2 
  100 CONTINUE 
!                                                                       
!                                                                       
!***********************************************************************
!***********************************************************************
!.....SECTION (1):                                                      
!     SW VS. PRES   (VALUE CALCULATED ON EACH CALL TO UNSAT)            
!     CODING MUST GIVE A VALUE TO SATURATION, SW.                       
!                                                                       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     THREE PARAMETER MODEL OF VAN GENUCHTEN(1980)                      
      SWRM1 = 1.E0 - SWRES 
      AAPVN = 1.E0 + (AA * ( - PRES) ) **VN 
      VNF = (VN - 1.E0) / VN 
      AAPVNN = AAPVN**VNF 
      SW = DBLE (SWRES + SWRM1 / AAPVNN) 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!***********************************************************************
!***********************************************************************
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      IF (IUNSAT - 2) 600, 1200, 1800 
!***********************************************************************
!***********************************************************************
!.....SECTION (2):                                                      
!     DSWDP VS. PRES, OR DSWDP VS. SW   (CALCULATED ONLY WHEN IUNSAT=1) 
!     CODING MUST GIVE A VALUE TO DERIVATIVE OF SATURATION WITH         
!     RESPECT TO PRESSURE, DSWDP.                                       
!                                                                       
  600 CONTINUE 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      DNUM = AA * (VN - 1.E0) * SWRM1 * (AA * ( - PRES) ) ** (VN - 1.E0) 
      DNOM = AAPVN * AAPVNN 
      DSWDP = DBLE (DNUM / DNOM) 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      GOTO 1800 
!***********************************************************************
!***********************************************************************
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!***********************************************************************
!***********************************************************************
!.....SECTION (3):                                                      
!     RELK VS. P, OR RELK VS. SW   (CALCULATED ONLY WHEN IUNSAT=2)      
!     CODING MUST GIVE A VALUE TO RELATIVE PERMEABILITY, RELK.          
!                                                                       
 1200 CONTINUE 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     GENERAL RELATIVE PERMEABILITY MODEL FROM VAN GENUCHTEN(1980)      
      SWSTAR = (SW - SWRES) / SWRM1 
      RELK = DBLE (SQRT (SWSTAR) * (1.E0 - (1.E0 - SWSTAR** (1.E0 / VNF)) ** (VNF) ) **2)                                                 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!***********************************************************************
!***********************************************************************
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
 1800 RETURN 
!                                                                       
      END SUBROUTINE UNSAT                          
