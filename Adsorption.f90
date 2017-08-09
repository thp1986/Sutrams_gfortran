!     SUBROUTINE        A  D  S  O  R  B        SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO CALCULATE VALUES OF EQUILIBRIUM SORPTION PARAMETERS FOR       
! ***  LINEAR, FREUNDLICH, AND LANGMUIR MODELS.                         
!                                                                       
      SUBROUTINE ADSORB () 
        use SutraMSPrecision 
        USE PARAMS 
        USE MODSOR 
        USE DIMS
        USE SutraZoneModule
        use SutraStorage, ONLY : CS1, CS2, CS3, SL, SR, U=>UITER
        use MSErrorHandler
        IMPLICIT NONE
        !LOCALS
          INTEGER (I4B) :: &
            imap
          INTEGER (I4B) :: &
            I, &
            NK
!          REAL (DP) :: &
!            CHCH, &
!            DCHI2, &
!            RH2, &
!            CHI2F
          REAL (DP), DIMENSION(NN) :: &
            CHCH, &
            DCHI2, &
            RH2, &
            CHI2F
!                                                                       
!.....NOTE THAT THE CONCENTRATION OF ADSORBATE, CS(I), IS GIVEN BY:     
!     CS(I) = SL(I)*U(I) + SR(I)                                        
!                                                                       
      WRITE ( * , * ) 'CALCULATING ADSORPTION FOR SPECIES [',TRIM(ADJUSTL(Val2Char(KSP))), ']' 
!                                                                       
!.....NO SORPTION                                                       
      IF (trim(ADSMOD(KSP)) .NE.'NONE') GOTO 450 
!      forall (i = 1:nn, nk = ksp:ksp) 
!!        CS1 (i, nk) = 0.D0 
!      end forall 
!      forall (i = 1:nn) 
!        CS2 (i) = 0.D0 
!        CS3 (i) = 0.D0 
!        SL (i)  = 0.D0 
!        SR (i)  = 0.D0 
!      end forall 
      do i = 1, NN
        CS1(i, KSP) = 0.D0 
        CS2(i)      = 0.D0 
        CS3(i)      = 0.D0 
        SL(i)       = 0.D0 
        SR(i)       = 0.D0 
      end do
      GOTO 2000 
!                                                                       
!.....LINEAR SORPTION MODEL                                             
  450 IF (trim(ADSMOD(KSP)) .NE.'LINEAR') GOTO 700 
!      forall (i = 1:nn, nk = ksp:ksp) 
!!        CS1 (i, nk) = CHI1 (nk) * RHOW0 
!        !imap = NodeMap(i)
!        CS1 (i, nk) = ProdSorp(NodeMap(i))%chi1(nk) * RHOW0 
!      end forall 
!      forall (i = 1:nn) 
!        !imap = NodeMap(i)
!        CS2 (i) = 0.D0 
!        CS3 (i) = 0.D0 
!!        SL (i) = CHI1 (KSP) * RHOW0 
!        SL (i) = ProdSorp(NodeMap(i))%chi1(KSP) * RHOW0 
!        SR (i) = 0.D0 
!      end forall 
      do i = 1, NN
!        CS1 (i, nk) = CHI1 (nk) * RHOW0 
        imap        = NodeMap(i)
        CS1(i, KSP) = ProdSorp(imap)%chi1(KSP) * RHOW0 
        CS2(i)      = 0.D0 
        CS3(i)      = 0.D0 
!        SL(i)       = CHI1(KSP) * RHOW0 
        SL(i)       = ProdSorp(imap)%chi1(KSP) * RHOW0 
        SR(i)       = 0.D0 
      end do 
      GOTO 2000 
!                                                                       
!.....FREUNDLICH SORPTION MODEL                                         
  700 IF (trim(ADSMOD(KSP)) .NE.'FREUNDLICH') GOTO 950 
!      !static variables for 'FREUNDLICH'                                
!      CHCH = CHI1 (KSP) / CHI2 (KSP) 
!      DCHI2 = 1.D0 / CHI2 (KSP) 
!      RH2 = RHOW0**DCHI2 
!      CHI2F = ( (1.D0 - CHI2 (KSP) ) / CHI2 (KSP) )
      !variables for 'FREUNDLICH' - VERSION 1.1
      forall (i = 1:nn)
        !imap     = NodeMap(i)
        CHCH(i)  = ProdSorp(NodeMap(i))%chi1(KSP) / ProdSorp(NodeMap(i))%chi2(KSP) 
        DCHI2(i) = 1.D0 / ProdSorp(NodeMap(i))%chi2(KSP) 
        RH2(i)   = RHOW0**DCHI2(i) 
        CHI2F(i) = ( (1.D0 - ProdSorp(NodeMap(i))%chi2(KSP) ) / ProdSorp(NodeMap(i))%chi2(KSP) )
      end forall
!                                                                       
      forall (i = 1:nn) 
        CS2(i) = 0.D0 
        CS3(i) = 0.D0 
        SR(i)  = 0.D0 
      end forall 
      !U(i,ksp) <= 0.0D0 -----------------------------------------------
      forall (i = 1:nn, nk = ksp:ksp, U(i, nk) <= 0.0d0) 
!        CS1 (i, nk) = CHCH * RH2 * 1.0D0 
        CS1(i, nk) = CHCH(i) * RH2(i) * 1.0D0 
      end forall 
      forall (i = 1:nn, U(i, KSP) <= 0.0d0) 
!        SL (i) = CHI1 (KSP) * RH2 * 1.0D0 
        !imap  = NodeMap(i)
        SL(i) = ProdSorp(NodeMap(i))%chi1(KSP) * RH2(i) * 1.0D0 
      end forall 
      !U(i,ksp) <= 0.0D0 -----------------------------------------------
!                                                                       
      !U(i,ksp) >  0.0D0 -----------------------------------------------
      forall (i = 1:nn, nk = ksp:ksp, U (i, nk) >0.0d0) 
!        CS1 (i, nk) = CHCH * RH2 * (U (i, nk) **CHI2F) 
        CS1(i, nk) = CHCH(i) * RH2(i) * (U (i, nk)**CHI2F(i)) 
      end forall 
      forall (i = 1:nn, U (i, KSP) >0.0d0) 
!        SL (i) = CHI1 (KSP) * RH2 * (U (i, KSP) **CHI2F) 
        !imap  = NodeMap(i)
        SL(i) = ProdSorp(NodeMap(i))%chi1(KSP) * RH2(i) * (U(i, KSP)**CHI2F(i)) 
      end forall 
      !U(i,ksp) >  0.0D0 -----------------------------------------------
                                                                        
      GOTO 2000 
!                                                                       
!.....LANGMUIR SORPTION MODEL                                           
  950 IF (trim(ADSMOD(KSP)) .NE.'LANGMUIR') GOTO 2000 
      forall (i = 1:nn, nk = ksp:ksp) 
!        CS1 (i, nk) = (CHI1 (nk) * RHOW0) / ( (1.D0 + CHI2 (nk) * RHOW0 * &
!        U (i, nk) ) * (1.D0 + CHI2 (nk) * RHOW0 * U (i, nk) ) )           
        !imap  = NodeMap(i)
        CS1(i, nk) = (ProdSorp(NodeMap(i))%chi1(nk) * RHOW0) / ( (1.D0 + ProdSorp(NodeMap(i))%chi2(nk) * RHOW0 * &
                      U(i, nk) ) * (1.D0 + ProdSorp(NodeMap(i))%chi2(nk) * RHOW0 * U(i, nk) ) )           
      end forall 
      !results below depend on results above (CS1)  --------------------
      forall (i = 1:nn) 
        CS2(i) = 0.D0 
        CS3(i) = 0.D0 
        SL(i) = CS1(i, KSP) 
!        SR(i) = CS1(i, KSP) * CHI2(KSP) * RHOW0 * U(i, KSP) * U(i,KSP)                                                              
        !imap  = NodeMap(i)
        SR(i) = CS1(i, KSP) * ProdSorp(NodeMap(i))%chi2(KSP) * RHOW0 * U(i, KSP) * U(i,KSP)                                                              
      end forall 
      !results above depend on results above (CS1) --------------------
!                                                                       
!.....RETURN TO CALLING ROUTINE                                         
 2000 RETURN 
      END SUBROUTINE ADSORB                         
