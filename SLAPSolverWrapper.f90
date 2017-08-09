
!     SUBROUTINE        S  L  P  W  R  P       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO SERVE AS A "WRAPPER" FOR THE SLAP SOLVERS, PERFORMING         
! ***  SOME PRELIMINARIES ON VECTORS AND MATRIX POINTERS BEFORE         
! ***  CALLING A SOLVER.                                                
!                                                                       

      LOGICAL FUNCTION SLPWRP(KKK, KPU, KSOLVR, A, R, XITER, B, NNP, IHALFB, &
                              MAXNP, MAXBW)
        USE SOLVI 
        USE ITSOLR 
        USE ITSOLI 
        USE PARAMS 
        USE FUNITS 
        USE ATSDATA 
        USE DIMX
        USE MSErrorHandler
        USE Triad2Column
        USE SutraStorage, ONLY : IWK, FWK, IPARM, FPARM, IA=>ITRI, JA=>JTRI
        use SutraMSPrecision
        USE ColumnStorage

        IMPLICIT NONE

        LOGICAL :: lrc 
        INTEGER (I4B) :: &
          NNP, &
          MAXNP, MAXBW, &
          IHALFB, &
          KKK, KPU, KSOLVR
        REAL (DP) :: &
          A (NELT) 
        REAL (DP) :: &
          XITER (NNP), R (NNP)
        REAL (DP) :: &
          B (NNNX) 

        !LOCAL VARIABLES
        CHARACTER(1) KPUTXT (2) 
        INTEGER (I4B) :: &
          K, N, &
          IUNIT, ISYM, ITRMX, ITOL, ITER, IERR, &
          NSAVE
        REAL (DP) :: &
          ERR, TOL

        !DATA
        DATA (KPUTXT (K) , K = 1, 2) / 'P', 'U' / 
        SAVE KPUTXT 
        EXTERNAL DSICCG, DSLUGM 
!
!
        MSErrorValue%cDataSet='SLW'
!
        lrc=.TRUE.
!                                                                       
!.......RESET MATRIX ARRAY POINTERS IF NECESSARY.                         
!                                                                       
        IF (KKK.EQ.0) then
        
          if( .not.lColumnStorage ) then
            !SET IA AND JA TO SAVED VALUES OF nIA AND nJA
            IA=nIA
            JA=nJA

            if(.not.MapTriad2Column(NELT,A)) call ErrorIO('SLPWRP:: Error mapping A from SLAP Triad format to SLAP Column format')

            if(lDebugData) then
              write (fLST,'(//a)') 'Column Format - Triad Converted to Column using SLAP Routines'
              write (fLST,'(a7,20(1x,i6),1000(:/7x,20(1x,i6)))') &
                'SLAP ia',ia(1:nelt)
              write (fLST,'(a7,20(1x,i6),1000(:/7x,20(1x,i6)))') &
                'SLAP ja',ja(1:nelt)
              write (fLST,'(a7,20(1x,i6),1000(:/7x,20(1x,i6)))') &
                'SLAP a', a(1:nelt)
            end if
          end if
        end if
!                                                                       
!.......COPY THE RHS VECTOR "R" INTO "B", THEN USE "R" AS THE SLAP        
!          SOLUTION VECTOR.  INITIALIZE IT FROM THE LATEST SUTRA SOLUTION.
!          WE'RE NOT USING XITER AS THE SOLUTION VECTOR BECAUSE IT MIGHT  
!          MESS UP SUBSEQUENT SUTRA CALCULATIONS IF WE CHANGE IT.         
!                                                                       
        forall (n = 1:nnp) 
          B (n) = R (n) 
          R (n) = XITER (n) 
        endforall 
!                                                                       
!.......SET ITERATIVE SOLVER PARAMETERS.                                  
!          IUNIT --> UNIT # ON WHICH TO WRITE SOLVER ERROR (0 = NONE)     
!          ISYM  --> 0 = FULL STORAGE; 1 = SYMMETRIC STORAGE              
!          ITRMX --> MAXIMUM # OF SOLVER ITERATIONS                       
!          ITOL  --> SPECIFIES TYPE OF CONVERGENCE CRITERION              
!          TOL   --> CONVERGENCE TOLERANCE                                
!          NSAVE --> # OF DIRECTION VECTORS (USED BY GMRES & ORTHOMIN)    
!                                                                       
        IF (KPU.EQ.1) THEN 
!..........SET PARAMETERS FOR ITERATIVE P SOLUTION.                       
           IUNIT = 0 
           ISYM = 0 
           ITRMX = ITRMXP 
           ITOL = ITOLP 
           TOL = TOLP 
           NSAVE = NSAVEP 
       ELSE 
!..........SET PARAMETERS FOR ITERATIVE U SOLUTION.                       
           IUNIT = 0 
           ISYM = 0 
           ITRMX = ITRMXU (KSP) 
           ITOL = ITOLU (KSP) 
           TOL = TOLU (KSP) 
           NSAVE = NSAVEU (KSP) 
        ENDIF 
!                                                                       
        IUNIT = fSMY 
!                                                                       
!.......SET ITER {NUMBER OF SOLVER ITERATIONS} TO 1 INITIALLY             
        ITER = 1 
!                                                                       
!.......CALL A SLAP SOLVER:                                               
!          DSICCG = CG WITH IC PRECONDITIONING,                           
!          DSLUGM = GMRES WITH ILU PRECONDITIONING,                       
!          DSLUOM = ORTHOMIN WITH ILU PRECONDITIONING.                    
!                                                                       
        IF (KSOLVR.EQ.1) THEN 
           CALL DSICCG (NNP, B, R, NELT, IA, JA, A, ISYM, ITOL, TOL,      &
           ITRMX, ITER, ERR, IERR, IUNIT, FWK, NWF, IWK, NWI)             
        ELSEIF (KSOLVR.EQ.2) THEN 
           CALL DSLUGM (NNP, B, R, NELT, IA, JA, A, ISYM, NSAVE, ITOL,    &
           TOL, ITRMX, ITER, ERR, IERR, IUNIT, FWK, NWF, IWK, NWI)        
        ELSE 
           CALL DSLUOM (NNP, B, R, NELT, IA, JA, A, ISYM, NSAVE, ITOL,    &
           TOL, ITRMX, ITER, ERR, IERR, IUNIT, FWK, NWF, IWK, NWI)        
        ENDIF 

!                                                                       
      !ATSDATA MODULE            
      IF (iter>imaxiter) imaxiter = iter 
!
      IF (IERR.EQ.0) THEN 
         WRITE ( *  , 555) KPUTXT (KPU), ITER, ERR 
         WRITE (fSMY, 555) KPUTXT (KPU), ITER, ERR 
  555 FORMAT    (1X,6x,A1, '-solution converged ', I5, ' iters  (Error ~ ', 1PE8.1, ')')                           
      ELSE 
         WRITE ( *  , 557) KPUTXT (KPU), IERR, ITER, ERR 
         WRITE (fSMY, 557) KPUTXT (KPU), IERR, ITER, ERR 
  557 FORMAT    (                                                       &
     &      //3X,62('*'),                                               &
     &      /5X, A1, '-MATRIX SOLUTION TERMINATED ',                    &
     &      'WITH A SLAP SOLVER ERROR:'                                 &
     &      /8X, 'Error flag.....................IERR = ', I3,          &
     &      /8X, 'Number of solver iterations....ITER = ', I5,          &
     &      /8X, 'Error estimate..................ERR = ', 1PE8.1,      &
     &      /5X, '(For details see documentation in SLAP solver ',      &
     &      'source code)',                                             &
     &      /3X,62('*')// )                                             
         WRITE (fLST, 558) 
  558 FORMAT   (////////11X,'SIMULATION TERMINATED DUE TO ',            &
     &      'ERROR DURING MATRIX SOLUTION',                             &
     &      /11X,'********** ********** *** ** ',                       &
     &      '***** ****** ****** ********')                             
         lrc=.FALSE.
        ENDIF 
        SLPWRP=lrc                                                                        
      RETURN 
      END FUNCTION SLPWRP                         
