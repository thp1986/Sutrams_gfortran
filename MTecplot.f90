!
      MODULE  TECPLOT
        use DIMS
        use DIMX
        use TIMES
        use PARAMS 
        use OBS
        use CONTRL
        use FUNITS
        use MSErrorHandler
        USE SutraMSPrecision

        implicit none

        LOGICAL                           :: TecplotAll=.false.
        LOGICAL                           :: ONCETPN=.false.,ONCETPE=.false.
        LOGICAL                           :: LTpP,LTpU,LTpSW,LTpRho,LTpTotalHead,LTpVelocity, &
                                             LTpPrintVelocity, LTpPrintVelocityAtNodes
        LOGICAL                           :: lInitializeBinary = .false.
        LOGICAL                           :: lObservationsWritten = .false.
        character (len=80)                :: cTPNode, cTPElem, cTPObs
        INTEGER (I4B)                     :: NTPNP=0,NTPEP=0
        INTEGER (I4B)                     :: NodeFileHandle=0,ElementFileHandle=0,OnHandle=0
        INTEGER (I4B), allocatable        :: iConnection(:,:)
        real (DP), allocatable :: &
          VEXe(:), VEYe(:), VEZe(:) 
        real (DP), allocatable :: &
          VLXe(:), VLYe(:), VLZe(:) 
        real (DP), allocatable :: &
          TpNodeVx(:), TpNodeVy(:), TpNodeVz(:) 

        public :: &
               TecplotAll, &
               NTPNP,NTPEP, &
               cTPNode, cTPElem, cTPObs, &
               LTpP,LTpU,LTpSW,LTpRho,LTpTotalHead,LTpVelocity,LTpPrintVelocity, LTpPrintVelocityAtNodes, &
                                                         lObservationsWritten, &
               !for Tecplot Element Data
               AllocateTecplot, &
               PrintTpVelocity, &         !logical function to limit calls to MakeVelocityVector subroutine
               MakeVelocityVector, &      !generic subroutine for setting velocity vector
               PrintTecplotNodeData, &    !generic function for output of SUTRA element data
               PrintTecplotElementData    !generic function for output of SUTRA element data

        contains
          !function to allocate storage for Tecplot arrays
          logical function AllocateTecplot()
            use SutraStorage, ONLY: IN
            use TotalStorage
            implicit none
            integer (I4B) :: &
              II, III, &
              L

            lOk=.false.

            !allocate Tecplot velocity arrays
            if (fTPE>0) then
              if(.not.allocated(VEXe)) then
                allocate (VEXe(NE), VEYe(NE), VEZe(NE), &
                          VLXe(NE), VLYe(NE), VLZe(NE), &
                          stat=ios)
                if (ios/=0) call ErrorIO('Tecplot::AllocateTecplot Error allocating nodal arrays for Tecplot output')
                !Calculate memory requirements
                ios=AddMemory(MemoryIndex('TPL'),'Vec',NE) !VEXe(NE)
                ios=AddMemory(MemoryIndex('TPL'),'Vec',NE) !VEYe(NE)
                ios=AddMemory(MemoryIndex('TPL'),'Vec',NE) !VEZe(NE)
                ios=AddMemory(MemoryIndex('TPL'),'Vec',NE) !VLXe(NE)
                ios=AddMemory(MemoryIndex('TPL'),'Vec',NE) !VLYe(NE)
                ios=AddMemory(MemoryIndex('TPL'),'Vec',NE) !VLZe(NE)
              end if
            end if
            !allocate space of Tecplot velocity at nodes
            if (LTpVelocity) then
              if(.not.allocated(TPNodeVx)) then
                allocate (TpNodeVx(NN), TpNodeVy(NN), stat=ios)
                if (ios/=0) &
                  call ErrorIO('Tecplot::AllocateTecplot Error allocating nodal '//&
                             & 'x- and y-velocity arrays for Tecplot output')
                !Calculate memory requirements for x- and y- nodal velocity arrays
                ios=AddMemory(MemoryIndex('TPL'),'Vec',NN) !TpNodeVx(NN)
                ios=AddMemory(MemoryIndex('TPL'),'Vec',NN) !TpNodeVy(NN)
                if (N48==8) then
                  allocate (TpNodeVz(NN), stat=ios)
                  if (ios/=0) &
                    call ErrorIO('Tecplot::AllocateTecplot Error allocating nodal '//&
                               & 'z-velocity arrays for Tecplot output')
                  !Calculate memory requirements for z-nodal velocity array
                  ios=AddMemory(MemoryIndex('TPL'),'Vec',NN) !TpNodeVz(NN)
                end if
              end if
            end if

            !allocate storage for temporary Tecplot array for error output and
            !  allocate iconnection and initialize for Tecplot node files
            if(.not.allocated(iConnection)) &
              allocate(iConnection(N48,NE), &
                       stat=ios)
            if(ios/=0) &
              call ErrorIO('Tecplot::AllocateTecplot Error allocating TempTecplot '//&
                         & 'value for errors and iConnection for nodal output files')
            !Calculate memory requirements
            ios=AddMemory(MemoryIndex('TPL'),'Irr',N48*NE) !iConnection(N48,NE)
            iConnection=0
            do L = 1, NE 
              do II = 1, N48 
                III = II + (L - 1) * N48 
                iConnection (II,L) = IN (III) 
              end do 
            end do 

            lOk=.true.

            AllocateTecplot=lOk
            return
          end function AllocateTecplot



          logical function PrintTpVelocity()
            use UserSpecifiedOutputTime, only : DPTIME
            implicit none
            real (DP) :: EPSILON=1E-10
            LTpPrintVelocity=.FALSE.
            if ( (IT <= 1) .OR. (IT == ITMAX) ) goto 9998
            if (fOTM < 1) then
              !TECPLOT PRINTING                                                  
              IF (NTPNP.LT.1) NTPNP = ITMAX 
              if(MOD (IT, NTPNP) .EQ.0) goto 9998
            else
               if (ABS(TSEC-DPTIME).LT.EPSILON) goto 9998
            end if
            goto 9999
09998&
            LTpPrintVelocity=.TRUE.
09999&
            PrintTpVelocity=.TRUE.
            return
          end function PrintTpVelocity


          !function to fill nodal velocity vectors
          subroutine MakeVelocityVector(Node,Vx,Vy,Vz)
            implicit none
            integer (I4B), intent(in)              :: Node
            real (DP),     intent(in)              :: Vx
            real (DP),     intent(in)              :: Vy
            real (DP),     intent(in), optional    :: Vz
            lOk=.false.
            TPNodeVx(Node)=Vx
            TPNodeVy(Node)=Vy
            if( present(Vz) ) TPNodeVz(Node)=Vz
            lOk=.true.

            return
          end subroutine MakeVelocityVector

!...........FUNCTION TO OUTPUT NODE DATA TO TECPLOT FILES
          logical function PrintTecplotNodeData() 
            USE SutraStorage, ONLY: PVEC, UVEC, SW, RHO, IN, X, Y, Z
            IMPLICIT NONE
!.............locals
            CHARACTER (LEN = 128) :: CTPVAR 
            CHARACTER (LEN =  80) :: COUTF, CTITLE 
            CHARACTER (LEN =  18) :: CTPV1 
            CHARACTER (LEN =  10) :: CINT1, CINT2 
            CHARACTER (LEN =  20) :: cDupList
            INTEGER (I4B)         :: nElementType
            INTEGER (I4B) :: &
              I, II, &
              J, &
              K, &
              L, &
              N
            INTEGER (I4B)          :: IUNIT
            REAL (DP), ALLOCATABLE :: TotalHead(:)

!.............OUTPUT FORMAT STATEMENTS
        100 FORMAT('VARIABLES = "X", "Y"',A) 
        105 FORMAT(',"',A,'"') 
        110 FORMAT('ZONE T="N_',A,'" N=',A,', E=',A,', F=FEBLOCK, ET=QUADRILATERAL')
        120 FORMAT('ZONE T="N_',A,'" N=',A,', E=',A,', F=FEBLOCK, ET=BRICK')
        130 FORMAT('ZONE T="N_',A,'" N=',A,', E=',A,', F=FEBLOCK, ET=QUADRILATERAL, D=(1,2,FECONNECT)')
        140 FORMAT('ZONE T="N_',A,'" N=',A,', E=',A,', F=FEBLOCK, ET=BRICK, D=(1,2,3,FECONNECT) ')
        200 FORMAT(10(1PE15.7,1X)) 
        300 FORMAT(8(I10,1X)) 

!.............CODE
            lOk   = .false.
            MSErrorValue%cDataSet = 'TPN'
            IUNIT = fTPN 
            IF (IUNIT.LT.1) GOTO 9999 
!.............WRITE BLANK LINE TO SEPARATE DATA
            WRITE (IUNIT, * ) ' ' 
!.............CREATE TEXT EQUIVALENTS OF PROBLEM DIMENSIONS
            WRITE (CINT1, '(I10)') NN 
            WRITE (CINT2, '(I10)') NE 
!.............WRITE APPROPRIATE HEADING                                         
            IF (ONCETPN) THEN 
               IF (N48.EQ.4) THEN 
                  WRITE (IUNIT, 130) TRIM (ADJUSTL (cTime()) ), CINT1, CINT2 
               ELSE 
                  WRITE (IUNIT, 140) TRIM (ADJUSTL (cTime()) ), CINT1, CINT2 
               ENDIF 
            ELSE 
!................Generate list of variables to write to Tecplot node file
               CTPVAR = ''
               CTPV1 = '' 
               IF (N48.EQ.8) THEN
                 WRITE (CTPV1, 105) 'Z' 
                 CTPVAR = TRIM (ADJUSTL (CTPVAR) ) //TRIM (ADJUSTL (CTPV1) ) 
               ENDIF 
               IF (LTpP) THEN 
                 WRITE (CTPV1, 105) 'P' 
                 CTPVAR = TRIM (ADJUSTL (CTPVAR) ) //TRIM (ADJUSTL (CTPV1) ) 
               ENDIF 
               IF (LTpU) THEN 
                 DO K = 1, NSPE 
                   WRITE (CTPV1, 105) TRIM (ADJUSTL (SPNAME (K) ) ) 
                   CTPVAR = TRIM (ADJUSTL (CTPVAR) ) //TRIM (ADJUSTL (CTPV1) ) 
                 ENDDO 
               ENDIF 
               IF (LTpRho) THEN
                 WRITE (CTPV1, 105) 'Rho' 
                 CTPVAR = TRIM (ADJUSTL (CTPVAR) ) //TRIM (ADJUSTL (CTPV1) ) 
               ENDIF 
               IF (LTpTotalHead) THEN
                 WRITE (CTPV1, 105) 'TotalHead' 
                 CTPVAR = TRIM (ADJUSTL (CTPVAR) ) //TRIM (ADJUSTL (CTPV1) ) 
               ENDIF 
               IF (LTpSW) THEN 
                 WRITE (CTPV1, 105) 'SW' 
                 CTPVAR = TRIM (ADJUSTL (CTPVAR) ) //TRIM (ADJUSTL (CTPV1) ) 
               ENDIF 
               IF (LTpVelocity) then
                 WRITE (CTPV1, 105) 'VX' 
                 CTPVAR = TRIM (ADJUSTL (CTPVAR) ) //TRIM (ADJUSTL (CTPV1) ) 
                 WRITE (CTPV1, 105) 'VY' 
                 CTPVAR = TRIM (ADJUSTL (CTPVAR) ) //TRIM (ADJUSTL (CTPV1) ) 
                 IF (N48 == 8) THEN
                   WRITE (CTPV1, 105) 'VZ' 
                   CTPVAR = TRIM (ADJUSTL (CTPVAR) ) //TRIM (ADJUSTL (CTPV1) ) 
                 END IF
               ENDIF 

!...............GENERATE APPROPRIATE TITLE for ASCII file
              CTITLE = 'TITLE = "SUTRA: FE VOLUME QUADRILATERAL DATA"' 
              if (N48==8) CTITLE = 'TITLE = "SUTRA: FE VOLUME BRICK DATA"'
!...............Write title and list of variables to ASCII file
              WRITE (IUNIT, '(1X,A)') CTITLE 
              WRITE (IUNIT, 100) TRIM (ADJUSTL (CTPVAR) ) 
!...............Write header for zone to ASCII file
              IF (N48.EQ.4) THEN 
                 WRITE (IUNIT, 110) TRIM (ADJUSTL (CTIME()) ), CINT1, CINT2 
              ELSE 
                 WRITE (IUNIT, 120) TRIM (ADJUSTL (CTIME()) ), CINT1, CINT2 
              ENDIF 
            ENDIF 
            
!.............WRITE DATA LINES TO FILE                                          
!.............COORDINATES - STATIC DATA
            IF (ONCETPN) GOTO 1000 
            WRITE (IUNIT, 200) (X (I), I = 1, NN) 
            WRITE (IUNIT, 200) (Y (I), I = 1, NN) 
            IF (N48.EQ.8) WRITE (IUNIT, 200) (Z (I), I = 1, NN) 
!.............TRANSIENT DATA                                                    
     1000   CONTINUE 
            IF (LTpP) THEN 
               WRITE (IUNIT, 200) (PVEC (I), I = 1, NN) 
            ENDIF 
            IF (LTpU) THEN 
               DO K = 1, NSPE 
                  WRITE (IUNIT, 200) (UVEC (I, K), I = 1, NN) 
               ENDDO 
            ENDIF 
            IF (LTpRho) THEN 
               WRITE (IUNIT, 200) (RHO (I), I = 1, NN) 
            ENDIF 
            IF (LTpTotalHead) THEN 
               allocate (TotalHead(NN),stat=ios)
               if(ios/=0)  &
                 call ErrorIO('Tecplot::PrintTecplotNodeData Error allocating real '//&
                            & 'vector for single-precision total head Tecplot output')
               if(.not.CalculateTotalHead(NN,TotalHead)) goto 9999
               WRITE (IUNIT, 200) (TotalHead(I), I = 1, NN) 
               deallocate (TotalHead,stat=ios)
               if(ios/=0)  &
                 call ErrorIO('Tecplot::PrintTecplotNodeData Error deallocating real '//&
                            & 'vector for single-precision total head Tecplot output')
            ENDIF 
            IF (LTpSW) THEN 
               WRITE (IUNIT, 200) (SW (I), I = 1, NN) 
            ENDIF 
            IF (LTpVelocity) THEN 
               WRITE (IUNIT, 200) (TpNodeVx(I), I = 1, NN) 
               WRITE (IUNIT, 200) (TpNodeVy(I), I = 1, NN) 
               IF (N48==8) WRITE (IUNIT, 200) (TpNodeVz(I), I = 1, NN) 
            ENDIF 
!.............WRITE CONNECTIVITY TO DATA FILE                                   
            IF (ONCETPN) GOTO 9999 
!.............WRITE BLANK LINE TO SEPARATE DATA                                 
            WRITE (IUNIT, * ) ' ' 
            DO L = 1, NE 
              WRITE (IUNIT, 300) (iConnection (II,L), II = 1, N48) 
            ENDDO 
!.............SET LOGICAL FLAG TO INDICATE INITIAL NODAL DATA WRITTEN
            IF (.NOT.TecplotAll) ONCETPN = .TRUE.
!
!.............RETURN TO CALLING ROUTINE                                         
9999 &
            lOk=.true.
            PrintTecplotNodeData = lOk
            RETURN 
          end function PrintTecplotNodeData


!...........FUNCTION TO OUTPUT ELEMENT DATA TO TECPLOT FILES
          logical function PrintTecplotElementData() 
            USE SutraStorage, ONLY : VMAG, VANG1, VANG2, IN, X, Y, Z
            IMPLICIT NONE
!.............LOCAL VARIABLES
            CHARACTER (LEN = 80) :: COUTF, CTITLE 
            CHARACTER (LEN = 10) :: CINT1, CINT2 
            INTEGER (I4B)        :: IINe(8)
            INTEGER (I4B)        :: IUNIT, L, LLe, IIe, IIIe
            REAL (DP) :: &
              DKTM2e, RN48e, &
              CENTRXe, CENTRYe, CENTRZe, &
              VECTRXe, VECTRYe, VECTRZe, &
              VA1e, VA2e, CVA2e 
            CHARACTER (LEN = 128) :: CTPVAR 
            CHARACTER (LEN =  20) :: cDupList
            INTEGER (I4B) :: nElementType

!.............OUTPUT FORMATS
        100 FORMAT('VARIABLES = "',A,'"',5(:',',1X,'"',A,'"')) 
        110 FORMAT('ZONE T="VE_',A,'" I=',A,', F=BLOCK')
        120 FORMAT('ZONE T="VE_',A,'" I=',A,', F=BLOCK')
        130 FORMAT('ZONE T="VE_',A,'" I=',A,', F=BLOCK, D=(1,2)')
        140 FORMAT('ZONE T="VE_',A,'" I=',A,', F=BLOCK, D=(1,2,3) ')
        200 FORMAT(10(1PE15.7,1X)) 

!.............CODE
            lOK=.false.
            MSErrorValue%cDataSet='TPE'
            IUNIT = fTPE 
!.............RETURN IF NOT CREATING THE ELEMENTAL TECPLOT FILE
            IF (IUNIT.LT.1) GOTO 9999 

            DKTM2e = DBLE (IABS (KTYPE) - 2) 
            VEXe = 0.0D+00 
            VEYe = 0.0D+00 
            VEZe = 0.0D+00 
            VLXe = 0.0D+00 
            VLYe = 0.0D+00 
            VLZe = 0.0D+00 
!.............WRITE BLANK LINE TO SEPARATE DATA
            WRITE (IUNIT, * ) ' ' 
!.............CREATE TEXT EQUIVALENT OF PROBLEM DIMENSIONS
            WRITE (CINT1, '(I10)') NE 
!.............WRITE APPROPRIATE HEADING
            IF (ONCETPE) THEN 
              IF (N48.EQ.4) THEN 
                WRITE (IUNIT, 130) TRIM (ADJUSTL (cTime()) ), CINT1 
              ELSE 
                WRITE (IUNIT, 140) TRIM (ADJUSTL (cTime()) ), CINT1 
              ENDIF 
            ELSE 
!................GENERATE TITLE
              CTITLE = 'TITLE = "SUTRA: POINT VELOCITY DATA"' 
              WRITE (IUNIT, '(1X,A)') CTITLE 
              IF (N48.EQ.4) THEN 
                WRITE (IUNIT, 100) 'X', 'Y', 'VX', 'VY' 
                WRITE (IUNIT, 110) TRIM (ADJUSTL (cTime()) ), CINT1 
              ELSE 
                WRITE (IUNIT, 100) 'X', 'Y', 'Z', 'VX', 'VY', 'VZ' 
                WRITE (IUNIT, 120) TRIM (ADJUSTL (cTime()) ), CINT1 
              ENDIF 
            ENDIF 

!.............VELOCITY CALCULATION ROUTINE FROM OUTK6                           
!.............THE VELOCITY DATA FOR THIS TIME STEP                              
            RN48e = 1D0 / DBLE (N48) 
            DO 50 L = 1, NE 
              IF (ONCETPE) GOTO 45 
               CENTRXe = 0D0 
               CENTRYe = 0D0 
               CENTRZe = 0D0 
               DO 40 IIe = 1, N48 
                  IIIe = IIe + (L - 1) * N48 
                  IINe (IIe) = IN (IIIe) 
                  CENTRXe = CENTRXe + X (IINe (IIe) ) 
                  CENTRYe = CENTRYe + Y (IINe (IIe) ) 
                  CENTRZe = CENTRZe + Z (IINe (IIe) ) 
      00040    CONTINUE
               CENTRXe = CENTRXe * RN48e 
               CENTRYe = CENTRYe * RN48e 
               CENTRZe = CENTRZe * RN48e 
!................SPECIFIC TO OUTTPE SUBROUTINE (COORDINATES)
               VEXe (L) = CENTRXe 
               VEYe (L) = CENTRYe 
               VEZe (L) = CENTRZe 
!................RETURN TO OUTK6 CODE
      00045    VA1e = 0.017453292D0 * VANG1 (L) 
               LLe = MIN(L, NEX) 
               VA2e = 0.017453292D0 * VANG2 (LLe) * DKTM2e 
               CVA2e = DCOS (VA2e) 
               VECTRXe = VMAG (L) * DCOS (VA1e) * CVA2e 
               VECTRYe = VMAG (L) * DSIN (VA1e) * CVA2e 
               VECTRZe = VMAG (L) * DSIN (VA2e) 
!................SPECIFIC TO OUTTPE SUBROUTINE (COORDINATES)
               VLXe (L) = VECTRXe 
               VLYe (L) = VECTRYe 
               VLZe (L) = VECTRZe 
!.............END OF ELEMENT LOOP
      00050 CONTINUE
!.............WRITE DATA LINES TO FILE                                          
!.............COORDINATES - STATIC DATA                                         
            IF (ONCETPE) GOTO 1000 
            WRITE (IUNIT, 200) (VEXe (L), L = 1, NE) 
            WRITE (IUNIT, 200) (VEYe (L), L = 1, NE) 
            IF (N48.EQ.8) WRITE (IUNIT, 200) (VEZe (L), L = 1, NE) 
!.............VELOCITY DATA                                                     
       1000 CONTINUE 
            WRITE (IUNIT, 200) (VLXe (L), L = 1, NE) 
            WRITE (IUNIT, 200) (VLYe (L), L = 1, NE) 
            IF (N48.EQ.8) WRITE (IUNIT, 200) (VLZe (L), L = 1, NE) 
!.............SET LOGICAL FLAG TO INDICATE INITIAL NODAL DATA WRITTEN
            IF (.NOT.TecplotAll) ONCETPE = .TRUE. 
            lOk=.true.
!.............RETURN TO CALLING ROUTINE                                         
       9999 PrintTecplotElementData=lOk
            RETURN 
          end function PrintTecplotElementData



!..........PRIVATE FUNCTIONS
!..........GENERIC ROUTINE FOR CALCULATING THE TOTAL HEAD
          logical function CalculateTotalHead(NN,TotalHead)
            use GRAVEC
            use SutraStorage, only : Y, Z, RHO, PVEC, UVEC
            IMPLICIT NONE
!.............PASSED VARIABLES            
            integer (I4B), intent(in)               :: NN
            real (DP), intent(inout), dimension(NN) :: TotalHead
!.............LOCAL VARIABLES            
            integer (I4B) :: k, n
            real (DP) :: Gravity
            real (DP) :: ElevationHead
            real (DP) :: RhoV
            CalculateTotalHead=.false.
            Gravity = -GRAVY
            if (N48 == 8) Gravity = -GRAVZ
            do n=1,NN
              select case (N48)
                case (4)
                  ElevationHead = Y(n)
                case (8)
                  ElevationHead = Z(n)
              end select
              if (TSEC<=TSTART) then
                RhoV = RHOW0 
                do k = 1, nspe 
                  RhoV = RhoV + DRWDU(k) * ( UVEC(n,k) - URHOW0(k) )                                                             
                end do 
              else
                RhoV = RHO(n)
              end if
              TotalHead(n)=( PVEC(n) / (RhoV * Gravity) ) + ElevationHead
            end do
            CalculateTotalHead=.true.
!.............RETURN TO CALLING ROUTINE
09999&
            return
          end function CalculateTotalHead

      END MODULE  TECPLOT
