!     MODULE            COMMONITEMS            SUTRA-MS VERSION 2004.1
!                                                                       
! ***  
! *** PURPOSE :                                                         
! ***  MODULES CONTAINING COMMON PROGRAM VARIABLES AND USED.
! ***  TO REPLACE COMMON STATEMENTS USED IN STANDARD
! ***  VERSION OF SUTRA.
! ***  
! ***  MODULES ALSO CONTAIN THE FOLLOWING SUBROUTINES:
! ***    CTIME -
! ***      FUNCTION THAT CONVERTS THE SIMULATION TIME (TSEC)
! ***      TO A NICE FORMAT TEXT STRING WITH UNITS OF
! ***      SEC, MIN, HR, DAYS, YEARS.  THE TIME UNITS ARE
! ***      CONCATENATED TO THE END OF THE STRING AND DEPEND
! ***      ON THE MAGNITUDE OF TSEC.
! ***    GENERICTIMER -
! ***      FUNCTION THAT ALLOWS SPECIFIC SUBROUTINES TO BE
! ***      TIMED.  USED TO DETERMIN NUMERICAL OVERHEAD ASSOCIATED
! ***      WITH SPECIFIC PORTIONS OF SUTRA-MS (I.E. INTEGRATION)
! ***    ALLOCATESTANDARDSUTRAOBSERVATIONS -
! ***      FUNCTION THAT ALLOCATED THE REQUIRED STORAGE FOR
! ***      THE STANDARD SUTRA OBSERVATION FUNCTIONALITY
! ***  
      MODULE ITERAT 
        USE SutraMSPrecision 
        REAL (DP) ::RPM, RPMAX 
        REAL (DP), ALLOCATABLE::RUMAX (:), RUM (:) 
        INTEGER (I4B) ::ITER, ITRMAX, IPWORS 
        INTEGER (I4B), ALLOCATABLE::IUWORS (:) 

        PUBLIC RPM, RPMAX, RUM, RUMAX, ITER, ITRMAX, IPWORS, IUWORS 
      END MODULE ITERAT 
!                                                                       
      MODULE PARAMS 
        USE SutraMSPrecision
        INTEGER (I4B), PARAMETER :: SPNameLen=20
        REAL (DP), PARAMETER :: SigmaAir=0.026  !W / (m degC)
        LOGICAL :: LVolAvgLambda = .true.
        CHARACTER (SPNameLen), ALLOCATABLE::SPNAME (:) 
        REAL (DP) ::COMPFL, COMPMA, CW, CS 
        REAL (DP) ::RHOS, SIGMAS 
        REAL (DP) ::RHOW0, BVISC0 
        REAL (DP), ALLOCATABLE::VISC0 (:) 
        REAL (DP), ALLOCATABLE::PRODF1 (:), PRODS1 (:) 
        REAL (DP), ALLOCATABLE::PRODF0 (:), PRODS0 (:) 
        REAL (DP), ALLOCATABLE::DRWDU (:), SIGMAW (:), URHOW0 (:) 
        REAL (DP), ALLOCATABLE::CHI1 (:), CHI2 (:) 
        REAL (DP), ALLOCATABLE::ATSPMULT (:) 
        INTEGER (I4B) :: NESP 
        INTEGER (I4B), ALLOCATABLE :: NSOU (:), NUBC (:) 
        INTEGER (I4B) :: KSP 
        PUBLIC :: &
                  SPNameLen, &
                  SigmaAir, &
                  SPNAME, COMPFL, COMPMA, CW, CS, RHOS, SIGMAS, RHOW0, BVISC0,&
                  VISC0, PRODF1, PRODS1, PRODF0, PRODS0, DRWDU, SIGMAW, URHOW0,     &
                  CHI1, CHI2, ATSPMULT, NESP, NSOU, NUBC, KSP, &
                  LVolAvgLambda
      END MODULE PARAMS 
!                                                                       
      MODULE CONTRL 
        USE SutraMSPrecision 
        REAL (DP) ::GNUP 
        REAL (DP), ALLOCATABLE::GNUU (:) 
        REAL (DP) ::UP, DTMULT, DTMAX 
        INTEGER (I4B) ::ME, ISSFLO, ISSTRA, ITCYC 
        INTEGER (I4B) ::NPCYC, NUCYC 
        INTEGER (I4B) ::NPRINT, IREAD, ISTORE, NOUMAT 
        INTEGER (I4B) ::IUNSAT, KTYPE 
        PUBLIC GNUP, GNUU, UP, DTMULT, DTMAX, ME, ISSFLO, ISSTRA, ITCYC,   &
        NPCYC, NUCYC, NPRINT, IREAD, ISTORE, NOUMAT, IUNSAT, KTYPE        
      END MODULE CONTRL 
!                                                                       
      MODULE SOLVI 
        USE SutraMSPrecision 
        INTEGER (I4B) :: KSOLVP, NN1, NN2, NN3 
        INTEGER (I4B), ALLOCATABLE::KSOLVU (:)
        INTEGER (I4B), ALLOCATABLE :: U_ML(:)
        REAL (DP), ALLOCATABLE :: StartUTime(:) 
!        LOGICAL :: lColumnStorage=.false.
!        LOGICAL :: lColumnStorage=.true.
        PUBLIC :: &
               KSOLVP, NN1, NN2, NN3, KSOLVU, U_ML, StartUTime !, &
               !lColumnStorage
      END MODULE SOLVI 
!                                                                       
      MODULE ITSOLI 
        USE SutraMSPrecision 
        INTEGER (I4B) ::ITRMXP, ITOLP, NSAVEP 
        INTEGER (I4B), ALLOCATABLE::ITRMXU (:), ITOLU (:), NSAVEU (:) 
        PUBLIC ITRMXP, ITOLP, NSAVEP, ITRMXU, ITOLU, NSAVEU 
      END MODULE ITSOLI 
!                                                                       
      MODULE ITSOLR 
        USE SutraMSPrecision 
        REAL (DP) ::TOLP 
        REAL (DP), ALLOCATABLE::TOLU (:) 
        PUBLIC TOLP, TOLU 
      END MODULE  
!                                                                       
      MODULE JCOLS
        USE SutraMSPrecision
        INTEGER (I4B), PARAMETER :: K5MAX=10  !8 CHANGED TO 10 - VERSION 1.1 
        INTEGER (I4B), PARAMETER :: K6MAX=7 
        CHARACTER (LEN = 1) ::K5SYM (K5MAX) 
        CHARACTER (LEN = 1), ALLOCATABLE::K5COL (:) 
        CHARACTER (LEN = 25) ::VARNK5 (K5MAX) 
        CHARACTER (LEN = 15), ALLOCATABLE::CSMSK5 (:) 
        CHARACTER (LEN = 2) ::K6SYM (K6MAX) 
        CHARACTER (LEN = 2), ALLOCATABLE::K6COL (:) 
        CHARACTER (LEN = 25) ::VARNK6 (K6MAX) 
        CHARACTER (LEN = 15), ALLOCATABLE::COLTK5 (:) 
        INTEGER ::NCOLMX 
        INTEGER ::NCOLPR, LCOLPR 
        INTEGER ::NCOLS5, NCOLS6 
        INTEGER, ALLOCATABLE::J5COL (:), J6COL (:) 
        !..VERSION 1.1 - ADDED H AND V
        DATA K5SYM / 'N', 'X', 'Y', 'Z', 'P', 'U', 'S', 'R', 'H', 'V'/
        DATA VARNK5 / 'NODE NUMBER', 'X-COORDINATE', 'Y-COORDINATE', 'Z-COORDINATE', &
                      'PRESSURE/HEAD', 'CONCENTRATION/TEMPERATURE', &
                      'SATURATION', 'DENSITY', 'TOTAL HEAD', 'NODAL VELOCITY' /                                                             
        DATA K6SYM / 'E', 'X', 'Y', 'Z', 'VX', 'VY', 'VZ' / 
        DATA VARNK6 / 'ELEMENT NUMBER', 'X-COORDINATE OF CENTROID', 'Y-COORDINATE OF CENTROID', &
                      'Z-COORDINATE OF CENTROID', 'X-VELOCITY', 'Y-VELOCITY', 'Z-VELOCITY' /                                        
        PUBLIC :: &
                  NCOLMX, NCOLPR, LCOLPR, &
                  NCOLS5, NCOLS6, J5COL, J6COL, &
                  K5SYM,  K5COL,  K5MAX, VARNK5, CSMSK5, &
                  K6SYM,  K6COL,  K6MAX, VARNK6, COLTK5               
      END MODULE JCOLS 
!                                                                       
      MODULE MODSOR 
        USE SutraMSPrecision 
        CHARACTER (LEN = 10), ALLOCATABLE::ADSMOD (:) 
        PUBLIC ADSMOD 
      END MODULE MODSOR 

      MODULE DIMS
        USE SutraMSPrecision
        IMPLICIT NONE
        INTEGER (I4B) :: &
          NN, &
          NE, &
          NIN, &
          NBI, &
          NCBI, &
          NB, &
          NBHALF, &
          NPBC, &
          MNUBC, &
          NSOP, &
          MNSOU, &
          NBCN, &
          NSPE                                           
      END MODULE DIMS

      MODULE DIMX
        USE SutraMSPrecision
        IMPLICIT NONE
        INTEGER (I4B) :: &
          NBIX, &
          NWI, &
          NWF, &
          NWL, &
          NIPARM, &
          NFPARM, &
          NELT, &
          NNNX, &
          NEX, &
          N48
      END MODULE DIMX

      MODULE TIMES
        use MSErrorHandler
        USE SutraMSPrecision
        IMPLICIT NONE
        integer (I4B), parameter :: cTimerTextLen=128
        REAL (DP) :: &
          DELT, &
          TSEC, &
          TMIN, &
          THOUR, &
          TDAY, &
          TWEEK, &
          TMONTH, &
          TYEAR, &
          TMAX, &
          DELTP, &
          DELTU, &
          DLTPM1, &
          DLTUM1, &
          TSTART      
        INTEGER (I4B) :: &
          IT, &
          ITMAX
        logical :: lTimerInitialized=.false.
        real (DP) :: startTimer, endTimer
        character (cTimerTextLen) :: cTimerText
        public :: &
          DELT, TSEC, TMIN, THOUR, TDAY, TWEEK, TMONTH, TYEAR, TMAX, &
          DELTP, DELTU, DLTPM1, DLTUM1, &
          TSTART, &
          cTime, GenericTimer
        contains
          character (len=25) function cTime(tInput)
          USE SutraMSPrecision
          implicit none
          real (DP), optional :: tInput
          !locals
          character (len=15)  :: cNumber
          character (len=8)   :: cTimeUnit
          real (DP)           :: tValue
          real (DP)           :: ontime
          !first line of code
            tValue=TSEC
            if (present(tInput)) tValue=tInput
            ontime=tValue
            cTimeUnit='SEC'
            if(ontime>=60.and.ontime<(60*60)) then
              ontime=ontime/60
              cTimeUnit='MIN'
            else if(ontime>=(60*60).and.ontime<(60*60*24)) then
              ontime=ontime/(60*60)
              cTimeUnit='HRS'
            else if(ontime>=(60*60*24).and.ontime<(60*60*24*365.25)) then
              ontime=ontime/(60*60*24)
              cTimeUnit='DAYS'
            else if(ontime>=(60*60*24*365.25).and.ontime<(60*60*24*365.25*1000)) then
              ontime=ontime/(60*60*24*365.25)
              cTimeUnit='YEARS'
            else if(ontime>=(60*60*24*365.25*1000).and.ontime<(60*60*24*365.25*1000000)) then
              ontime=ontime/(60*60*24*365.25*1000)
              cTimeUnit='K YEARS'
            else if(ontime>=(60*60*24*365.25*1000000)) then
              ontime=ontime/(60*60*24*365.25*1000000)
              cTimeUnit='Ma YEARS'
            end if
            WRITE (cNumber, '(f15.4)') ontime
            cTime=trim(adjustl(cNumber))//' '//trim(adjustl(cTimeUnit))
            return
          end function cTime

          !Generic timer routine
          logical function GenericTimer(cTimerText)
            character (len=cTimerTextLen), intent (in), optional :: cTimerText
            character (len=cTimerTextLen) :: cTempTimerText
            lOk=.false.
            cTempTimerText='Subroutine Time'
            if (present(cTimerText)) cTempTimerText=cTimerText
            if(lTimerInitialized) then
              call CPU_TIME(endTimer)
              write ( * , '(1x,6x,a,1x,a,1x,a)') &
                           trim(adjustl(cTempTimerText))//' (',trim(adjustl(cTime(endTimer - startTimer))),')'
              lTimerInitialized=.false.
            else
              call CPU_TIME(startTimer)
              lTimerInitialized=.true.
            end if
            lOk=.true.
            GenericTimer=lOk
            return
          end function GenericTimer

      END MODULE TIMES

      MODULE KPRINT
        USE SutraMSPrecision
        IMPLICIT NONE
        INTEGER (I4B) :: &
          KNODAL, &
          KELMNT, &
          KINCID, &
          KPLOTP, &
          KPLOTU, &
          KVEL, &
          KBUDG
      END MODULE KPRINT

      MODULE OBS
        USE SutraMSPrecision
        IMPLICIT NONE
        LOGICAL :: &
          lTIMOBSRESET=.FALSE.
        CHARACTER (LEN=128), ALLOCATABLE :: &
          COBSNAME(:)
        INTEGER (I4B) :: &
          NOBSN, &
          NTOBSN, &
          NOBCYC, &
          ITCNT, &
          ONOBS=1
        REAL (DP), ALLOCATABLE :: &
           TIMOBS(:)
        PUBLIC :: &
          lTIMOBSRESET, &
          COBSNAME, &
          NOBSN, NTOBSN, NOBCYC, ITCNT, ONOBS, &
          TIMOBS, &
          AllocateStandardSutraObservations
        contains
          logical function AllocateStandardSutraObservations()
            use DIMS
            use CONTRL
            use TIMES
            use FUNITS
            use MSErrorHandler
            use TotalStorage
            implicit none
            !locals
            integer (I4B) :: &
              JT=0, JTMAX, KT=0, &
              NTOBS, &
              n
            real (DP) :: &
              DELTK, TS
            real (DP), allocatable :: TEMP(:)
            
            ios=.false.
            allocate(Temp(ITMAX+1),stat=ios) !dimension Temp to max number of timesteps + INITIAL
            if (ios /= 0) goto 9999
            TS=TSTART
            DELTK=DELT
!         10 CONTINUE 
            do while (JT.LT.ITMAX.AND.TS.LT.TMAX)
              JT = JT + 1 
              IF(MOD(JT, ITCYC).EQ.0.AND.JT.GT.1) DELTK = DELTK * DTMULT 
              IF(DELTK.GT.DTMAX) DELTK = DTMAX 
              TS = TS + DELTK 
              IF(MOD(JT, NOBCYC).EQ.0.OR.JT.EQ.1)THEN
                KT = KT + 1 
                Temp(KT)=TS
              END IF
              !IF(JT.LT.ITMAX.AND.TS.LT.TMAX) GOTO 10
            end do
            JTMAX = JT 
            IF(JTMAX.GT.1.AND.MOD(JT, NOBCYC).NE.0)THEN
              KT = KT + 1 
              Temp(KT)=TMAX
            END IF
            NTOBS = KT
            NTOBSN = NTOBS + 2 
            allocate( &
                      TIMOBS(NTOBS), &
                      stat=ios)
            if (ios /= 0) goto 9999
            if(.not. allocated(COBSNAME)) &
            allocate( &
                      COBSNAME(NOBSN), &
                      stat=ios)
            if (ios /= 0) goto 9999
            !calculate memory requirements
            ios=AddMemory(MemoryIndex('OBS'),'Vec',                 NTOBS) !TIMOBS(NTOBS)
            ios=AddMemory(MemoryIndex('OBS'),'Chr',len(COBSNAME(1))*NOBSN) !COBSNAME(NOBSN)

            !set TIMEOBS to Temp
            TIMOBS=Temp(1:NTOBS)
            !initialize COBSNAME
            if (fSOB==0) then
              do n=1,NOBSN
                COBSNAME(n)='Observation '//trim(adjustl(Val2Char(n)))
              end do
            end if

            !write Observation times to fLST
            write (fLST,'(//a/100(10(:f15.2,1x)/))') &
              'Simulation data will be output at '//trim(adjustl(Val2Char(NOBSN)))//' locations at the following times', &
              (TIMOBS(n),n=1,NTOBS)

            ios=.true.
09999&
            AllocateStandardSutraObservations=lOk
            if(allocated(Temp)) deallocate(Temp)
            return
          end function AllocateStandardSutraObservations

      END MODULE OBS 

      MODULE GRAVEC
        USE SutraMSPrecision
        IMPLICIT NONE
        REAL (DP) :: &
          GRAVX, &
          GRAVY, &
          GRAVZ 
      END MODULE GRAVEC

      MODULE SOLVC
        IMPLICIT NONE
        CHARACTER (LEN=10), DIMENSION(0:10) :: &
          SOLWRD = &
           (/ 'DIRECT  ', &
              'CG      ', &
              'GMRES   ', &
              'ORTHOMIN', &
              'UNKNOWN ', &
              '        ', &
              '        ', &
              '        ', &
              '        ', &
              '        ', &
              '        ' /)
        CHARACTER (LEN=40), DIMENSION(0:10) :: &
            SOLNAM = &
             (/ 'BANDED GAUSSIAN ELIMINATION (DIRECT)', &
                'IC-PRECONDITIONED CONJUGATE GRADIENT', &
                'ILU-PRECONDITIONED GMRES            ', &
                'ILU-PRECONDITIONED ORTHOMIN         ', &
                'Hook for additional solver          ', &
                '                                    ', &
                '                                    ', &
                '                                    ', &
                '                                    ', &
                '                                    ', &
                '                                    ' /)
      END MODULE SOLVC

      MODULE SOLVN
        USE SutraMSPrecision
        IMPLICIT NONE
        INTEGER (I4B) :: &
          NSLVRS = 5
      END MODULE SOLVN

      MODULE PLT1
        USE SutraMSPrecision
        IMPLICIT NONE
        LOGICAL :: &
          ONCEK5, &
          ONCEK6, &
          ONCEK7 
      END MODULE PLT1

      MODULE SOURCEITEMS
        USE SutraMSPrecision
        IMPLICIT NONE
        integer (I4B) :: &
          IQSOPT, &
          IQSOUT, &
          IPBCT, &
          IUBCT, &
          IBCT
      END MODULE SOURCEITEMS


                                                                      
