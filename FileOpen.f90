!     SUBROUTINE        F  O  P  E  N          SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  OPENS FILES FOR SUTRA SIMULATION.                                
! ***  OPENS ERROR OUTPUT FILE, READS FILE NUMBERS AND NAMES,           
! ***  CHECKS FOR EXISTENCE OF INPUT FILES, AND WRITES ERROR MESSAGES.  
!                                                                       
      SUBROUTINE FOPEN () 
      USE FUNITS 
      USE SOLVI
      !MS specific modules
      USE SutraMSPrecision 
      USE MSErrorHandler
      USE ColumnStorage
      USE TECPLOT, ONLY: TecplotAll
      !
      IMPLICIT NONE
      !LOCAL VARIABLES
        CHARACTER (LEN = 3) :: CK
        CHARACTER (LEN=80)  :: FN
        CHARACTER (LEN=20) :: cStatus
        LOGICAL :: LOLD
        LOGICAL :: lOpened
        INTEGER :: &
          NLSKIP, &
          IU, &
          KERR, &
          NFILE, &
          NF
        integer (kind=4)   :: narg
        integer (kind=4)   :: nlen
        character (len=80) :: ctmp
        integer (kind=4)   :: i
!.....INITIALIZE fLST to fSMY in case fLST is not specified in SUTRA.FIL
        fLST = fSMY
!
!.....COMMENT NEXT TWO LINES OUT IF USING COMPAQ COMPILER WITH DFLIB
      lColumnStorage=.true.
      write (*,*) 'Using AP Column Storage Routine'
!
!.....GET COMMAND LINE ARGUMENTS FOR COMPAQ COMPILER
!     Comment out for all other compilers
      narg = iargc() - 1
      i=0
      do i = 1, narg
        call getarg(i,ctmp)
        select case(trim(adjustl(ctmp)))
          case ("-d", "-D")  !command line argument for output of special debug information to fLST
            lDebugData=.true.
            write (*,*) 'Special Debug Data Output to *.LST File'
          case ("-f", "-F")  !command line argument for output of special debug information to fLST
            cSutraFil=trim(adjustl(ctmp(3:LEN_TRIM(ctmp))))
            write (*,*) 'Processing simulation "SUTRA.FIL"...',trim(adjustl(cSutraFil))
          case ("-tri","-TRI")
            lColumnStorage = .false.
            write (*,*) 'Using SLAP TRIAD STORAGE FORMAT'
          case ("-tecplotall","-TECPLOTALL")
            TecplotAll = .true.
            write (*,*) 'Full Tecplot dataset will be created'
          case default
            write (*,*) 'Unknown command line option ['//ctmp(1:2)//']'
        end select
      end do
!     Comment out for all other compilers
!                                                                       
!.....OPEN FILE UNIT CONTAINING UNIT NUMBERS AND FILE ASSIGNMENTS       
      MSErrorValue%cDataSet='  0'
      IU = fSutraFil 
      FN = cSutraFil
      cStatus=cOld
      if(.not.OpenFormattedFile(IU,FN,trim(adjustl(cStatus)))) goto 9999
      if(ios/=0) call ErrorIO('Error opening '//trim(adjustl(FN)))
!                                                                       
!.....READ FILE CONTAINING UNIT NUMBERS AND FILE ASSIGNMENTS            
      NFILE = 0 
  100 READ (fSutraFil, *, END = 200) CK, IU, FN
      cStatus=cUnknown
      IU = FindOpenUnit()
      select case (CK) 
         case ('SMY')
            fSMY = IU
            cSMY = FN
         case ('K1','INP')
            fINP = IU
            cINP = FN
            cStatus=cOld
         case ('K2','ICS')
            fICS = IU
            cICS = FN
            cStatus=cOld
         case ('K3','LST')
            fLST = IU
            cLST = FN
         case ('K4','RST')
            fRST = IU
            cRST = FN
         case ('K5','NOD')
            fNOD = IU
            cNOD = FN
         case ('K6','ELE')
            fELE = IU
            cELE = FN
         case ('K7','OBS')
            fOBS = IU
            cOBS = FN
         case ('TBC')
            fTBC = IU
            cTBC = FN
            cStatus=cOld
         case ('OTM')
            fOTM = IU
            cOTM = FN
            cStatus=cOld
!.......AUTOMATIC TIMESTEPPING                                          
         case ('ATS')
            fATS = IU
            cATS = FN
            cStatus=cOld
!.......SPECIFIED OBSERVATION LOCATIONS - optional
         case ('SOB')
            fSOB = IU
            cSOB = FN
            cStatus=cOld
!.......ZONE FILE INSTEAD OF DATASET 14B AND 15B
         case ('ZON')
            fZON = IU
            cZON = FN
            cStatus=cOld
!.......TECPLOT NODAL OUTPUT - VERSION 1.1
         case ('TPN')
           fTPN = IU
           cTPN = FN
!.......TECPLOT ELEMENTAL OUTPUT - VERSION 1.1
         case ('TPE')
           fTPE = IU
           cTPE = FN
!.......BUDGET OUTPUT 
         case ('BUD')
           fBUD = IU
           cBUD = FN
         case default
            goto 100
      end select
      !open file
      if(.not.OpenFormattedFile(IU,FN,trim(adjustl(cStatus)))) goto 9999
      goto 100
  200 continue 

!.....OPEN FILE UNIT FOR ERROR MESSAGES -
!       if not defined in SUTRA.fil and default file is being used
      if (index(cSMY,'SUTRA.SMY')>0) then
        inquire (unit=fSMY, OPENED=lOpened)
        if (.not. lOpened) then
          cStatus=cUnknown
          if(.not.OpenFormattedFile(fSMY,cSMY,trim(adjustl(cStatus)))) goto 9999
        end if
      end if
!                                                                       
!.....CLOSE SUTRA.FIL                                                   
09999&
      CLOSE (fSutraFil) 
!.....RETURN TO CALLING ROUTINE                                         
      RETURN 
!                                                                       
      END SUBROUTINE FOPEN                          
