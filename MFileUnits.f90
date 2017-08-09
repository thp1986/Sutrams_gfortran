! ***  
! *** PURPOSE :                                                         
! ***  MODULE CONTAINING COMMON PROGRAM VARIABLES AND USED
! ***  TO REPLACE FUNITS COMMON STATEMENTS USED IN STANDARD
! ***  VERSION OF SUTRA.
! ***  
! ***  MODULE ALSO CONTAIN THE FOLLOWING SUBROUTINES:
! ***    FILEOPENUNIT - 
! ***      FUNCTION FOR DETERMINING THE NEXT
! ***      AVAILABLE FORTRAN UNIT NUMBER. 
! ***      ELIMINATES THE NEED TO SPECIFY UNIQUE
! ***      UNIT NUMBERS FOR FILES IN SUTRA.FIL
! ***    OPENFORMATTEDFILE -
! ***      FUNCTION THAT HANDLES THE OPENING OF A
! ***      UNFORMATTED FILE. USED FOR ALL FILE OPENING.
! ***    CLOSEFORMATTEDFILE - 
! ***      FUNCTION THAT HANDLES THE CLOSING OF A 
! ***      UNFORMATTED FILE.  USED FOR ALL FILE CLOSING.
! ***  

      MODULE  FUNITS
        USE SutraMSPrecision
        CHARACTER (LEN=20), parameter     :: cOld='OLD'
        CHARACTER (LEN=20), parameter     :: cUnknown='UNKNOWN'
        logical                           :: l_onunit=.true.
        LOGICAL                           :: lDebugData=.false.
        CHARACTER (LEN=80)                :: &
                                             cSMY, cSutraFil, &
                                             cINP, cICS, cLST, cRST, &
                                             cNOD, cELE, cOBS, &
                                             cTBC, cGTP, &
                                             cOTM, cATS, cZON, cSOB, &
                                             cTPN, cTPE  !VERSION 1.1
        INTEGER (I4B)                     :: fSMY, fSutraFil
        INTEGER (I4B)                     :: fINP, fICS, fLST, fRST
        INTEGER (I4B)                     :: fNOD, fELE, fOBS
        INTEGER (I4B)                     :: fTBC, fGTP
        INTEGER (I4B)                     :: fOTM, fATS, fZON, fSOB
        INTEGER (I4B)                     :: fTPN, fTPE  !VERSION 1.1
        INTEGER (I4B), PARAMETER          :: NKFLE=20
        INTEGER (I4B)                     :: n_onunit=19
        PUBLIC  :: &
                cOld, cUnknown, &
                lDebugData, &
                cSMY, cSutraFil, &
                cINP, cICS, cLST, cRST, &
                cNOD, cELE, cOBS, &
                cTBC, cGTP, &
                cOTM, cATS, cZON, cSOB, &
                fSMY,fSutraFil,&
                fINP,fICS,fLST,fRST,fNOD,fELE,fOBS,&
                fTBC,fGTP,fOTM,fATS,fZON,fSOB, &
                NKFLE, &
                n_onunit, &
                FindOpenUnit, &
                OpenFormattedFile, &
                CloseFormattedFile
        contains
          integer function FindOpenUnit()
            logical :: lOpened
	    !find open file unit
            l_onunit=.true.
            do while (l_onunit)
              n_onunit=n_onunit+1
              inquire(unit=n_onunit,OPENED=lOpened)
              if(.not.lOpened) l_onunit=.false.
            end do
            FindOpenUnit=n_onunit
            return
          end function FindOpenUnit
          logical function OpenFormattedFile(IU,FN,cStatus)
            character *(*), intent(in)      :: cStatus
            character (len=80), intent(in)  :: FN
            integer (I4B), intent(in)       :: IU
            OpenFormattedFile=.false.
            open (UNIT = IU, FILE = FN, STATUS = cStatus, FORM = 'FORMATTED', IOSTAT = ios)                              
            if(ios/=0) call ErrorIO('Error opening '//trim(adjustl(FN)))
            OpenFormattedFile=.true.
09999&
            return
          end function OpenFormattedFile
          integer function CloseFormattedFile(IU)
            integer (I4B), intent(in) :: IU
            logical :: lOpened
            CloseFormattedFile=0
            inquire (unit = IU, opened = lOpened) 
            if (.not.lOpened) goto 9999
            close (IU) 
            CloseFormattedFile=1
09999&
            return
          end function CloseFormattedFile

      END MODULE  FUNITS

