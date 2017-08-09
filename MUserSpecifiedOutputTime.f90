!     MODULE        USERSPECIFIEDOUTPUTTIME    SUTRA-MS VERSION 2004.1
                                                                       
! ***
! *** PURPOSE :                                                         
! ***  MODULE THAT IS USED TO READ USER SPRCIFIED OUTPUT TIMES
! ***  (READUSERSPECIFIEDOUTPUTTIME).  THE MODULE GIVES GLOBAL ACCESS
! ***  TO THE VARIABLE DPTIME WHICH IS THE NEXT USER SPECIFIED OUTPUT
! ***  TIME.
! ***

      !OTM package
      module UserSpecifiedOutputTime
        use SutraMSPrecision
        !MS specific modules
        !use MSPrecision 
        implicit none
        real (DP) :: DPTIME
        public :: &
          DPTIME, &
          ReadUserSpecifiedOutputTime

        contains
          logical function ReadUserSpecifiedOutputTime()
          USE FUNITS 
          USE TIMES
          USE MSErrorHandler
          implicit none
          !Locals
          character (len=20) :: &
            CTERM
          integer (I4B) :: &
            NLSKIP
!
          lOk=.false.
          !read next output print time                                       
          ! DPTIME only controls output to nodal and elemental output files   
          MSErrorValue%cDataSet='OTM'
          CALL SKPCOM (fOTM, NLSKIP) 
          READ (fOTM, *, IOSTAT = ios) CTERM
          if(ios/=0) call ErrorIO('Error specifying next user specified printing time (DPTIME)')
          IF(TRIM(CTERM)=='END'.OR.TRIM(CTERM)=='end') THEN
            DPTIME=TMAX+1.0D0
            GOTO 9998
          END IF
          READ (CTERM, *, IOSTAT = ios) DPTIME 
          if(ios/=0) call ErrorIO('Internal Error converting CTERM ['//trim(adjustl(CTERM))//'] to DPTIME')
          WRITE (fLST,'(//A,1X,G15.7,1X,A)') '*** NEXT PRINT TIME :',DPTIME,' ***'
09998&
          lOk=.true.
09999&
          ReadUserSpecifiedOutputTime=lOk
          !return to calling routine                                         
          return
          end function ReadUserSpecifiedOutputTime

      end module UserSpecifiedOutputTime
