!     MODULE            A  T  S  D  A  T  A    SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  MODULE CONTAINING VARIABLES, AND INPUT FUNCTION
! ***  FOR AUTOMATIC TIMESTEP (ATS) OPTION.
!                                                                       
      MODULE  ATSDATA
        use SutraMSPrecision
        IMPLICIT NONE
        LOGICAL :: LATS = .false.
        LOGICAL :: LContinue
        LOGICAL :: LATSAllowTimeStepIncrease = .true.
        LOGICAL :: LATSFailConvergence=.false.
        INTEGER (I4B) :: imaxtarget,imaxiter
        INTEGER (I4B) :: &
          NConvMax=100000, &
          NMaxCount
        REAL     (DP) :: dtmin

        public :: &
          LATS, &
          LContinue, &
          LATSAllowTimeStepIncrease, &
          LATSFailConvergence, &
          imaxtarget,imaxiter, &
          NConvMax, NMaxCount, &
          dtmin, &
          MAutomaticTimeStep


        contains
          !Logical function to setup Automatic Time Step Routine
          logical function MAutomaticTimeStep()
            use FUNITS
            use SutraMSPrecision
            !use MSPrecision 
            use MSErrorHandler
            implicit none
            INTEGER (I4B) :: &
              NLSKIP
!
            lOk=.false.
            !ATS DATASET 1
            !READ AUTOMATIC TIMESTEP DATA
            !MAXIMUM ITERATION TARGET VALUE (imaxtarget), 
            !MINIMUM ALLOWED TIMESTEP (dtmin), and
            !whether the timestep should be rerun (LContinue)
            MSErrorValue%cDataSet='ATS'
            CALL SKPCOM(fATS, NLSKIP)
            READ(fATS,*,IOSTAT=ios) imaxtarget,dtmin,LContinue
            if(ios/=0) call ErrorIO('Error specifying DATASET 1 of ATS package')
            write(fLST,'(//1x,"ATS Package DATA"/1x,"imaxtarget ",i10,&
                                              /1x,"dtmin      ",g15.7,&
                                              /1x,"LContinue  ",l10)')           &
                                              imaxtarget,dtmin,LContinue
            !ATS DATASET 2
            !READ ONLY IF LContinue is true
            !MAXIMUM NUMBER OF TIMES A TIMESTEP SHOULD BE RERUN
            !IF imaxtarget IS EXCEEDED
            IF(.not.LContinue) then
              CALL SKPCOM(fATS, NLSKIP)
              READ(fATS,*,IOSTAT=ios) NConvMax
              if(ios/=0) call ErrorIO('Error specifying DATASET 2 of ATS package')
              write(fLST,'(1x,"NConvMax   ",i10//)')           &
                         NConvMax
            END IF
            !SUCCESSFUL READ OF ATD FILE 
            !RETURN TO CALLING ROUTINE
            lOk=.true.
00999 &
            MAutomaticTimeStep=lOk
            LATS = .true.
            RETURN
          end function MAutomaticTimeStep

      END MODULE  ATSDATA
