
!     SUBROUTINE     S  U  T  R  A  E  N  D    SUTRA-MS VERSION 2004.1
! ***
! *** PURPOSE :                                                         
! ***  SUBROUTINE THAT IS USED TO TERMINATE SUTRA-MS.  THIS SUBROUTINE
! ***  REPLACES ALL STOP COMMANDS THAT WERE PRESENT IN SUTRA.
! ***


    SUBROUTINE SutraEnd()
      
      USE FUNITS
      USE MSErrorHandler
      
      IMPLICIT NONE

        if(LNormal) then
          write(*,'(1x,4x,a)') 'Normal Model Termination'
          if(fLST>0) write(fLST,'(1x,4x,a)') 'Normal Model Termination'
        else
          write(*,'(1x,4x,a)') 'Abnormal Model Termination'
          write(*,'(1x,4x,a)') 'See listing file for error details'
          if(fLST>0) then
            write(fLST,'(1x,4x,a)') 'Abnormal Model Termination'
          end if
        end if

        !close all open files
        CALL CLOSEFILE()

        !only STOP termination in program
        STOP

    end subroutine SutraEnd