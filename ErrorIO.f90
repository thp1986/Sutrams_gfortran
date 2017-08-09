!     SUBROUTINE        E  R  R  O  R  I  O    SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  STANDARD INTERFACE FOR OUTPUTTING ERROR INFORMATION TO THE
! ***  SCREEN AND THE LISTING (LST) FILE.
! ***  SUBROUTINE ALSO TERMINATES THE PROGRAM.
!                                                                       
!                                                                       
      subroutine ErrorIO(cError)
        USE FUNITS
        USE MSErrorHandler
        implicit none
        character *(*), intent(in)       :: cError
        character (len=40) :: cErrText
        character (len=16), dimension(5) :: &
          cEType= (/'Unknown         ',        &
                    'End of Record   ',  &
                    'End of File     ',    &
                    'Normal          ',         &
                    'Fortran Run-Time'/)
          !Get error text
          !Call GERROR (cErrText)
          !write error message
          write(*,100) MSErrorValue%cDataSet,cError(1:min(60,len_trim(cError)))
          if(fLST>0) write(fLST,105) MSErrorValue%cDataSet,cError
          select case (ios)
            case(:-3)
              write(*,110) cEType(1),ios,cErrText
              if(fLST>0) write(fLST,110) cEType(1),ios,cErrText
            case(-2)
              write(*,110) cEType(2),ios,cErrText
              if(fLST>0) write(fLST,110) cEType(2),ios,cErrText
            case(-1)
              write(*,110) cEType(3),ios,cErrText
              if(fLST>0) write(fLST,110) cEType(3),ios,cErrText
            !should not execute this case if have entered ErrorIO subroutine
            case(0)
              write(*,110) cEType(4),ios
              if(fLST>0) write(fLST,110) cEType(4),ios
            case(1:58)
              write(*,110) cEType(5),ios,cErrText
              if(fLST>0) write(fLST,110) cEType(5),ios,cErrText
            case(59)
              write(*,110) cEType(5),ios,cErrText
              if(fLST>0) write(fLST,110) cEType(5),ios,cErrText
            case(60:)
              write(*,110) cEType(5),ios,cErrText
              if(fLST>0) write(fLST,110) cEType(5),ios,cErrText
          end select
100&
          format(//, &
                 1x,'*** Input Data Error ***',/,1x,4x,'Data Set ',a3,/, &
                 1x,4x,'Specific Error Information:',/, &
                 1x,4x,a,' . . . ')
105&
          format(//, &
                 1x,'*** Input Data Error ***',/,1x,4x,'Data Set ',a3,/, &
                 1x,4x,'Specific Error Information:',/, &
                 1x,4x,a)
110&
          format(1x,4x,'Error is of type',1x,a16,1x,'[',i3,']',/&
                 1x,4x,a)
          !terminate program
          call SutraEnd()
        return
      end subroutine ErrorIO