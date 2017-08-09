!     MODULE            MSUTRAMSERRORHANDLER   SUTRA-MS VERSION 2004.1
!                                                                       
! ***  
! *** PURPOSE :                                                         
! ***  MODULE CONTAINING ALL OF THE ERROR HANDLING FUNCTIONALITY
! ***  USED IN SUTRA-MS.
! ***
! ***  MODULES ALSO CONTAIN THE FOLLOWING SUBROUTINES:
! ***    VAL2CHAR - 
! ***      OVERLOADED FUNCTION THAT CONVERTS
! ***      INT, SP, AND DP VARIABLES TO A TEXT
! ***      STRING
! ***

      module MSErrorHandler
        USE SutraMSPrecision

        type MSErrorType
          character (len=3) :: cDataSet
        end type

        type(MSErrorType) :: MSErrorValue

        !overloaded function
        interface Val2Char
          module procedure Int2Char,Real2Char,Dbl2Char
        end interface Val2Char

        logical             :: LOk
        logical             :: LNormal
        integer (I4B)       :: ios
        public MSErrorValue, LOk, LNormal, ios, Val2Char

        contains
          character(len=16) function Int2Char(ival)
            implicit none
            integer (I4B) :: ival
            write(Int2Char,'(i16)') ival
            return
          end function Int2Char

          character(len=16) function Real2Char(rval)
            implicit none
            real (SP) :: rval
            write(Real2Char,'(g15.7)') rval
            return
          end function Real2Char

          character(len=16) function Dbl2Char(dval)
            implicit none
            real (DP) :: dval
            write(Dbl2Char,'(g15.7)') dval
            return
          end function Dbl2Char


      end module MSErrorHandler


