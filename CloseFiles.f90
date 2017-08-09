!     SUBROUTINE   C  L  O  S  E  F  I  L  E   SUTRA-MS VERSION 2004.1

      SUBROUTINE CLOSEFILE 
      USE FUNITS 
      use SutraMSPrecision
      implicit none
      INTEGER (I4B) :: i
!.....fSMY
      if(fSMY/=fLST) &
        i = CloseFormattedFile(fSMY)
!.....fSutraFil                                                                
      i = CloseFormattedFile(fSutraFil)
!.....fINP                                                                
      i = CloseFormattedFile(fINP)
!.....fICS                                                                
      i = CloseFormattedFile(fICS)
!.....fRST                                                                
      i = CloseFormattedFile(fRST)
!.....fNOD                                                                
      i = CloseFormattedFile(fNOD)
!.....fELE                                                                
      i = CloseFormattedFile(fELE)
!.....fOBS                                                                
      i = CloseFormattedFile(fOBS)
!.....fTBC                                                                
      i = CloseFormattedFile(fTBC)
!.....SPECIFIED OUTPUT PRINT TIMES
      i = CloseFormattedFile(fOTM)
!.....AUTOMATIC TIME STEP FILE
      i = CloseFormattedFile(fATS)
!.....ZONE FILE
      i = CloseFormattedFile(fZON)
!.....SPECIFIED OBSERVATION LOCATION FILE
      i = CloseFormattedFile(fSOB)
!.....TECPLOT NODAL AND ELEMENTAL FILES - VERSION 1.1
      i = CloseFormattedFile(fTPN)
      i = CloseFormattedFile(fTPE)
!
!.....fLST...Closed last to allow all possible errors to be written
      i = CloseFormattedFile(fLST)
!.....RETURN TO CALLING ROUTINE                                         
      RETURN 
      END SUBROUTINE CLOSEFILE                      
