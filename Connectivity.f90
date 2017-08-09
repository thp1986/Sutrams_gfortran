!     SUBROUTINE        C  O  N  N  E  C       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO READ ,ORGANIZE, AND CHECK DATA ON NODE INCIDENCES.            
!                                                                       
      SUBROUTINE CONNEC (IN) 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE KPRINT
      !MS specific modules
      use SutraMSPrecision
      USE MSErrorHandler
      implicit none
      CHARACTER (LEN=10) :: CDUM10 * 10 
      integer (I4B) :: &
        IN(NIN), &
        IIN(8)
      !Locals
        integer (I4B) :: &
          IPIN=0, &
          NLSKIP
        integer (I4B) :: &
          II, &
          III, &
          L, &
          LL, &
          M, &
          M1, &
          M8
!                                                                       
      IF (KINCID.EQ.0) WRITE (fLST, 1) 
    1 FORMAT(1H1////11X,'M E S H   C O N N E C T I O N   D A T A'//     &
     &   16X,'PRINTOUT OF NODAL INCIDENCES CANCELLED.')                 
      IF (KINCID.EQ. + 1) WRITE (fLST, 2) 
    2 FORMAT(1H1////11X,'M E S H   C O N N E C T I O N   D A T A',      &
     &   ///11X,'**** NODAL INCIDENCES ****'///)                        
!                                                                       
!.....INPUT DATASET 22 AND CHECK FOR ERRORS                             
      MSErrorValue%cDataSet=' 22'
      CALL SKPCOM (fINP, NLSKIP) 
      READ (fINP,*,iostat=ios) CDUM10 
      if(ios/=0) call ErrorIO('Error specifying line 1')
      IF (trim(CDUM10).NE.'INCIDENCE') THEN 
         call ErrorIO('DATASET 22 must begin with the word "INCIDENCE"')
      ENDIF 
      DO 1000 L = 1, NE 
         READ (fINP,*,iostat=ios) LL
         if(ios/=0) call ErrorIO('Error specifying element number - error on line ['//trim(Val2Char(L))//'] of Data Set 22')
         if(LL>NE)  call ErrorIO('Error: Specified element number ['//trim(Val2Char(LL))//'] exceeds NE ['//trim(Val2Char(NE))//']')
         backspace(fINP)
         READ (fINP,*,iostat=ios) LL, (IIN (II), II = 1, N48) 
         if(ios/=0) call ErrorIO('Error specifying incidence data - error on line ['//trim(Val2Char(L))//'] of Data Set 22')
!.....PREPARE NODE INCIDENCE LIST FOR MESH, IN.                         
         DO 5 II = 1, N48 
            III = II + (L - 1) * N48 
    5    IN (III) = IIN (II) 
         IF (IABS (LL) .EQ.L) GOTO 500 
         call ErrorIO('Incidence data for element ['//trim(Val2Char(LL))//'] is not in numerical order in the data set')
!                                                                       
!                                                                       
  500    M1 = (L - 1) * N48 + 1 
         M8 = M1 + N48 - 1 
         IF (KINCID.EQ.0) GOTO 1000 
         WRITE (fLST, 650) L, (IN (M), M = M1, M8) 
  650 FORMAT(11X,'ELEMENT',I9,5X,' NODES AT : ',6X,'CORNERS ',          &
     &   5(1H*),8I9,1X,5(1H*))                                          
!                                                                       
 1000 END DO 
!                                                                       
!.....return to calling routine                                                                       
 5000 RETURN 
      END SUBROUTINE CONNEC                         
