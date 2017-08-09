subroutine ReadSorption14C()
  use CONTRL
  use DIMS
  use FUNITS
  use MODSOR 
  use PARAMS 
  use SutraMSPrecision
  !MS specific modules
  use MSErrorHandler
  use SutraStorage
  use SutraZoneModule
  implicit none
  !local variables
  integer (I4B) :: NLSKIP
  integer (I4B) :: I, II
  integer (I4B) :: K
  !functions
  !input formats
  !output formats
2000 FORMAT(////11X,'N O D A L  A D S O R P T I O N   P A R A M E T E R S' &
     &        /11X,'      NODE',1X,'   SPECIES',                          &
     &          1X,'           CHI1',1X,'           CHI2')
2010 FORMAT( 11X,2(I10,1X),2(1PD15.4,1X))
  !code
!..FOR ENERGY TRANSPORT:
  if ( ME == 1 ) goto 9999
!..DETERMINE IF ADSORPTION ACTIVE FOR AT LEAST ONE SPECIES
  I = 0
  do K = 1, NSPE
    if (ADSMOD(K) /= 'NONE      ') then
      I = 1
      exit
    end if
  end do
  if ( I == 0 ) goto 9999
!...read nodal data for nodes that differ from the default values specified in Data Set 11
  MSErrorValue%cDataSet='14C'
  CALL SKPCOM (fINP, NLSKIP) 
  I = 1
  do
     READ(fINP,*,iostat=ios) II
     if (ios /= 0) call ErrorIO('Error specifying node number - error on line ['//trim(Val2Char(I))//'] of Data Set 14C')
     if (II > NN)  call ErrorIO('Error: Specified node number ['//trim(Val2Char(II))//'] exceeds NN ['//trim(Val2Char(NN))//']')
     if (II < 1)   exit
     backspace(fINP)
     READ(fINP,*,iostat=ios) II, (ProdSorp(II)%chi1(K),ProdSorp(II)%chi2(K),K=1,NSPE)
     if (ios /= 0) call ErrorIO('Error specifying ProdSorp sorption data - error on line ['//trim(Val2Char(I))//'] of Data Set 14C')
     !adjust if neccesary
      select case (ME)
!.........FOR SOLUTE ONLY TRANSPORT:                                        
        case(-1)
          !no adjustment
        case(0)
!...........FOR SOLUTE AND ENERGY TRANSPORT
          if ( K == NESP ) then                                   
            ProdSorp(II)%chi1(K) = 0.0D0 
            ProdSorp(II)%chi2(K) = 0.0D0 
          end if 
!.........FOR ENERGY TRANSPORT:
        case(1)
          ProdSorp(II)%chi1(K) = 0.0D0 
          ProdSorp(II)%chi2(K) = 0.0D0 
        case default
          call ErrorIO('Error:  ME must be -1, 0, or 1 ['//trim(Val2Char(ME))//']')
      end select
     !write data to fLST
     !write header
     if ( I == 1 ) then
       WRITE(fLST,2000)
     end if        
     !write data
     do K = 1, NSPE
       WRITE(fLST,2010) II,K,ProdSorp(II)%chi1(K),ProdSorp(II)%chi2(K)
     end do
     I = I + 1
  end do
!...return to calling routine
09999 &
  return
end subroutine ReadSorption14C

subroutine ReadSorption14D()
  use CONTRL
  use DIMS
  use FUNITS 
  use PARAMS 
  use SutraMSPrecision
  !MS specific modules
  use MSErrorHandler
  use SutraStorage
  use SutraZoneModule
  implicit none
  !local variables
  integer (I4B) :: NLSKIP
  integer (I4B) :: I, II
  integer (I4B) :: K
  !functions
  !input formats
  !output formats
2000 FORMAT(////11X,'N O D A L  P R O D U C T I O N   P A R A M E T E R S' &
     &        /11X,'      NODE',1X,'   SPECIES',                          &
     &          1X,'         PRODF0',1X,'         PRODS0', &
     &          1X,'         PRODF1',1X,'         PRODS1')
2010 FORMAT( 11X,2(I10,1X),4(1PD15.4,1X))
  !code
!..DETERMINE IF PRODUCTION ACTIVE FOR AT LEAST ONE SPECIES
  I = 0
  do K = 1, NSPE
    if ( PRODF0(K)/=0.0D0 ) then
      I = 1
      exit
    end if
    if ( PRODS0(K)/=0.0D0 ) then
      I = 1
      exit
    end if
    if ( PRODF1(K)/=0.0D0 ) then
      I = 1
      exit
    end if
    if ( PRODS1(K)/=0.0D0 ) then
      I = 1
      exit
    end if
  end do
  if ( I == 0 ) goto 9999
!...read nodal data for nodes that differ from the default values specified in Data Set 12
  MSErrorValue%cDataSet='14D'
  CALL SKPCOM (fINP, NLSKIP) 
  I = 1
  do
     READ(fINP,*,iostat=ios) II
     if (ios /= 0) call ErrorIO('Error specifying node number - error on line ['//trim(Val2Char(I))//'] of Data Set 14D')
     if (II > NN)  call ErrorIO('Error: Specified node number ['//trim(Val2Char(II))//'] exceeds NN ['//trim(Val2Char(NN))//']')
     if (II < 1)   exit
     backspace(fINP)
     READ(fINP,*,iostat=ios) II, (ProdSorp(ii)%prodf0(K),ProdSorp(ii)%prods0(K), &
                                  ProdSorp(ii)%prodf1(K),ProdSorp(ii)%prods1(K),K=1,NSPE)
     if (ios /= 0) call ErrorIO('Error specifying ProdSorp production data - error on line &
['//trim(Val2Char(I))//'] of Data Set 14D')
     !adjust if neccesary
      select case (ME)
!.........FOR SOLUTE ONLY TRANSPORT:                                        
        case(-1)
          !no adjustment
        case(0)
!...........FOR SOLUTE AND ENERGY TRANSPORT
          if ( K == NESP ) then                                   
            ProdSorp(II)%prodf1(K) = 0.0D0 
            ProdSorp(II)%prods1(K) = 0.0D0 
          end if 
!.........FOR ENERGY TRANSPORT:
        case(1)
          ProdSorp(II)%prodf1(K) = 0.0D0 
          ProdSorp(II)%prods1(K) = 0.0D0 
        case default
          call ErrorIO('Error:  ME must be -1, 0, or 1 ['//trim(Val2Char(ME))//']')
      end select
     !write data to fLST
     !write header
     if ( I == 1 ) then
       WRITE(fLST,2000)
     end if        
     !write data
     do K = 1, NSPE
       WRITE(fLST,2010) II,K,ProdSorp(II)%prodf0(K),ProdSorp(II)%prods0(K), &
                             ProdSorp(II)%prodf1(K),ProdSorp(II)%prods1(K)
     end do
     I = I + 1
  end do
!...return to calling routine
09999 &
  return
end subroutine ReadSorption14D
