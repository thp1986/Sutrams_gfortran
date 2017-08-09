!     MODULE            TRIAD2COLUMN           SUTRA-MS VERSION 2004.1
                                                                       
! ***
! *** PURPOSE :                                                         
! ***  MODULE THAT IS CONVERT THE TRIAD SOLVER STORAGE FORMAT TO THE
! ***  COLUMN SOLVER STORAGE FORMAT.  THIS MODULE WAS CREATED BECAUSE
! ***  THE STANDARD VERSION OF SUTRA USED THE TRIAD SOLVER STORAGE 
! ***  METHOD WHICH IS CONVERTED TO THE COLUMN STORAGE FORMAT INTERNALLY
! ***  WITHIN THE SLAP SOLVER SUBROUTINES BUT DOES NOT RESTORE THE TRIAD
! ***  FORMAT AFTER CONVERGENCE.  THIS RESULTED IN THE NEED TO RECALCULATE
! ***  THE TRIAD STORAGE IA AND JA POINTERS AND SUBSEQUENT CONVERSION
! ***  TO THE COLUMN STORAGE FORMAT WITHIN SLAP EACH TIMESTEP.  THIS 
! ***  RECALCULATION AND CONVERSION HAD A SIGNIFICANT NUMERICAL OVERHEAD
! ***  ASSOCIATED WITH IT.  THIS MODULE CALCULATES A MAPPING ARRAY FOR
! ***  CONVERTING THE TRIAD FORMAT TO A COLUMN FORMAT BEFORE THE FIRST
! ***  TIMESTEP, USES IT BEFORE SENDING THE MATRIX TO THE SOLVER EACH
! ***  TIMESTEP, AND SIGNIFICANTLY REDUCES THE NUMERICAL OVERHEAD ASSOCIATED
! ***  WITH PREPARING THE Ax=b MATRICES AND VECTORS FOR SOLUTION.
! ***
! ***  THE MODULE CONTAINS THE FUNCTIONS ALLOCATETRIAD2COLUMN TO ALLOCATE
! ***  STORAGE REQUIRED FOR THE MODULE AND MAPTRIAD2COLUMN WHICH IS USED
! ***  TO APPLY THE MAPPING ARRAY BEFORE CALLING SLAP EACH TIMESTEP.
! ***  THE MODULE ALSO GIVE GLOBAL ACCESS TO THE NIA, NJA, AND NAMAP
! ***  MAPPING VECTORS AND ARRAYS.  THE MAPPING ARRAYS ARE GENERATED USING
! ***  THE DS2YWS AND QS2I1D FUNCTIONS WHICH WERE TAKEN DIRECTLY FROM THE 
! ***  SLAP SOLVER LIBRARY AND IS USED BY SLAP TO PERFORM THIS SAME OPERATION.
! ***  THE DS2YWS AND QS2I1D HAVE BEEN MODIFIED SLIGHTLY FROM THE ORIGINAL
! ***  VERSIONS IN ORDER TO CONVERT THEM FROM SUBROUTINES TO FUNCTIONS. 
! ***

module Triad2Column
    use MSErrorHandler
    use SutraMSPrecision

    INTEGER (I4B), ALLOCATABLE :: &
      nIA(:),   &
      nJA(:),   &
      nAMAP(:)
    real (DP), allocatable :: &
      rAMAP(:), &
      tA(:)

    public &
      nIA, &
      nJA, &
      nAMAP, &
      DS2Yws, &
      AllocateTriad2Column, &
      MapTriad2Column

    contains

    !allocate storage for Triad2Column modifications
    logical function AllocateTriad2Column()
      use DIMX, ONLY: NELT
      use TotalStorage
      implicit none
      lOk=.false.
      IF(.NOT.ALLOCATED(nIA)) &
        allocate( &
                  nIA(NELT), &
                  nJA(NELT), &
                  nAMAP(NELT), &
                  rAMAP(NELT), &
                  tA(NELT), &
                  stat=ios)
        if(ios/=0) goto 9999
        !calculate memory requirements
        ios=AddMemory(MemoryIndex('TRI'),'Int',NELT) !nIA(NELT)
        ios=AddMemory(MemoryIndex('TRI'),'Int',NELT) !nJA(NELT)
        ios=AddMemory(MemoryIndex('TRI'),'Int',NELT) !nAMAP(NELT)
        ios=AddMemory(MemoryIndex('TRI'),'Vec',NELT) !rAMAP(NELT)
        ios=AddMemory(MemoryIndex('TRI'),'Vec',NELT) !tA(NELT)

      lOk=.true.
09999&
      AllocateTriad2Column=lOk
      return
    end function AllocateTriad2Column

    !Map A array in 
    logical function MapTriad2Column(NELT,A)
      use SOLVI
      use SutraMSPrecision
      implicit none
      integer (I4B), INTENT(IN) :: NELT
      real (DP), INTENT(INOUT)  :: A(NELT)
      !local
      logical                   :: lOk
      integer (I4B)             :: i
      
      lOk=.false.

!      if(.not.allocated(tA)) allocate(tA(NELT))      
      do i=1,NELT
        tA(i)=A(nAMAP(i))
      end do
      A=tA
!      deallocate(tA)

      lOk=.true.

09999&
      MapTriad2Column=lOk
      return

    end function MapTriad2Column

    !----------------------------------------------------------------------!
    ! Modified routine to convert from the triad format used in SUTRA      !
    ! for the SLAP solver to a column format which is used internally      !
    ! in the SLAP solver.                                                  !
    !                                                                      !
    ! Have converted the original subroutines to logical functions for     !
    ! improved error handling and have added the SUTRA funits module       !
    ! for error messages. Have also added the module f90type in order      !
    ! to consistently define variable kinds                                !
    !                                                                      !
    !                                                                      !
    !                                                                      !
    !----------------------------------------------------------------------!

    logical function DS2Yws (N, NELT, IA, JA, A, ISYM)
    !***BEGIN PROLOGUE  DS2Y
    !***PURPOSE  SLAP Triad to SLAP Column Format Converter.
    !            Routine to convert from the SLAP Triad to SLAP Column
    !            format.
    !***LIBRARY   SLATEC (SLAP)
    !***CATEGORY  D1B9
    !***TYPE      DOUBLE PRECISION (SS2Y-S, DS2Y-D)
    !***KEYWORDS  LINEAR SYSTEM, SLAP SPARSE
    !***AUTHOR  Seager, Mark K., (LLNL)
    !             Lawrence Livermore National Laboratory
    !             PO BOX 808, L-60
    !             Livermore, CA 94550 (510) 423-3141
    !             seager@llnl.gov
    !***DESCRIPTION
    !
    ! *Usage:
    !     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
    !     DOUBLE PRECISION A(NELT)
    !
    !     CALL DS2Y( N, NELT, IA, JA, A, ISYM )
    !
    ! *Arguments:
    ! N      :IN       Integer
    !         Order of the Matrix.
    ! NELT   :IN       Integer.
    !         Number of non-zeros stored in A.
    ! IA     :INOUT    Integer IA(NELT).
    ! JA     :INOUT    Integer JA(NELT).
    ! A      :INOUT    Double Precision A(NELT).
    !         These arrays should hold the matrix A in either the SLAP
    !         Triad format or the SLAP Column format.  See "Description",
    !         below.  If the SLAP Triad format is used, this format is
    !         translated to the SLAP Column format by this routine.
    ! ISYM   :IN       Integer.
    !         Flag to indicate symmetric storage format.
    !         If ISYM=0, all non-zero entries of the matrix are stored.
    !         If ISYM=1, the matrix is symmetric, and only the lower
    !         triangle of the matrix is stored.
    !
    ! *Description:
    !       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
    !       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
    !       Column format.  The user can hand this routine either of the
    !       of these data structures.  If the SLAP Triad format is give
    !       as input then this routine transforms it into SLAP Column
    !       format.  The way this routine tells which format is given as
    !       input is to look at JA(N+1).  If JA(N+1) = NELT+1 then we
    !       have the SLAP Column format.  If that equality does not hold
    !       then it is assumed that the IA, JA, A arrays contain the
    !       SLAP Triad format.
    !
    !       =================== S L A P Triad format ===================
    !       This routine requires that the  matrix A be   stored in  the
    !       SLAP  Triad format.  In  this format only the non-zeros  are
    !       stored.  They may appear in  *ANY* order.  The user supplies
    !       three arrays of  length NELT, where  NELT is  the number  of
    !       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
    !       each non-zero the user puts the row and column index of that
    !       matrix element  in the IA and  JA arrays.  The  value of the
    !       non-zero   matrix  element is  placed  in  the corresponding
    !       location of the A array.   This is  an  extremely  easy data
    !       structure to generate.  On  the  other hand it   is  not too
    !       efficient on vector computers for  the iterative solution of
    !       linear systems.  Hence,   SLAP changes   this  input    data
    !       structure to the SLAP Column format  for  the iteration (but
    !       does not change it back).
    !
    !       Here is an example of the  SLAP Triad   storage format for a
    !       5x5 Matrix.  Recall that the entries may appear in any order.
    !
    !           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
    !                              1  2  3  4  5  6  7  8  9 10 11
    !       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
    !       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
    !       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
    !       | 0  0  0 44  0|
    !       |51  0 53  0 55|
    !
    !       =================== S L A P Column format ==================
    !
    !       This routine  requires that  the matrix A  be stored in  the
    !       SLAP Column format.  In this format the non-zeros are stored
    !       counting down columns (except for  the diagonal entry, which
    !       must appear first in each  "column")  and are stored  in the
    !       double precision array A.   In other words,  for each column
    !       in the matrix put the diagonal entry in  A.  Then put in the
    !       other non-zero  elements going down  the column (except  the
    !       diagonal) in order.   The  IA array holds the  row index for
    !       each non-zero.  The JA array holds the offsets  into the IA,
    !       A arrays  for  the  beginning  of each   column.   That  is,
    !       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
    !       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
    !       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
    !       Note that we always have  JA(N+1) = NELT+1,  where N is  the
    !       number of columns in  the matrix and NELT  is the number  of
    !       non-zeros in the matrix.
    !
    !       Here is an example of the  SLAP Column  storage format for a
    !       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
    !       column):
    !
    !           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
    !                              1  2  3    4  5    6  7    8    9 10 11
    !       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
    !       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
    !       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
    !       | 0  0  0 44  0|
    !       |51  0 53  0 55|
    !
    !***REFERENCES  (NONE)
    !***ROUTINES CALLED  QS2I1D
    !***REVISION HISTORY  (YYMMDD)
    !   871119  DATE WRITTEN
    !   881213  Previous REVISION DATE
    !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
    !   890922  Numerous changes to prologue to make closer to SLATEC
    !           standard.  (FNF)
    !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
    !   910411  Prologue converted to Version 4.0 format.  (BAB)
    !   910502  Corrected C***FIRST EXECUTABLE STATEMENT line.  (FNF)
    !   920511  Added complete declaration section.  (WRB)
    !   930701  Updated CATEGORY section.  (FNF, WRB)
    !***END PROLOGUE  DS2Y
    !     .. Modules
          USE funits
    !     .. Scalar Arguments ..
          integer (I4B), intent(IN)                     :: ISYM, N, NELT
    !     .. Array Arguments ..
          real (DP), intent(INOUT), dimension(NELT)     :: A
          integer (I4B), intent(INOUT), dimension(NELT) :: IA(NELT), JA(NELT)
    !     .. Local Scalars ..
          real (DP)                                     :: TEMP
          integer (I4B)                                 :: I, IBGN, ICOL, IEND, ITEMP, J
    !     .. Local logical
          logical                                       :: lrc  
    !     .. External Subroutines ..
    !***FIRST EXECUTABLE STATEMENT  DS2Y
    !
    !         Initialize DS2Yws to false
          DS2Yws = .false.
    !
    !         Check to see if the (IA,JA,A) arrays are in SLAP Column
    !         format.  If it's not then transform from SLAP Triad.
    !
          IF( JA(N+1).EQ.NELT+1 ) then
            write(fLST,*) 'JA(N+1) [',JA(N+1),'] = NELT+1 [',(NELT+1)
            goto 999
          end if
    !
    !         Sort into ascending order by COLUMN (on the ja array).
    !         This will line up the columns.
    !
         lrc = QS2I1Dws( JA, IA, A, NELT, 1 )
         if(.not.lrc) then
           write(fLST,*) 'Failed Sorting function QS2I1Dws'
           goto 999
         end if
    !
    !         Loop over each column to see where the column indices change
    !         in the column index array ja.  This marks the beginning of the
    !         next column.
    !
    !VD$R NOVECTOR
          JA(1) = 1
          DO 20 ICOL = 1, N-1
             DO 10 J = JA(ICOL)+1, NELT
                IF( JA(J).NE.ICOL ) THEN
                   JA(ICOL+1) = J
                   GOTO 20
                ENDIF
     10      CONTINUE
     20   CONTINUE
          JA(N+1) = NELT+1
    !
    !         Mark the n+2 element so that future calls to a SLAP routine
    !         utilizing the YSMP-Column storage format will be able to tell.
    !
          JA(N+2) = 0
    !
    !         Now loop through the IA array making sure that the diagonal
    !         matrix element appears first in the column.  Then sort the
    !         rest of the column in ascending order.
    !
          DO 70 ICOL = 1, N
             IBGN = JA(ICOL)
             IEND = JA(ICOL+1)-1
             DO 30 I = IBGN, IEND
                IF( IA(I).EQ.ICOL ) THEN
    !
    !              Swap the diagonal element with the first element in the
    !              column.
    !
                   ITEMP = IA(I)
                   IA(I) = IA(IBGN)
                   IA(IBGN) = ITEMP
                   TEMP = A(I)
                   A(I) = A(IBGN)
                   A(IBGN) = TEMP
                   GOTO 40
                ENDIF
     30      CONTINUE
     40      IBGN = IBGN + 1
             IF( IBGN.LT.IEND ) THEN
                DO 60 I = IBGN, IEND
                   DO 50 J = I+1, IEND
                      IF( IA(I).GT.IA(J) ) THEN
                         ITEMP = IA(I)
                         IA(I) = IA(J)
                         IA(J) = ITEMP
                         TEMP = A(I)
                         A(I) = A(J)
                         A(J) = TEMP
                      ENDIF
     50            CONTINUE
     60         CONTINUE
             ENDIF
     70   CONTINUE
    !
    !         Successful completion of function
          DS2Yws = .true.
    !
    !     Return to calling routine
          00999 RETURN
    !------------- LAST LINE OF DS2Y FOLLOWS ----------------------------
          end function DS2Yws


          logical function QS2I1Dws (IA, JA, A, N, KFLAG)
    !***BEGIN PROLOGUE  QS2I1D
    !***SUBSIDIARY
    !***PURPOSE  Sort an integer array, moving an integer and DP array.
    !            This routine sorts the integer array IA and makes the same
    !            interchanges in the integer array JA and the double pre-
    !            cision array A.  The array IA may be sorted in increasing
    !            order or decreasing order.  A slightly modified QUICKSORT
    !            algorithm is used.
    !***LIBRARY   SLATEC (SLAP)
    !***CATEGORY  N6A2A
    !***TYPE      DOUBLE PRECISION (QS2I1R-S, QS2I1D-D)
    !***KEYWORDS  SINGLETON QUICKSORT, SLAP, SORT, SORTING
    !***AUTHOR  Jones, R. E., (SNLA)
    !           Kahaner, D. K., (NBS)
    !           Seager, M. K., (LLNL) seager@llnl.gov
    !           Wisniewski, J. A., (SNLA)
    !***DESCRIPTION
    !     Written by Rondall E Jones
    !     Modified by John A. Wisniewski to use the Singleton QUICKSORT
    !     algorithm. date 18 November 1976.
    !
    !     Further modified by David K. Kahaner
    !     National Bureau of Standards
    !     August, 1981
    !
    !     Even further modification made to bring the code up to the
    !     Fortran 77 level and make it more readable and to carry
    !     along one integer array and one double precision array during
    !     the sort by
    !     Mark K. Seager
    !     Lawrence Livermore National Laboratory
    !     November, 1987
    !     This routine was adapted from the ISORT routine.
    !
    !     ABSTRACT
    !         This routine sorts an integer array IA and makes the same
    !         interchanges in the integer array JA and the double precision
    !         array A.
    !         The array IA may be sorted in increasing order or decreasing
    !         order.  A slightly modified quicksort algorithm is used.
    !
    !     DESCRIPTION OF PARAMETERS
    !        IA - Integer array of values to be sorted.
    !        JA - Integer array to be carried along.
    !         A - Double Precision array to be carried along.
    !         N - Number of values in integer array IA to be sorted.
    !     KFLAG - Control parameter
    !           = 1 means sort IA in INCREASING order.
    !           =-1 means sort IA in DECREASING order.
    !
    !***SEE ALSO  DS2Y
    !***REFERENCES  R. C. Singleton, Algorithm 347, An Efficient Algorithm
    !                 for Sorting With Minimal Storage, Communications ACM
    !                 12:3 (1969), pp.185-7.
    !***ROUTINES CALLED  XERMSG
    !***REVISION HISTORY  (YYMMDD)
    !   761118  DATE WRITTEN
    !   890125  Previous REVISION DATE
    !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
    !   890922  Numerous changes to prologue to make closer to SLATEC
    !           standard.  (FNF)
    !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
    !   900805  Changed XERROR calls to calls to XERMSG.  (RWC)
    !   910411  Prologue converted to Version 4.0 format.  (BAB)
    !   910506  Made subsidiary to DS2Y and corrected reference.  (FNF)
    !   920511  Added complete declaration section.  (WRB)
    !   920929  Corrected format of reference.  (FNF)
    !   921012  Corrected all f.p. constants to double precision.  (FNF)
    !***END PROLOGUE  QS2I1D
    !VD$R NOVECTOR
    !VD$R NOCONCUR
    !     .. Modules
          USE funits
    !     .. Scalar Arguments ..
          integer (I4B), intent(IN)                   :: KFLAG, N
    !     .. Array Arguments ..
          real (DP), intent(INOUT), dimension(N)      :: A(N)
          integer (I4B), intent(INOUT), dimension (N) :: IA(N), JA(N)
    !     .. Local Scalars ..
          real (DP)                                   :: R, TA, TTA
          integer (I4B)                               :: I, IIT, IJ, IT, &
                                                         J, JJT, JT,     &
                                                         K, KK, L, M, NN
    !     .. Local Arrays ..
          integer (I4B), dimension(21)                :: IL, IU
    !     .. External Subroutines ..
    !     .. Intrinsic Functions ..
          INTRINSIC ABS, INT
    !***FIRST EXECUTABLE STATEMENT  QS2I1D
    !
    !         Initialize QS2I1Dws to false
          QS2I1Dws = .false.
    !
          NN = N
          IF (NN.LT.1) THEN
    !         CALL XERMSG ('SLATEC', 'QS2I1D',
    !     $      'The number of values to be sorted was not positive.', 1, 1)
    !         RETURN
             write (fLST,*) 'The number of values to be sorted was not positive.'
             goto 999
          ENDIF
    !      IF( N.EQ.1 ) RETURN
          IF( N.EQ.1 ) then
             write (fLST,*) 'N = 1'
             goto 998
          end if
          KK = ABS(KFLAG)
          IF ( KK.NE.1 ) THEN
    !         CALL XERMSG ('SLATEC', 'QS2I1D',
    !     $      'The sort control parameter, K, was not 1 or -1.', 2, 1)
    !         RETURN
             write (fLST,*) 'The sort control parameter, K, was not 1 or -1.'
             goto 999
          ENDIF
    !
    !     Alter array IA to get decreasing order if needed.
    !
          IF( KFLAG.LT.1 ) THEN
             DO 20 I=1,NN
                IA(I) = -IA(I)
     20      CONTINUE
          ENDIF
    !
    !     Sort IA and carry JA and A along.
    !     And now...Just a little black magic...
          M = 1
          I = 1
          J = NN
          R = .375D0
     210  IF( R.LE.0.5898437D0 ) THEN
             R = R + 3.90625D-2
          ELSE
             R = R-.21875D0
          ENDIF
     225  K = I
    !
    !     Select a central element of the array and save it in location
    !     it, jt, at.
    !
          IJ = I + INT ((J-I)*R)
          IT = IA(IJ)
          JT = JA(IJ)
          TA = A(IJ)
    !
    !     If first element of array is greater than it, interchange with it.
    !
          IF( IA(I).GT.IT ) THEN
             IA(IJ) = IA(I)
             IA(I)  = IT
             IT     = IA(IJ)
             JA(IJ) = JA(I)
             JA(I)  = JT
             JT     = JA(IJ)
             A(IJ)  = A(I)
             A(I)   = TA
             TA     = A(IJ)
          ENDIF
          L=J
    !
    !     If last element of array is less than it, swap with it.
    !
          IF( IA(J).LT.IT ) THEN
             IA(IJ) = IA(J)
             IA(J)  = IT
             IT     = IA(IJ)
             JA(IJ) = JA(J)
             JA(J)  = JT
             JT     = JA(IJ)
             A(IJ)  = A(J)
             A(J)   = TA
             TA     = A(IJ)
    !
    !     If first element of array is greater than it, swap with it.
    !
             IF ( IA(I).GT.IT ) THEN
                IA(IJ) = IA(I)
                IA(I)  = IT
                IT     = IA(IJ)
                JA(IJ) = JA(I)
                JA(I)  = JT
                JT     = JA(IJ)
                A(IJ)  = A(I)
                A(I)   = TA
                TA     = A(IJ)
             ENDIF
          ENDIF
    !
    !     Find an element in the second half of the array which is
    !     smaller than it.
    !
      240 L=L-1
          IF( IA(L).GT.IT ) GO TO 240
    !
    !     Find an element in the first half of the array which is
    !     greater than it.
    !
      245 K=K+1
          IF( IA(K).LT.IT ) GO TO 245
    !
    !     Interchange these elements.
    !
          IF( K.LE.L ) THEN
             IIT   = IA(L)
             IA(L) = IA(K)
             IA(K) = IIT
             JJT   = JA(L)
             JA(L) = JA(K)
             JA(K) = JJT
             TTA   = A(L)
             A(L)  = A(K)
             A(K)  = TTA
             GOTO 240
          ENDIF
    !
    !     Save upper and lower subscripts of the array yet to be sorted.
    !
          IF( L-I.GT.J-K ) THEN
             IL(M) = I
             IU(M) = L
             I = K
             M = M+1
          ELSE
             IL(M) = K
             IU(M) = J
             J = L
             M = M+1
          ENDIF
          GO TO 260
    !
    !     Begin again on another portion of the unsorted array.
    !
      255 M = M-1
          IF( M.EQ.0 ) GO TO 300
          I = IL(M)
          J = IU(M)
      260 IF( J-I.GE.1 ) GO TO 225
          IF( I.EQ.J ) GO TO 255
          IF( I.EQ.1 ) GO TO 210
          I = I-1
      265 I = I+1
          IF( I.EQ.J ) GO TO 255
          IT = IA(I+1)
          JT = JA(I+1)
          TA =  A(I+1)
          IF( IA(I).LE.IT ) GO TO 265
          K=I
      270 IA(K+1) = IA(K)
          JA(K+1) = JA(K)
          A(K+1)  =  A(K)
          K = K-1
          IF( IT.LT.IA(K) ) GO TO 270
          IA(K+1) = IT
          JA(K+1) = JT
          A(K+1)  = TA
          GO TO 265
    !
    !     Clean up, if necessary.
    !
      300 IF( KFLAG.LT.1 ) THEN
             DO 310 I=1,NN
                IA(I) = -IA(I)
     310     CONTINUE
          ENDIF
    !
    !     Successful completion of function
          00998 QS2I1Dws = .true.
    !
    !     Return to calling routine
    !
          00999 RETURN
    !------------- LAST LINE OF QS2I1D FOLLOWS ----------------------------
          end function QS2I1Dws

end module Triad2Column