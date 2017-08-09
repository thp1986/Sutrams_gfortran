    module ColumnStorage
      use SutraMSPrecision
      use DIMS
      use DIMX
      USE FUNITS
      USE CONTRL
      USE TotalStorage
      USE MSErrorHandler
      implicit none

!.....DEFINE DERIVED TYPE LNKLST (LINKED LIST) WITH TWO COMPONENTS:      
!        NODNUM (NODE NUMBER) AND NENT (POINTER TO NEXT ENTRY).
      TYPE LNKLST
         INTEGER (I4B) :: NODNUM
         TYPE (LNKLST), POINTER :: NENT
      END TYPE LNKLST

!.....DECLARE DENT, DENTPV, DENTPI, AND DENTNW AS GENERAL-PURPOSE
!        POINTERS OF TYPE LNKLST.
      TYPE (LNKLST), POINTER :: DENT, DENTPV, DENTPI, DENTNW

!.....DEFINE DERIVED TYPE IPOINT WITH ONE COMPONENT: A POINTER, PL,
!        OF TYPE LNKLST.
      TYPE IPOINT
         TYPE (LNKLST), POINTER :: PL
      END TYPE IPOINT

!.....DECLARE HLIST, AN ARRAY OF POINTERS THAT WILL POINT TO THE HEAD
!        OF THE LINKED LIST OF NEIGHBORS FOR EACH NODE.
      TYPE (IPOINT), ALLOCATABLE :: HLIST(:)

      INTEGER (I4B), ALLOCATABLE :: LLIST(:)

      !misc items
!      LOGICAL :: lColumnStorage=.false.
      LOGICAL :: lColumnStorage=.true.

      public :: &
        lColumnStorage, &
        SetSolverDimensions, &
        ColAllocateSolverStorage, &
        MakeITRIJTRI


      contains

      !generic function for setting/resetting solver dimensions
      logical function SetSolverDimensions()
        use SOLVI
        use ITSOLI
!        logical       :: lOk
        integer (I4B) :: &
          NELTA, &
          NWIP, &
          NWFP, &
          NWIU, &
          NWFU, &
          KMXSOLVU, &
          NMXSAVEU

        lOk=.false.
!........SET DIMENSIONS FOR ITERATIVE SOLVER(S)
         NCBI = 1
         NBIX = NBI
         NELTA = NELT
         !KSOLVR = KSOLVP
         !NSAVE = NSAVEP
         CALL DIMWRK(KSOLVP, NSAVEP, NN, NELTA, NWIP, NWFP)               
!........MAXIMUM DIMENSION OF U ITERATIVE SOLUTION                      
         KMXSOLVU = MAXVAL(KSOLVU)
         NMXSAVEU = MAXVAL(NSAVEU)
         CALL DIMWRK(KMXSOLVU, NMXSAVEU, NN, NELTA, NWIU, NWFU)
         !KSOLVR = KSOLVU                                                 
         !NSAVE = NSAVEU                                                  
         !CALL DIMWRK(KSOLVR, NSAVE, NN, NELTA, NWIU, NWFU)               
         NWI = MAX(NWIP, NWIU)                                           
         NWF = MAX(NWFP, NWFU)                                           
        lOk=.true.
09999&
        SetSolverDimensions=lOk
        return
      end function SetSolverDimensions


      !generic function for allocating storage for iterative solvers
      !called from SutraStorage::AllocateMainStorage
      logical function ColAllocateSolverStorage()
        use SutraStorage
!        logical       :: lOk
        integer (I4B) :: ialloerr
        integer (I4B) :: NNP1

        lOk=.false.
        NNP1=NELT
        if(lColumnStorage) NNP1=NN+1
        Allocate( &
                  PMAT(NELT, NCBI), &
                  UMAT(NELT,NCBI), &
                  ITRI(NELT), &
                  JTRI(NNP1), &
                  NBI27(NBIX), &
                  IWK(NWI), &
                  IPARM(NIPARM, 2), &
                  FWK(NWF), &
                  FPARM(NFPARM, 2), &
                  B(NNNX), &
                  stat=ios)
        if (ios/=0) goto 9999
        !Calculate memory requirements
        ios=AddMemory(MemoryIndex('SOL'),'Arr',NELT*NCBI) !PMAT(NELT, NCBI)
        ios=AddMemory(MemoryIndex('SOL'),'Arr',NELT*NCBI) !UMAT(NELT,NCBI)
        ios=AddMemory(MemoryIndex('SOL'),'Int',     NELT) !ITRI(NELT)
        ios=AddMemory(MemoryIndex('SOL'),'Int',     NNP1) !JTRI(NNP1)
        ios=AddMemory(MemoryIndex('SOL'),'Int',     NBIX) !NBI27(NBIX)
        ios=AddMemory(MemoryIndex('SOL'),'Int',      NWI) !IWK(NWI)
        ios=AddMemory(MemoryIndex('SOL'),'Irr', NIPARM*2) !IPARM(NIPARM, 2)
        ios=AddMemory(MemoryIndex('SOL'),'Vec',      NWF) !FWK(NWF)
        ios=AddMemory(MemoryIndex('SOL'),'Arr', NFPARM*2) !FPARM(NFPARM, 2)
        ios=AddMemory(MemoryIndex('SOL'),'Vec',     NNNX) !B(NNNX)

        lOk=.true.
09999&
        ColAllocateSolverStorage=lOk
        return
      end function ColAllocateSolverStorage

      !function to set IA JA arrays in column storage format
      logical function MakeITRIJTRI()
        use SutraStorage
!        logical :: lOk
        integer (I4B) :: &
          I, K, L, &
          IL, IC, &
          JL, JC
!        integer (I4B) :: &				!MT: Commented out due to error occur during compiling by gfortran (conflict with "N48" in module "dimx")
!          N48
        integer (I4B) :: &
          JTRISTRT


        lOk=.false.
!
!.....SET UP POINTER ARRAYS ITRI AND JTRI THAT SPECIFY MATRIX STRUCTURE IN
!        "SLAP COLUMN" FORMAT.  FOR EACH NODE, CONSTRUCT A LINKED LIST
!        OF NEIGHBORING NODES.  HLIST(K) POINTS TO THE HEAD OF THE LIST
!        FOR NODE K.  THEN, TRANSFER THE LISTS TO ARRAYS ITRI AND JTRI.
!
!.....ALLOCATE HLIST AND LLIST, AND INITITRILIZE LIST LENGTHS TO ZERO.
      ALLOCATE(LLIST(NN), HLIST(NN))
      DO 490 I=1,NN
         ALLOCATE(HLIST(I)%PL)
         LLIST(I) = 0
  490 CONTINUE
!.....3-D MESH                                                          
      IF(IABS(KTYPE).EQ.3)THEN 
         N48 = 8 
!.....2-D MESH                                                          
      ELSE 
         N48 = 4 
      ENDIF 
!.....LOOP THROUGH INCIDENCE LIST.
      DO 500 L=1,NE
      DO 500 IL=1,N48
         IC = IN((L-1)*N48+IL)
      DO 500 JL=1,N48
         JC = IN((L-1)*N48+JL)
!........INSERT NEIGHBOR JC IN LIST FOR NODE IC IN ASCENDING ORDER.
!           (IF DUPLICATE OR SELF-NEIGHBOR, SKIP IT.)
         IF (JC.EQ.IC) THEN
!...........SKIP SELF-NEIGHBOR.
            GOTO 500
         ELSE IF (LLIST(IC).EQ.0) THEN
!...........PLACE FIRST LIST ENTRY AT HEAD.
            HLIST(IC)%PL%NODNUM = JC
            GOTO 498
         ELSE
!...........INSERT INTO LIST, OR SKIP IF DUPLICATE.
            ALLOCATE(DENTPV)
            DENTPI => DENTPV
            DENTPV%NENT => HLIST(IC)%PL
            DO 495 K=1,LLIST(IC)
               DENT => DENTPV%NENT
               IF (JC.EQ.DENT%NODNUM) THEN
                  DEALLOCATE(DENTPI)
                  GOTO 500
               ELSE IF (JC.LT.DENT%NODNUM) THEN
                  ALLOCATE(DENTNW)
                  DENTNW%NODNUM = JC
                  DENTNW%NENT => DENT
                  IF (K.EQ.1) THEN
                     HLIST(IC)%PL => DENTNW
                  ELSE
                     DENTPV%NENT => DENTNW
                  END IF
                  DEALLOCATE(DENTPI)
                  GOTO 498
               END IF
               DENTPV => DENT
  495       CONTINUE
!...........APPEND TO TAIL.
            ALLOCATE(DENTNW)
            DENTNW%NODNUM = JC
            DENT%NENT => DENTNW
            DEALLOCATE(DENTPI)
         END IF
  498    LLIST(IC) = LLIST(IC) + 1
  500 CONTINUE
!.....COMPUTE THE ARRAY DIMENSION NELT AND ALLOCATE ARRAY ITRI.
      NELT = 0
      DO 600 I=1,NN
  600    NELT = NELT + LLIST(I) + 1

      MSErrorValue%cDataSet='CS1'
      if (.not.SetSolverDimensions()) &
        call ErrorIO('SetSolverDimensions:: Could not set solver storage sizes')

      MSErrorValue%cDataSet='CS2'
      if(.not.ColAllocateSolverStorage()) &
        call ErrorIO('ColAllocateSolverStorage:: Could not allocate storage for iterative solvers (are they already allocated??)')

!.....TRANSFER THE LINKED LISTS TO ARRAYS ITRI AND JTRI IN SLAP COLUMN
!        FORMAT.  DEALLOCATE POINTERS AS THEY ARE TRANSFERRED.
      JTRISTRT = 1
      DO 660 I=1,NN
         JTRI(I) = JTRISTRT
         ITRI(JTRISTRT) = I
         DENT => HLIST(I)%PL
         DO 650 K=1,LLIST(I)
            ITRI(JTRISTRT + K) = DENT%NODNUM
            DENTPV => DENT
            DENT => DENT%NENT
            DEALLOCATE(DENTPV)
  650    CONTINUE
         JTRISTRT = JTRISTRT + LLIST(I) + 1
  660 CONTINUE
      JTRI(NN + 1) = NELT + 1
      DEALLOCATE(HLIST, LLIST)
!
  700 CONTINUE

      if(lDebugData) then
        write (fLST,'(//a)') 'Column Format'
        write (fLST,'(a2,20(1x,i6),1000(:/2x,20(1x,i6)))') &
          '  ',(i,i=1,NELT)
        write (fLST,'(a2,20(1x,i6),1000(:/2x,20(1x,i6)))') &
          'IA',(ITRI(i),i=1,NELT)
        write (fLST,'(a2,20(1x,i6),1000(:/2x,20(1x,i6)))') &
          'JA',(JTRI(i),i=1,NN+1)
        write (fLST,*) ' '
      end if

        lOk=.true.
09999&
        MakeITRIJTRI=lOk
        return
      end function MakeITRIJTRI

    end module ColumnStorage