      MODULE SolverStorage
! ***
! ***  CONTAINS FUNCTION TO ALLOCATE STORAGE FOR SOLVER
! ***  (ALLOCATESOLVERSTORAGE).
! ***  
! ***  ALLOCATE SOLVERSTORAGE IS USED FOR ALL ITERATIVE SOLVERS.
! ***
        USE SOLVI 
        USE SutraMSPrecision
        USE ColumnStorage

        IMPLICIT NONE

        PUBLIC :: &
               AllocateSolverStorage 
        CONTAINS

          !generic function for allocating storage for iterative solvers
          !called from SutraStorage::AllocateMainStorage
          logical function AllocateSolverStorage()
            use SutraStorage
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
            AllocateSolverStorage=lOk
            return
          end function AllocateSolverStorage

      END MODULE SolverStorage
