!     MODULE            SUTRASTORAGE           SUTRA-MS VERSION 2004.1
                                                                       
! ***
! *** PURPOSE :                                                         
! ***  MODULE THAT DEFINES ALL ALLOCATABLE VARIABLES FOR SUTRA-MS AND
! ***  CONTAINS FUNCTION TO ALLOCATE INITIAL AND PRIMARY STORAGE FOR
! ***  SUTRA-MS (ALLOCATEINITIALSTORAGE AND ALLOCATEMAINSTORAGE).
! ***

     module SutraStorage
       USE DIMS
       USE DIMX
       USE OBS
       USE FUNITS
       USE ATSDATA
       USE PARAMS
       USE TotalStorage
       USE MSErrorHandler
       use SutraMSPrecision 
       implicit none
       !SUTRA MS Type definitions
         !Specified Pressure structure
           type tSpecifiedP
            integer (I4B)             :: node
            real (DP)                 :: P
            real (DP), allocatable    :: U(:)
           end type tSpecifiedP
         !Specified Concentration/Temperature structure
           type tSpecifiedU
            integer (I4B) :: node
            real (DP)     :: U
           end type tSpecifiedU
         !Specified U wrapper type
           type tMultiSpeciesBC
            type (tSpecifiedU), allocatable :: SpecifiedU(:)
           end type tMultiSpeciesBC
         !Specified Pressure class
           type (tSpecifiedP), allocatable     :: SpecifiedPBC(:)
         !Specified Multi-species class for specified U
           type (tMultiSpeciesBC), allocatable :: MultiSpeciesBC(:)
       !standard SUTRA dimensions
         integer (I4B), allocatable :: &
           IQSOP(:),  IPBC (:),    IN (:)
         real (DP), allocatable :: &
             QIN(:),    PBC(:), QPLITR(:), &
               X(:),      Y(:),      Z(:), &
              SW(:),  DSWDP(:),    RHO(:),  &
            PVEL(:), &
            VMAG(:),  VANG1(:),  VANG2(:), &
          GXSI(:,:), GETA(:,:), GZET(:,:), &
             VOL(:), PMAT(:,:),   PVEC(:), UMAT(:,:), &
             PM1(:),  PITER(:),   RCIT(:), RCITM1(:), & 
             VISCO(:)         								!MT 
       !SMS DIMENSIONS
         integer(I4B), allocatable :: &
           IQSOU(:,:), IUBC(:,:), &
           NOCONU(:), &
           IOBS(:), ITOBS(:), &
           NBI27(:) , &
           IWK(:), IPARM(:,:), &
           ITRI(:), JTRI(:), &
           MIOFF(:)
         real(DP), allocatable :: &
           UITER(:,:), UM1(:,:), UM2(:,:),       &
           SL(:), SR(:), CS1(:,:), CS2(:), CS3(:), UIN(:,:), &
           QUIN(:,:), UVEC(:,:), CVEC(:), CITER(:), &
           UBC(:,:), &
           RU(:), &
           FWK(:), FPARM(:,:), &
           B(:)
       !SMS DIMENSIONS                                                   

       !ARRAYS TO STORE DATA FROM LAST SUCCESSFUL ITERATION              
       !ONLY ALLOCATED IF THE AUTOMATIC TIMESTEP ALGORITHM IS USED       
         REAL (DP), ALLOCATABLE :: &
             PVEC_(:), &
             UVEC_(:,:), &
             PM1_(:), &
             UM1_(:,:), &
             RCIT_(:), &
             CS1_(:,:), &
             SW_(:), &
             QPLITR_(:), &            
             VISCO_(:)

         public &
           IQSOP, IPBC, IN, &

           QIN, PBC, QPLITR, &
           X, Y, Z, SW, DSWDP, RHO,  &
           PVEL, &
           VMAG, VANG1, VANG2, &
           GXSI, GETA, GZET, &
           VOL, PMAT, PVEC, UMAT, &
           PM1, PITER, RCIT, RCITM1, &

           IQSOU, IUBC, &
           NOCONU, &
           IOBS, ITOBS, &
           NBI27 , &
           IWK, IPARM, &
           ITRI, JTRI, &
           MIOFF, &

           UITER, UM1, UM2,       &
           SL, SR, CS1, CS2, CS3, UIN, &
           QUIN, UVEC, CVEC, CITER, &
           UBC, &
           RU, &
           FWK, FPARM, &
           B, &

           PVEC_, &
           UVEC_, &
           PM1_, &
           UM1_, &
           RCIT_, &
           CS1_, &
           SW_, &
           QPLITR_, &
           VISCO_, &
!
           SpecifiedPBC, &
           MultiSpeciesBC, &
!
           AllocateInitialStorage, &
           AllocateMainStorage

           contains
             !function to allocate storage for 
             !ITERAT, PARAMS, CONTRL, SOLVI, ITSOLI, ITSOLR, JCOLS, and MODSPR modules
             logical function AllocateInitialStorage()
               use SutraMSPrecision
               use ITERAT 
               use PARAMS 
               use CONTRL 
               use SOLVI 
               use ITSOLR 
               use ITSOLI 
               use JCOLS 
               use MODSOR 
               IMPLICIT NONE
               !local
               REAL (DP)     :: DZERO = 0.0D0 
               REAL (DP)     :: DONE = 1.0D0 
               INTEGER (I4B) :: NZERO = 0 
               lOk=.false.

                NCOLMX = 10 + NSPE-1 !for JCOLS - 9 CHANGED TO 10 FOR VERSION 1.1
!                NCOLMX = 11 + NSPE-1 !for JCOLS - 9 CHANGED TO 10 FOR VERSION 1.1	!HYDCON 01/6/2017
                !ALLOCATE MODULE ARRAYS
                !*** Add additional allocated variables in memory calculations below
                allocate ( &
                          !ITERAT
                          RUMAX(NSPE), &
                          RUM(NSPE), & 
                          IUWORS(NSPE), & 
                          !PARAMS
                          SPNAME(NSPE) ,&
                          VISC0(NSPE), & 
                          PRODF1(NSPE), PRODS1(NSPE), & 
                          PRODF0(NSPE), PRODS0(NSPE), & 
                          DRWDU(NSPE), SIGMAW(NSPE), URHOW0(NSPE), & 
                          CHI1(NSPE), CHI2(NSPE), & 
                          ATSPMULT(NSPE), & 
                          NSOU(NSPE), NUBC(NSPE), & 
                          !CONTRL
                          GNUU(NSPE), & 
                          !SOLVI
                          KSOLVU(NSPE), & 
                          U_ML(NSPE), &
                          StartUTime(NSPE), &
                          !ITSOLI
                          ITRMXU(NSPE), ITOLU(NSPE), NSAVEU(NSPE), & 
                          !ITSOLR
                          TOLU(NSPE), & 
                          !JCOLS                                                             
                          CSMSK5(NCOLMX), &               !char
                          K5COL(NCOLMX), &                !char
                          K6COL(NCOLMX), &                !char
                          J5COL(NCOLMX), J6COL(NCOLMX), & !int
                          COLTK5(7 + NSPE-1), &           !char
                          !MODSPR                                                            
                          ADSMOD(NSPE), &                 !char
                          stat=ios)
                if (ios /= 0) goto 9999
                !*** Add additional allocated variables in memory calculations below

                !Calculate Memory requirements
                !ITERAT
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !RUMAX(NSPE)
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !RUM(NSPE) 
                ios=AddMemory(MemoryIndex('INI'),'Int',NSPE) !IUWORS(NSPE) 
                !PARAMS
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !VISC0(NSPE) 
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !PRODF1(NSPE)
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !PRODS1(NSPE) 
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !PRODF0(NSPE)
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !PRODS0(NSPE) 
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !DRWDU(NSPE)
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !SIGMAW(NSPE)
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !URHOW0(NSPE) 
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !CHI1(NSPE)
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !CHI2(NSPE) 
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !ATSPMULT(NSPE) 
                ios=AddMemory(MemoryIndex('INI'),'Int',NSPE) !NSOU(NSPE)
                ios=AddMemory(MemoryIndex('INI'),'Int',NSPE) !NUBC(NSPE) 
                !CONTRL
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !GNUU(NSPE) 
                !SOLVI
                ios=AddMemory(MemoryIndex('INI'),'Int',NSPE) !KSOLVU(NSPE) 
                ios=AddMemory(MemoryIndex('INI'),'Int',NSPE) !U_ML(NSPE)
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !StartUTime(NSPE)
                !ITSOLI
                ios=AddMemory(MemoryIndex('INI'),'Int',NSPE) !ITRMXU(NSPE)
                ios=AddMemory(MemoryIndex('INI'),'Int',NSPE) !ITOLU(NSPE)
                ios=AddMemory(MemoryIndex('INI'),'Int',NSPE) !NSAVEU(NSPE) 
                !ITSOLR
                ios=AddMemory(MemoryIndex('INI'),'Vec',NSPE) !TOLU(NSPE) 
                !JCOLS                                                             
                ios=AddMemory(MemoryIndex('INI'),'Chr',len(CSMSK5(1))*NCOLMX) !CSMSK5(NCOLMX)
                ios=AddMemory(MemoryIndex('INI'),'Chr',len( K5COL(1))*NCOLMX) !K5COL(NCOLMX)
                ios=AddMemory(MemoryIndex('INI'),'Chr',len( K6COL(1))*NCOLMX) !K6COL(NCOLMX)
                ios=AddMemory(MemoryIndex('INI'),'Int',NCOLMX) !J5COL(NCOLMX)
                ios=AddMemory(MemoryIndex('INI'),'Int',NCOLMX) !J6COL(NCOLMX)
                ios=AddMemory(MemoryIndex('INI'),'Chr',len(COLTK5(1))*NCOLMX) !COLTK5(7 + NSPE-1)
                !MODSPR                                                            
                ios=AddMemory(MemoryIndex('INI'),'Chr',len(ADSMOD(1))*NCOLMX) !ADSMOD(NSPE)

                !initialize all arrays to 0 or blank, &
                ! except for atspmult which is initialized to one (1)
                SPNAME = '                    ' 
                RUMAX = DZERO 
                RUM = DZERO 
                IUWORS = NZERO 
                VISC0 = DZERO 
                PRODF1 = DZERO 
                PRODS1 = DZERO 
                PRODF0 = DZERO 
                PRODS0 = DZERO 
                DRWDU = DZERO 
                SIGMAW = DZERO 
                URHOW0 = DZERO 
                CHI1 = DZERO 
                CHI2 = DZERO 
                ATSPMULT = DONE 
                NSOU = NZERO 
                NUBC = NZERO 
                GNUU = DZERO 
                KSOLVU = NZERO 
                ITRMXU = NZERO 
                ITOLU = NZERO 
                NSAVEU = NZERO 
                TOLU = DZERO 
                U_ML = NZERO
                StartUTime = -1D0
          !.....JCOLS                                                             
                CSMSK5 = '                         ' 
                K5COL = ' ' 
                K6COL = '  ' 
                J5COL = NZERO 
                J6COL = NZERO 
                ADSMOD = '          ' 
                COLTK5 = '               ' 

               lOk=.true.
09999&
               AllocateInitialStorage=lOk
               return
             end function AllocateInitialStorage


             !function for allocating main SUTRA arrays
             logical function AllocateMainStorage()
               integer (I4B) :: k, n
               lOk=.false.
!...............Allocate Dynamic Storage
                !standard SUTRA dimensions
                !*** Add additional allocated variables in memory calculations below
                allocate( &
                  IQSOP(NSOP), IPBC (NBCN), IN (NIN), &
!
                  QIN(NN), PBC(NBCN), QPLITR(NBCN), &
                  X(NN), Y(NN), Z(NN), SW(NN), DSWDP(NN), RHO(NN),  &
                  PVEL(NN), &
                  VMAG(NE), VANG1(NE), VANG2(NEX), &
                  GXSI(NE, N48), GETA(NE, N48), GZET(NEX, N48), &
                  VOL(NN), PVEC(NN), &
                  PM1(NN), PITER(NN), RCIT(NN), RCITM1(NN), &
				  VISCO(NN), &										
                !start SMS DIMENSIONS
                  IQSOU(MNSOU, NSPE), IUBC(NBCN, NSPE), &
                  NOCONU(NSPE), &
                  IOBS(NOBSN), ITOBS(NTOBSN), &
                  MIOFF(27), &
!
                  UITER(NN, NSPE), UM1(NN, NSPE), UM2(NN, NSPE),       &
                  SL(NN), SR(NN), CS1(NN, NSPE), CS2(NN), CS3(NN), UIN(NN,NSPE), &
                  QUIN(NN, NSPE), UVEC(NN, NSPE), CVEC(NN), CITER(NN), &
                  UBC(NBCN, NSPE), &
                  RU(NSPE),&
                !end SMS DIMENSIONS                                                   
                    stat=ios) 
                if(ios/=0) goto 9999
                !*** Add additional allocated variables in memory calculations below

                !Calculate Memory requirements
                ios=AddMemory(MemoryIndex('SUT'),'Int',      NSOP) !IQSOP(NSOP)
                ios=AddMemory(MemoryIndex('SUT'),'Int',      NBCN) !IPBC (NBCN)
                ios=AddMemory(MemoryIndex('SUT'),'Int',       NIN) !IN (NIN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !QIN(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',      NBCN) !PBC(NBCN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',      NBCN) !QPLITR(NBCN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !X(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !Y(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !Z(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !SW(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !DSWDP(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !RHO(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !PVEL(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NE) !VMAG(NE)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NE) !VANG1(NE)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',       NEX) !VANG2(NEX)
                ios=AddMemory(MemoryIndex('SUT'),'Arr',    NE*N48) !GXSI(NE, N48)
                ios=AddMemory(MemoryIndex('SUT'),'Arr',    NN*N48) !GETA(NE, N48)
                ios=AddMemory(MemoryIndex('SUT'),'Arr',    NN*N48) !GZET(NEX, N48)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !VOL(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !PVEC(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !PM1(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !PITER(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !RCIT(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !RCITM1(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !VISCO(NN)               
                !start SMS DIMENSIONS
                ios=AddMemory(MemoryIndex('SUT'),'Irr',MNSOU*NSPE) !IQSOU(MNSOU, NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Irr', NBCN*NSPE) !IUBC(NBCN, NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Int',      NSPE) !NOCONU(NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Int',     NOBSN) !IOBS(NOBSN)
                ios=AddMemory(MemoryIndex('SUT'),'Int',    NTOBSN) !ITOBS(NTOBSN)
                ios=AddMemory(MemoryIndex('SUT'),'Int',        27) !MIOFF(27)
                ios=AddMemory(MemoryIndex('SUT'),'Arr',   NN*NSPE) !UITER(NN, NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Arr',   NN*NSPE) !UM1(NN, NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Arr',   NN*NSPE) !UM2(NN, NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !SL(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !SR(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Arr',   NN*NSPE) !CS1(NN, NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !CS2(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !CS3(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Arr',   NN*NSPE) !UIN(NN,NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Arr',   NN*NSPE) !QUIN(NN, NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Arr',   NN*NSPE) !UVEC(NN, NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !CVEC(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',        NN) !CITER(NN)
                ios=AddMemory(MemoryIndex('SUT'),'Arr', NBCN*NSPE) !UBC(NBCN, NSPE)
                ios=AddMemory(MemoryIndex('SUT'),'Vec',      NSPE) !RU(NSPE)

                !Allocate arrays for automatic time step storage, of required
                if (fATS.GT.0 .and. .not.LContinue) then
                  !*** Add additional allocated variables in memory calculations below
                  allocate ( &
                              PVEC_(NN), &
                              UVEC_(NN, NSPE), &
                              PM1_(NN), &
                              UM1_(NN, NSPE), &
                              RCIT_(NN), &
                              CS1_(NN, NSPE), &
                              SW_(NN), &
                              QPLITR_(NBCN), &
							  VISCO_(NN), &
                              stat=ios ) 
                  !*** Add additional allocated variables in memory calculations below
                  if(ios/=0) goto 9999
                  !calculate memory requirements
                  ios=AddMemory(MemoryIndex('ATS'),'Vec',        NN) !PVEC_(NN)
                  ios=AddMemory(MemoryIndex('ATS'),'Arr',   NSPE*NN) !UVEC_(NN, NSPE)
                  ios=AddMemory(MemoryIndex('ATS'),'Vec',        NN) !PM1_(NN)
                  ios=AddMemory(MemoryIndex('ATS'),'Arr',   NSPE*NN) !UM1_(NN, NSPE)
                  ios=AddMemory(MemoryIndex('ATS'),'Vec',        NN) !RCIT_(NN)
                  ios=AddMemory(MemoryIndex('ATS'),'Arr',   NSPE*NN) !CS1_(NN, NSPE)
                  ios=AddMemory(MemoryIndex('ATS'),'Vec',        NN) !SW_(NN)
                  ios=AddMemory(MemoryIndex('ATS'),'Vec',      NBCN) !QPLITR_(NBCN)
                  ios=AddMemory(MemoryIndex('ATS'),'Vec',      NN)   !VISCO_(NN)				  
                end if
                !allocate storage for specified pressure BC data
                allocate (SpecifiedPBC(NPBC),stat=ios)
                if(ios/=0) goto 9999
                !update memory required for specified pressure data
                ios=AddMemory(MemoryIndex('SUT'),'Int',   NPBC) !MultiSpeciesBC()%node
                ios=AddMemory(MemoryIndex('SUT'),'Vec',   NPBC) !MultiSpeciesBC()%P
                do n=1,NPBC
                  allocate (SpecifiedPBC(n)%U(NSPE),stat=ios)
                  if(ios/=0) goto 9999
                end do
                ios=AddMemory(MemoryIndex('SUT'),'Vec', NSPE*NPBC) !MultiSpeciesBC()%U
                !allocate storage for multispecies BC data
                allocate (MultiSpeciesBC(NSPE),stat=ios)
                if(ios/=0) goto 9999
                do k=1,NSPE
                  if(NUBC(k)<1) cycle
                  allocate(MultiSpeciesBC(k)%SpecifiedU(NUBC(k)),stat=ios)
                  if(ios/=0) goto 9999
                  !update memory required for specified concentration data
                  ios=AddMemory(MemoryIndex('SUT'),'Int',   NUBC(k)) !MultiSpeciesBC()%node
                  ios=AddMemory(MemoryIndex('SUT'),'Vec',   NUBC(k)) !MultiSpeciesBC()%U
                end do

               lOk=.true.
09999&
               AllocateMainStorage=lOk
               return
             end function AllocateMainStorage

     end module SutraStorage
