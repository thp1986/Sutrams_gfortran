!     MODULE            TOTALSTORAGE           SUTRA-MS VERSION 2004.1
                                                                       
! ***
! *** PURPOSE :                                                         
! ***  MODULE THAT IS USED TO CALCULATE ALLOCATABLE STORAGE USED BY
! ***  EACH SUBROUTINE/FUNCTION (ADDMEMORY), CALCULATE THE TOTAL 
! ***  ALLOCATABLE STORAGE USED BY SUTRA-MS (CALCULATETOTALSTORAGE),
! ***  AND DETERMINE THE APPROPRIATE DATA STRUCTURE (MEMORY STRUCTURE)
! ***  FOR A PARTICULAR OPTION USING A SPECIFIED SET OF KEYWORDS
! ***  (MEMORYINDEX).
! ***

    module TotalStorage
      USE MSErrorHandler
      use SutraMSPrecision
      implicit none
      integer (I4B), parameter :: SizeDP=8, SizeSP=4, SizeInt=4, SizeChar=1
      integer (I4B), parameter :: nEntries=10
      integer (I4B), parameter :: IDLen=3
      integer (I4B), parameter :: TitleLen=34
      type :: tMemory
        character (len=IDLen)    :: ID
        character (len=TitleLen) :: Title
        logical :: IsUsed
        integer (I4B) :: nID
        integer (I4B) :: Arr
        integer (I4B) :: Vec
        integer (I4B) :: Int
        integer (I4B) :: Irr
        integer (I4B) :: Chr
      end type tMemory
      integer (I4B) :: &
        SumMemory=0
      type (tMemory), dimension(nEntries) :: Memory = &
        (/ tMemory('INI', '             Initial SUTRA Storage', .false.,  1, 0, 0, 0, 0, 0), &
           tMemory('SUT', '             Primary SUTRA Storage', .false.,  2, 0, 0, 0, 0, 0), &
           tMemory('ATS', '       Automatic Time Step Storage', .false.,  3, 0, 0, 0, 0, 0), &
           tMemory('SOL', '              SUTRA solver Storage', .false.,  4, 0, 0, 0, 0, 0), &
           tMemory('OBS', '       Standard SUTRA Obs. Storage', .false.,  5, 0, 0, 0, 0, 0), &
           tMemory('SOB', '   Specified Obs. Location Storage', .false.,  6, 0, 0, 0, 0, 0), &
           tMemory('ZON', '  Hydraulic Parameter Zone Storage', .false.,  7, 0, 0, 0, 0, 0), &
           tMemory('TRI', 'Triad to Column conversion Storage', .false.,  8, 0, 0, 0, 0, 0), &
           tMemory('TPL', '                   Tecplot Storage', .false.,  9, 0, 0, 0, 0, 0), &
           tMemory('TOT', '     Total Allocated SUTRA Storage',  .true., 10, 0, 0, 0, 0, 0) /)

      public :: &
        Memory, &
        MemoryIndex, &
        AddMemory, &
        CalculateTotalStorage

      contains
          integer function MemoryIndex(cID)
            character (len=IDLen), intent(in) :: cID
            !local
            integer (I4B) :: i
            lok=.false.
            MemoryIndex=0
            do i=1,nEntries
              if(index(cID,Memory(i)%ID)>0) then
                MemoryIndex=i
                exit
              end if
            end do
            if(MemoryIndex==0)  &
              call ErrorIO('TotalStorage::MemoryIndex Error could not match '//trim(adjustl(cID))//' in Memory structure')
            lok=.true.
09999&
            return
          end function MemoryIndex

          integer function AddMemory(Index,cType,Size)
            implicit none
            integer (I4B), intent(in)         :: Index
            character (len=IDLen), intent(in) :: cType
            integer (I4B), intent(in)         :: Size

            AddMemory=0
            Memory(Index)%IsUsed=.true.
            select case (cType)
              case('Vec')
                Memory(Index)%Vec=Memory(Index)%Vec+Size*SizeDP
              case('Arr')
                Memory(Index)%Arr=Memory(Index)%Arr+Size*SizeDP
              case('Int')
                Memory(Index)%Int=Memory(Index)%Int+Size*SizeInt
              case('Irr')
                Memory(Index)%Irr=Memory(Index)%Irr+Size*SizeInt
              case('Chr')
                Memory(Index)%Chr=Memory(Index)%Chr+Size*SizeChar
              case default
                AddMemory=1
                call ErrorIO('TotalStorage::AddMemory Error could not match '//trim(adjustl(cType))//'&
							 as a valid type in Memory structure')			!MT: split line due to error when compiling with gfortran 
            end select
            return
          end function AddMemory

          !Calculate and write total memory utilized
          logical function CalculateTotalStorage()
            use FUNITS
            integer (I4B) :: n
          
            lok=.false.

            !Calculate Total Storage
            Memory(MemoryIndex('TOT'))%Vec = sum(Memory(1:(nEntries-1))%Vec)
            Memory(MemoryIndex('TOT'))%Arr = sum(Memory(1:(nEntries-1))%Arr)
            Memory(MemoryIndex('TOT'))%Int = sum(Memory(1:(nEntries-1))%Int)
            Memory(MemoryIndex('TOT'))%Irr = sum(Memory(1:(nEntries-1))%Irr)
            Memory(MemoryIndex('TOT'))%Chr = sum(Memory(1:(nEntries-1))%Chr)

            !write summary information to fLST
            do n=1,nEntries
              if(.not. Memory(n)%IsUsed) cycle
              write (fLST,100) &
                  trim(adjustl(Memory(n)%Title)), &
                  dble(Memory(n)%Vec)/1d6, &
                  dble(Memory(n)%Arr)/1d6, &
                  dble(Memory(n)%Int)/1d6, &
                  dble(Memory(n)%Irr)/1d6, &
                  dble(Memory(n)%Chr)/1d6, &
                  dble(Memory(n)%Vec+Memory(n)%Arr+Memory(n)%Int+Memory(n)%Irr+Memory(n)%Chr)/1d6
            end do

0100 &
            FORMAT(//13X,A/&
                     13X, F10.4,1X,'MBytes OF REAL    VECTORS WERE ALLOCATED ',/, &
                     13X, F10.4,1X,'MBytes of REAL     ARRAYS WERE ALLOCATED ',/, &
                     13X, F10.4,1X,'MBytes OF INTEGER VECTORS WERE ALLOCATED ',/, &
                     13X, F10.4,1X,'MBytes OF INTEGER  ARRAYS WERE ALLOCATED ',/, &
                     13X, F10.4,1X,'MBytes OF CHARACTERS      WERE ALLOCATED ',/, &
                     13X, F10.4,1X,'Mbytes TOTAL MEMORY WAS USED',//)                          

            lok=.true.
  09999&            
            CalculateTotalStorage=lok
            return
          end function CalculateTotalStorage

    end module TotalStorage
