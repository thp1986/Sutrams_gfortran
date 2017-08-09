!     MODULE            S  O  B  I  O          SUTRA-MS VERSION 2004.1
                                                                       
! ***
! *** PURPOSE :                                                         
! ***  MODULE THAT IS USED TO ALLOCATE STORAGE AND READ SPECIFIED 
! ***  OBSERVATION LOCATION DATA (ALLO_RD_SOBDATA), CALCULATE
! ***  OBSERVATION NODE LOCATIONS FROM X, Y, AND Z DATA (CALCOBSNODE), 
! ***  AND CALCULATING THE CLOSET NODE TO THE SPECIFIED LOCATION
! ***  (CALCULATENODE).  THE CALCULATENODE FUNCTION IS GENERIC CAN
! ***  CAN BE USED BY ANY SUBROUTINE OR FUNCTIONS THAT HAS X, Y, AND Z,
! ***  DATA CONTAINED IN A STRUCTURE USING THE TXYZOBSDATA TYPE.
! ***

    module SobIO
      USE FUNITS
      USE DIMS
      USE OBS
      USE CONTRL
      use SutraMSPrecision
      !MS specific modules
      USE MSErrorHandler
      USE SutraStorage
      implicit none

      type tXYZObsData
        CHARACTER (LEN=128) :: &
          COBSNAME
        character (len=15) :: &
          cInMesh
        integer (I4B) :: &
            node
        real (DP) :: &
            x, &
            y, &
            z, &
            dist, &
            distz
      end type tXYZObsData

      type (tXYZObsData), allocatable :: SobIOData(:)

      integer (I4B) :: &
        SobObs, &
        SobCyc

      public :: &
        tXYZObsData, &
        Allo_Rd_SobData, &
        CalcObsNode, &
        CalculateNode

      contains
        logical function Allo_Rd_SobData(NOBS,NOBCYC)
          use TotalStorage
          integer(I4B), intent(out) :: NOBS
          integer(I4B), intent(out) :: NOBCYC
          integer(I4B) :: &
            n, &
            NLSKIP
          character (len=256) :: cLine
          character (len=128) :: cName

          Allo_Rd_SobData=.false.

          MSErrorValue%cDataSet='SB1'
          CALL SKPCOM (fSOB, NLSKIP) 
          read(fSOB,*,IOSTAT=ios) SobObs, SobCyc
          if(ios/=0) call ErrorIO('Error: Reading SobObs from SOB (specified observation locations) file')

          !set NOBS used in Main Routine to SobObs
          NOBS=SobObs

          !set NOBCYC used in Main Routine to SobCyc
          NOBCYC=SobCyc

          allocate ( &
                     SobIOData(SobObs), &
                     stat=ios)
          if(ios/=0) call ErrorIO('SOBIO::Allo_Rd_SobData Error allocating SobIOData and cInMesh')
          if(.not. allocated(COBSNAME)) &
            allocate( &
                       COBSNAME(NOBS), &
                       stat=ios)
          if (ios /= 0) goto 9999
          !Calculate memory requirements
          ios=AddMemory(MemoryIndex('SOB'),'Int',                           SobObs) !SobIOData(SobObs)%node
          ios=AddMemory(MemoryIndex('SOB'),'Vec',                           SobObs) !SobIOData(SobObs)%x
          ios=AddMemory(MemoryIndex('SOB'),'Vec',                           SobObs) !SobIOData(SobObs)%y
          ios=AddMemory(MemoryIndex('SOB'),'Vec',                           SobObs) !SobIOData(SobObs)%z
          ios=AddMemory(MemoryIndex('SOB'),'Vec',                           SobObs) !SobIOData(SobObs)%dist
          ios=AddMemory(MemoryIndex('SOB'),'Vec',                           SobObs) !SobIOData(SobObs)%distz
          ios=AddMemory(MemoryIndex('SOB'),'Chr',len(SobIOData(1)%cInMesh) *SobObs) !SobIOData(SobObs)%cInMesh
          ios=AddMemory(MemoryIndex('SOB'),'Chr',len(SobIOData(1)%COBSNAME)*SobObs) !SobIOData(SobObs)%COBSNAME

          !read x, y, z data
          MSErrorValue%cDataSet='SB2'
          CALL SKPCOM (fSOB, NLSKIP)
          
          do n=1,SobObs
              if (IABS(KTYPE)==3) then
                  read(fSOB,'(a256)',IOSTAT=ios) &
                    cLine
                  read(cLine,*,IOSTAT=ios) &
                    SobIOData(n)%x,SobIOData(n)%y,SobIOData(n)%z,cName
                  if(ios/=0) then
                    read(cLine,*,IOSTAT=ios) &
                      SobIOData(n)%x,SobIOData(n)%y,SobIOData(n)%z
                  else
                    SobIOData(n)%COBSNAME=cName
                    COBSNAME(n)=cName
                  end if
              else
                  read(fSOB,'(a256)',IOSTAT=ios) &
                    cLine
                  read(cLine,*,IOSTAT=ios) &
                    SobIOData(n)%x,SobIOData(n)%y,cName
                  if(ios/=0) then
                    read(cLine,*,IOSTAT=ios) &
                      SobIOData(n)%x,SobIOData(n)%y
                  else
                    SobIOData(n)%COBSNAME=cName
                    COBSNAME(n)=cName
                  end if
                  SobIOData(n)%z=0.0
              end if
              if(ios/=0) call ErrorIO('Error: Reading x,y,z from SOB entry ['//trim(Val2Char(n))//']')
          end do

09998&
          Allo_Rd_SobData=.true.
09999&
          return
        end function Allo_Rd_SobData



        logical function CalcObsNode()
          integer (I4B) :: &
            i

          CalcObsNode=.false.


          do i=1,SobObs
            if(.not.CalculateNode(SobIOData(i))) goto 9999
            !Set observation node
            IOBS(i)=SobIOData(i)%node
          end do

          !Write observation data to lst file
          WRITE (fLST, '(////11X,"O B S E R V A T I O N   N O D E S")') 
          WRITE (fLST, '(//11X,"OBSERVATIONS WILL BE MADE EVERY ",I5," TIME STEPS," &
                        /11X,"AS WELL AS ON THE FIRST AND LAST TIME STEP,"        &
                        /11X,"FOR A TOTAL OF ",I5," TIME STEPS.")') &
                NOBCYC, NTOBSN-2 
          WRITE (fLST, '(//,2x,&
                       "              X",1x,"              Y",1x,"              Z",1x,&
                       "      NODE",1x,"HORIZ. DISTANCE",1x," VERT. DISTANCE","  MESH LOCATION"," Observation Name",/,92("-"))') 
          do i=1,SobObs
            write(fLST,'(2x,3(f15.7,1x),i10,2(1x,f15.7),a15,1x,a)') &
                  SobIOData(i)%x,SobIOData(i)%y,SobIOData(i)%z, &
                  IOBS(i),SobIOData(i)%dist,SobIOData(i)%distz,SobIOData(i)%cInMesh,trim(adjustl(COBSNAME(i)))
            if(index(SobIOData(i)%cInMesh,'Outside')>0) &
              write(*,'(1x,a,1x,i10,1x,a)') 'SOB Location',i,'is located outside of the SUTRA mesh'
          end do

          !clean up
          deallocate(SobIOData)

          CalcObsNode=.true.
09999&
          return
        end function CalcObsNode

        logical function CalculateNode(XYZ2Node)
          type (tXYZObsData), intent(inout) :: XYZ2Node
          logical       :: &
            lpx, &
            lnx, &
            lpy, &
            lny, &
            lpz, &
            lnz
          integer (I4B) :: &
            minnode
          integer (I4B) :: &
            n, &
            i
          real (DP) :: &
            mindist, &
            minxydist, &
            minzdist, &
            dist, &
            distxy, &
            distz
          lOk=.false.

            XYZ2Node%cInMesh='   Outside Mesh'
            mindist  =1.0d12
            minzdist =1.0d12
            minnode=NN+1
            !reset logical variables
            lpx      = .false.
            lnx      = .false.
            lpy      = .false.
            lny      = .false.
            lpz      = .false.
            lnz      = .false.
            !find location relative to mesh and distance 
            do n=1,NN
              if(IABS(KTYPE)==3) then
                !check that at least observation point is within the mesh
                if((X(n)-XYZ2Node%x)>=0) lpx=.true.
                if((X(n)-XYZ2Node%x)<=0) lnx=.true.
                if((Y(n)-XYZ2Node%y)>=0) lpy=.true.
                if((Y(n)-XYZ2Node%y)<=0) lny=.true.
                if((Z(n)-XYZ2Node%z)>=0) lpz=.true.
                if((Z(n)-XYZ2Node%z)<=0) lnz=.true.
                !calculate horizontal and vertical distance
                dist   = sqrt((X(n)-XYZ2Node%x)**2+ &
                              (Y(n)-XYZ2Node%y)**2+ &
                              (Z(n)-XYZ2Node%z)**2)
                distxy = sqrt((X(n)-XYZ2Node%x)**2+ &
                              (Y(n)-XYZ2Node%y)**2)
                distz  = abs(Z(n)-XYZ2Node%z)
              else
                !check that at least observation point is within the mesh
                if((X(n)-XYZ2Node%x)>=0) lpx=.true.
                if((X(n)-XYZ2Node%x)<=0) lnx=.true.
                if((Y(n)-XYZ2Node%y)>=0) lpy=.true.
                if((Y(n)-XYZ2Node%y)<=0) lny=.true.
                distxy = sqrt((X(n)-XYZ2Node%x)**2+ &
                              (Y(n)-XYZ2Node%y)**2)
                distz = 0.0
                dist  = distxy
              end if
              if(dist < mindist) then
                if( distz <= distxy ) then
                  minnode   = n
                  mindist   = dist
                  minxydist = distxy
                  minzdist  = distz
                end if
              end if
            end do
            !check that location is inside mesh
            if(IABS(KTYPE)==3) then
              if(lpx.and.lnx.and.lpy.and.lny.and.lpz.and.lnz) XYZ2Node%cInMesh='    Inside Mesh'
            else
              if(lpx.and.lnx.and.lpy.and.lny) XYZ2Node%cInMesh='    Inside Mesh'
            end if
            XYZ2Node%node =minnode
            XYZ2Node%dist =minxydist
            XYZ2Node%distz=minzdist

          lOk=.true.
09999&
          CalculateNode=lOk
          return
        end function CalculateNode
      
    end module sobIO