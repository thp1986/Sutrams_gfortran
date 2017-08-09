!     MODULE            SUTRAZONEMODULE        SUTRA-MS VERSION 2004.1
                                                                       
! ***
! *** PURPOSE :                                                         
! ***  MODULE THAT IS USED TO PROCESS AQUIFER PARAMETER ZONE DATA
! ***  FOR SIMULATIONS USING THE ZONE (ZON) OPTION.  THE MODULE
! ***  CONTAINS STRUCTURES THAT ALLOCATE STORAGE AND READ DATA FOR
! ***  THE AQUIFER ZONES (RDZONEDATA) AND CALCULATES THE SPECIFIC
! ***  STORAGE VALUES FOR EACH NODE FROM THE SPECIFIED MATRIX AND
! ***  FLUID COMPRESSIBILITIES (MKZONESOP).
! ***

module SutraZoneModule

  use DIMS
  use FUNITS
  USE CONTRL
  USE MSErrorHandler

  implicit none

  integer (I4B) :: &
    NodeZones, &
    ElemZones

  integer (I4B), allocatable :: &
    NodeMap(:), &
    ElemMap(:)

  type tNodeData
    real (DP) :: &
      por,    &  !1
      compma, &  !2
      sop,    &  !3 calculated
      rhos       !4
  end type tNodeData
  
  !VERSION 1.1 SPATIAL DISTRIBUTION OF ADSORPTION AND PRODUCTION
  type tSorpProd
    real (DP), dimension(:), allocatable :: &
      prodf0, &  !1 VERSION 1.1
      prods0, &  !2 VERSION 1.1 
      prodf1, &  !3 VERSION 1.1 
      prods1, &  !4 VERSION 1.1 
      chi1,   &  !5 VERSION 1.1
      chi2       !6 VERSION 1.1
  end type tSorpProd

  type tElemData
    real (DP) :: &
      pmax, &     ! 1
      pmid, &     ! 2
      pmin, &     ! 3
      angle1, &   ! 4
      angle2, &   ! 5
      angle3, &   ! 6 
      anglex, &   ! 7
      pangl1, &   ! 8
      pangl2, &   ! 9
      pangl3, &   !10
      almax, &    !11
      almid, &    !12
      almin, &    !13
      at1max, &   !14
      at1mid, &   !15
      at1min, &   !16
      at2max, &   !17
      at2mid, &   !18
      at2min, &   !19
      permxx, &   !20 calculated
      permxy, &   !21 calculated
      permxz, &   !22 calculated
      permyx, &   !23 calculated
      permyy, &   !24 calculated
      permyz, &   !25 calculated
      permzx, &   !26 calculated
      permzy, &   !27 calculated
      permzz, &   !28 calculated
      lambdas     !29
  end type tElemData

  type (tNodeData), allocatable :: NodeData(:)
  type (tSorpProd), allocatable :: ProdSorp(:)  !VERSION 1.1
  type (tElemData), allocatable :: ElemData(:)

  public :: &
    NodeData, &
    ProdSorp, &
    ElemData, &
    RdZoneData, &
    MkZoneSOP

  contains

  logical function RdZoneData()
    use TotalStorage

!    logical :: LOk=.false.				!MT: Commented out due to error occur during compiling by gfortran (conflict with "lok" in mserrorhandler)
    integer (I4B) :: &
      NLSKIP
    integer (I4B) :: &
      i, &
      k, &
      n, &
      l
    integer (I4B) :: &
      nzn, &
      nze
    integer (I4B) :: &
      MSTRUC
    real (dp) :: &
      D2R, &
      ZERO, &
      radiax, &
      sina, &
      cosa, & 
      sina2, &
      cosa2 

    !if zone file is not defined
    if(fZON==0) then
      NodeZones=nn
      ElemZones=ne
      goto 1010
    end if

    MSErrorValue%cDataSet='ZH1'
    CALL SKPCOM (fZON, NLSKIP) 
    READ (fZON,*,iostat=ios) NodeZones, ElemZones 
    if(ios/=0) call ErrorIO('Error in SutraZoneModule::RdZoneData specifying NodeZones and ElemZones')

    !allocate storage
01010&
    allocate(NodeData(NodeZones), &
             ElemData(ElemZones), &
             NodeMap(nn), &
             ElemMap(ne), &
             stat=ios)
    if(ios/=0) call ErrorIO('SutraZoneModule::RdZoneData could not allocate arrays for hydraulic zones')
    !allocate storage for 
    allocate(ProdSorp(NodeZones),stat=ios)
    if(ios/=0) call ErrorIO('SutraZoneModule::RdZoneData could not allocate arrays for sorption and production parameter zones')
    do n = 1, NodeZones
      allocate(ProdSorp(n)%prodf0(NSPE), &
               ProdSorp(n)%prods0(NSPE), &
               ProdSorp(n)%prodf1(NSPE), &
               ProdSorp(n)%prods1(NSPE), &
               ProdSorp(n)%chi1(NSPE), &
               ProdSorp(n)%chi2(NSPE), &
               stat=ios)
      if(ios/=0) call ErrorIO(&
'SutraZoneModule::RdZoneData could not allocate species variables for sorption and production parameter zones')
    end do
    !Calculate memory requirements
    ios=AddMemory(MemoryIndex('ZON'),'Vec',      4*NodeZones) !NodeData(NodeZones)
    ios=AddMemory(MemoryIndex('ZON'),'Vec',     29*ElemZones) !ElemData(ElemZones), &
    ios=AddMemory(MemoryIndex('ZON'),'Int',               NN) !NodeMap(nn), &
    ios=AddMemory(MemoryIndex('ZON'),'Int',               NE) !ElemMap(ne)
    ios=AddMemory(MemoryIndex('ZON'),'Vec', 6*NSPE*NodeZones) !ProdSorp(NodeZones)

    if(fZON==0) then
      write(fLST,'(//,132("-"),/,1x,a,/,2(1x,a,1x,i6,1x,a,/),132("-"),/)') &
        'SUTRA-MS ZONE FILE is not specified',&
        'Number of   nodal zones = number of    nodes [NN = ',nn,']',&
        'Number of element zones = number of elements [NE = ',ne,']'
      goto 9998
    end if

    !write header for nodal zones
    write(fLST,'(//,1x,a,/,132("-"),/,1x,a10,1x,3(a15,1x))') &
      'NODAL ZONES',&
      '      Zone','       Porosity','Matrix Compres.',' Matrix Density'
    MSErrorValue%cDataSet='ZN1'
    CALL SKPCOM (fZON, NLSKIP) 
    nzn=0
01020&
    read (fZON,*,iostat=ios) i
    if (i==0) goto 1030
    if (i<0) goto 1020
    if (i>NodeZones) goto 1020
    nzn=nzn+1
    backspace(fZON)
    read (fZON,*,iostat=ios) i, &
      NodeData(i)%por,NodeData(i)%compma,NodeData(i)%rhos

    if(ios/=0) call ErrorIO('Error in SutraZoneModule::RdZoneData specifying nodal zone values for zone entry ['//trim(Val2Char(i))&
 //']')

    !write out zone data
    write(fLST,'(1x,i10,1x,3(1PD15.5,1x))') &
      i, &
      NodeData(i)%por,NodeData(i)%compma,NodeData(i)%rhos
    goto 1020
    !check that a nodal value has been defined for each nodal zone
01030 continue
    if(nzn/=NodeZones) &
      call ErrorIO('Error in SutraZoneModule::RdZoneData  - ['//trim(Val2Char(NodeZones))//'] node zones defined but only &
['//trim(Val2Char(nzn))//'] zones found in zone file')

! READ DATA FOR SORPTION AND PRODUCTION VALUES FOR EACH SPECIES      
    !write header for nodal zones
    write(fLST,'(//,1x,a,/,132("-"),/,1x,2(a10,1x),6(a15,1x))') &
      'NODAL ZONES',&
      '      Zone','   Species', &
      '           CHI1','           CHI2', &
      '         PRODF0','         PRODS0', & 
      '         PRODF1','         PRODS1' 
    MSErrorValue%cDataSet='ZN2'
    CALL SKPCOM (fZON, NLSKIP) 
    nzn=0
01040&
    read (fZON,*,iostat=ios) i
    if (i==0) goto 1050
    if (i<0) goto 1040
    if (i>NodeZones) goto 1040
    nzn=nzn+1
    backspace(fZON)
    read (fZON,*,iostat=ios) i, k, &
      ProdSorp(i)%chi1(k),   ProdSorp(i)%chi2(k),   &
      ProdSorp(i)%prodf0(k), ProdSorp(i)%prods0(k), &
      ProdSorp(i)%prodf1(k), ProdSorp(i)%prods1(k)

    if(ios/=0) call ErrorIO('Error in SutraZoneModule::RdZoneData specifying nodal zone values for zone entry &
['//trim(Val2Char(i))//']')

    !write out zone data
    write(fLST,'(1x,2(i10,1x),6(1PD15.5,1x))') &
      i, k, &
      ProdSorp(i)%chi1(k),   ProdSorp(i)%chi2(k),   &
      ProdSorp(i)%prodf0(k), ProdSorp(i)%prods0(k), &
      ProdSorp(i)%prodf1(k), ProdSorp(i)%prods1(k)
    goto 1040
01050 continue
    if(nzn/=NSPE*NodeZones) &
      call ErrorIO('Error in SutraZoneModule::RdZoneData  - &
['//trim(Val2Char(NSPE*NodeZones))//'] node zones defined for ProdSorp but only &
['//trim(Val2Char(nzn))//'] zones found in zone file')

    MSErrorValue%cDataSet='ZE1'
    CALL SKPCOM (fZON, NLSKIP) 
    nze=0
01120&
    read (fZON,*,iostat=ios) i
    if (i==0) goto 1130
    if (i<0) goto 1120
    if (i>ElemZones) goto 1120
    nze=nze+1
    backspace(fZON)
    !3D
    if(iabs(KTYPE)==3) then
      read (fZON,*,iostat=ios) i, &
        ElemData(i)%pmax,ElemData(i)%pmid,ElemData(i)%pmin,&
        ElemData(i)%angle1,ElemData(i)%angle2,ElemData(i)%angle3,&
        ElemData(i)%almax,ElemData(i)%almid,ElemData(i)%almin,&
        ElemData(i)%at1max,ElemData(i)%at1mid,ElemData(i)%at1min,&
        ElemData(i)%lambdas
      if(ios/=0) call ErrorIO('Error in SutraZoneModule::RdZoneData specifying element zone values for zone entry &
['//trim(Val2Char(i))//']')
!
!.....Set AT2 = AT1 - Preserve second transverse dispersivity to allow return to original
!     dispersivity model of Alden Provost
      ElemData(i)%at2max=ElemData(i)%at1max
      ElemData(i)%at2mid=ElemData(i)%at1mid
      ElemData(i)%at2min=ElemData(i)%at1min
 
!                                                                       
!.....ROTATE PERMEABILITY FROM MAX/MID/MIN TO X/Y/Z DIRECTIONS.         
!.....THIS SECTION OF CODE (THROUGH THE CALL TO "TENSOR") WAS LIFTED    
!.....FROM DAVE POLLOCK'S "indat13-dwp.f".                              
      D2R = 1.745329252D-2 
      ElemData(i)%pangl1 = D2R * ElemData(i)%angle1 
      ElemData(i)%pangl2 = D2R * ElemData(i)%angle2 
      ElemData(i)%pangl3 = D2R * ElemData(i)%angle3 
      ZERO = 0D0 
      MSTRUC = 1 
      CALL TENSOR (ElemData(i)%pmax, ZERO, ZERO, ZERO, &
            ElemData(i)%pmid, ZERO, ZERO, ZERO, ElemData(i)%pmin, &
            ElemData(i)%pangl1, ElemData(i)%pangl2, ElemData(i)%pangl3,&
            ElemData(i)%permxx, ElemData(i)%permxy, ElemData(i)%permxz,&
            ElemData(i)%permyx, ElemData(i)%permyy, ElemData(i)%permyz, &
            ElemData(i)%permzx, ElemData(i)%permzy, ElemData(i)%permzz, &
            MSTRUC)                 

      !write element zone data to lst file
      write(fLST,'(//,1x,a,i10,/,132("-"),/,1x,6(a15,1x))') &
        'ELEMENT ZONE',i,&
        'Maximum Permea.',' Middle Permea.','Minimum Permea.',&
        '        Angle 1','        Angle 2','        Angle 3'
      write(fLST,'(1x,6(1PD15.5,1x))') &
        ElemData(i)%pmax, ElemData(i)%pmid, ElemData(i)%pmin, &
        ElemData(i)%angle1, ElemData(i)%angle2, ElemData(i)%angle3

      write(fLST,'(/1x,6(a15,1x))') &
        'Long. Disp. Max','Long. Disp. Mid','Long. Disp. Min',&
        'Tran. Disp. Max','Tran. Disp. Mid','Tran. Disp. Min'
      write(fLST,'(1x,6(1PD15.5,1x))') &
        ElemData(i)%almax,ElemData(i)%almid,ElemData(i)%almin,&
        ElemData(i)%at1max,ElemData(i)%at1mid,ElemData(i)%at1min

      write(fLST,'(/1x,1(a15,1x))') &
        'Matrix T. Cond.'
      write(fLST,'(1x,6(1PD15.5,1x))') &
        ElemData(i)%lambdas

      !write 3x3 pemeability matrix
      write(fLST,'(/,3(1xa,3(a5),a,10x,a,1x,3(1PD15.5,1x),a,/)/,132("-"),/)') &
        '|','xx','xy','xz','|','|',ElemData(i)%permxx,ElemData(i)%permxy,ElemData(i)%permxz,'|',&
        '|','yx','yy','yz','|','|',ElemData(i)%permyx,ElemData(i)%permyy,ElemData(i)%permyz,'|',&
        '|','zx','zy','zz','|','|',ElemData(i)%permzx,ElemData(i)%permzy,ElemData(i)%permzz,'|'

    !2D
    else
      read (fZON,*,iostat=ios) i, &
        ElemData(i)%pmax,ElemData(i)%pmin,&
        ElemData(i)%anglex,&
        ElemData(i)%almax,ElemData(i)%almin,&
        ElemData(i)%at1max,ElemData(i)%at1min,&
        ElemData(i)%lambdas
      if(ios/=0) call ErrorIO('Error in SutraZoneModule::RdZoneData specifying element zone values for zone entry &
['//trim(Val2Char(i))//']')
!                                                                       
!.....ROTATE PERMEABILITY FROM MAXIMUM/MINIMUM TO X/Y DIRECTIONS        
      radiax = 1.745329D-2 * ElemData(i)%anglex
      sina = sin(radiax) 
      cosa = cos(radiax) 
      sina2 = sina * sina 
      cosa2 = cosa * cosa 
      ElemData(i)%permxx = ElemData(i)%pmax * cosa2 + ElemData(i)%pmin * sina2 
      ElemData(i)%permyy = ElemData(i)%pmax * sina2 + ElemData(i)%pmin * cosa2 
      ElemData(i)%permxy =(ElemData(i)%pmax - ElemData(i)%pmin) * sina * cosa 
      ElemData(i)%permyx = ElemData(i)%permxy
      ElemData(i)%pangl1 = radiax 


      !write element zone data to lst file
      write(fLST,'(//,1x,a,i10,/,132("-"),/,1x,3(a15,1x))') &
        'ELEMENT ZONE',i,&
        'Maximum Permea.','Minimum Permea.',&
        '        Angle 1'
      write(fLST,'(1x,3(1PD15.5,1x))') &
        ElemData(i)%pmax, ElemData(i)%pmin, ElemData(i)%anglex

      write(fLST,'(/1x,4(a15,1x))') &
        'Long. Disp. Max','Long. Disp. Min',&
        'Tran. Disp. Max','Tran. Disp. Min'
      write(fLST,'(1x,6(1PD15.5,1x))') &
        ElemData(i)%almax,ElemData(i)%almin,&
        ElemData(i)%at1max,ElemData(i)%at1min

      write(fLST,'(/1x,1(a15,1x))') &
        'Matrix T. Cond.'
      write(fLST,'(1x,6(1PD15.5,1x))') &
        ElemData(i)%lambdas

      !write 2x2 pemeability matrix
      write(fLST,'(/,2(1xa,2(a5),a,10x,a,1x,2(1PD15.5,1x),a,/)/,132("-"),/)') &
        '|','xx','xy','|','|',ElemData(i)%permxx,ElemData(i)%permxy,'|',&
        '|','yx','yy','|','|',ElemData(i)%permyx,ElemData(i)%permyy,'|'


    end if

    goto 1120

    !check that a element value has been defined for each nodal zone
01130 continue
    if(nze/=ElemZones) &
      call ErrorIO('Error in SutraZoneModule::RdZoneData  - &
['//trim(Val2Char(ElemZones))//'] element zones defined but only ['//trim(Val2Char(nze))//'] zones found in zone file')




09998&
    LOk=.true.
09999&
    RdZoneData=LOk
    return
  end function RdZoneData


  logical function MkZoneSOP(COMPFL)

    logical :: LOk=.false.
    integer (I4B) :: &
      i
    real (DP) :: &
      COMPFL

    if(fZON==0) goto 9998
    !write SOP header for nodal zones
    write(fLST,'(//,1x,a,/,132("-"),/,1x,a10,1x,4(a15,1x))') &
      'SOP for NODAL ZONES',&
      '      Zone','       Porosity','Matrix Compres.',' Fluid Compres.','SpecificStorage'

    do i=1,NodeZones
      !set sop for zone
      NodeData(i)%sop =(1.D0 - NodeData(i)%por ) * NodeData(i)%compma + NodeData(i)%por * COMPFL 

      !write out zone data
      write(fLST,'(1x,i10,1x,4(1PD15.5,1x))') &
        i, &
        NodeData(i)%por,NodeData(i)%compma,COMPFL,NodeData(i)%sop

    end do
09998&
    LOk=.true.
09999&
    MkZoneSOP=LOk
    return

  end function MkZoneSOP

end module SutraZoneModule