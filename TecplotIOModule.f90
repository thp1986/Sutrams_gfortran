    module TecplotIOInterface
      implicit none

      INTERFACE

        INTEGER*4 FUNCTION tecini &
         (Title, &
          Variables, &
          FName, &
          ScratchDir, &
          Debug, &
          VIsDouble)
          !MS$ATTRIBUTES STDCALL :: tecini
          !MS$ATTRIBUTES REFERENCE :: Title,Variables,FName
          !MS$ATTRIBUTES REFERENCE :: ScratchDir,Debug,VIsDouble
            character*(*) Title
            character*(*) Variables
            character*(*) FName
            character*(*) ScratchDir
            INTEGER*4 Debug
            INTEGER*4 VIsDouble
        END FUNCTION tecini

        INTEGER*4 FUNCTION teczne &
         (ZoneTitle, &
          IMx, &
          JMx, &
          KMx, &
          ZFormat, &
          DupList)
          !MS$ATTRIBUTES STDCALL :: teczne
          !MS$ATTRIBUTES REFERENCE :: ZoneTitle,IMx,JMx,KMx
          !MS$ATTRIBUTES REFERENCE :: ZFormat,DupList
          character*(*) ZoneTitle
          INTEGER*4 IMx
          INTEGER*4 JMx
          INTEGER*4 KMx
          character*(*) ZFormat
          character*(*) DupList
        END FUNCTION teczne

        INTEGER*4 FUNCTION tecdat &
         (N, &
          FieldData, &
          IsDouble)
          use SutraMSPrecision
          !MS$ATTRIBUTES STDCALL :: tecdat
          !MS$ATTRIBUTES REFERENCE :: N,FieldData,IsDouble
          INTEGER*4  N
          !set to DP if outputting double-precision binary data
          !set to SP if outputting single-precision binary data
          REAL (SP)  FieldData(*)
          INTEGER*4  IsDouble
        END FUNCTION tecdat

        INTEGER*4 FUNCTION tecnod &
         (NData)
          !MS$ATTRIBUTES STDCALL :: tecnod
          !MS$ATTRIBUTES REFERENCE :: NData
          INTEGER*4  NData(*)
        END FUNCTION tecnod

        INTEGER*4 FUNCTION teclab &
         (S)
          !MS$ATTRIBUTES STDCALL :: teclab
          !MS$ATTRIBUTES REFERENCE :: S
          character*(*) S
        END FUNCTION teclab

        INTEGER*4 FUNCTION tecusr &
         (S)
          !MS$ATTRIBUTES STDCALL :: tecusr
          !MS$ATTRIBUTES REFERENCE :: S
          character*(*) S
        END FUNCTION tecusr

        INTEGER*4 FUNCTION tecend &
	 (S)
          !MS$ATTRIBUTES STDCALL :: tecend
	  !MS$ATTRIBUTES REFERENCE :: S
        END FUNCTION tecend

        INTEGER*4 FUNCTION tecgeo &
         (XPos, &
          YPos, &
          ZPos, &
          PosCoordMode, &
          AttachToZone, &
          Zone, &
          Color, &
          FillColor, &
          IsFilled, &
          GeomType, &
          LinePattern, &
          PatternLength, &
          LineThickness, &
          NumEllipsePts, &
          ArrowheadStyle, &
          ArrowheadAttachment, &
          ArrowheadSize, &
          ArrowheadAngle, &
          Scope, &
          NumSegments, &
          NumSegPts, &
          XGeomData, &
          YGeomData, &
          ZGeomData, &
          mfc)
          !MS$ATTRIBUTES STDCALL :: tecgeo
          !MS$ATTRIBUTES REFERENCE :: XPos,YPos,ZPos,PosCoordMode
          !MS$ATTRIBUTES REFERENCE :: AttachToZone,Zone,Color,FillColor
          !MS$ATTRIBUTES REFERENCE :: IsFilled,GeomType,LinePattern
          !MS$ATTRIBUTES REFERENCE :: PatternLength,LineThickness
          !MS$ATTRIBUTES REFERENCE :: NumEllipsePts,ArrowheadStyle
          !MS$ATTRIBUTES REFERENCE :: ArrowheadAttachment,ArrowheadSize
          !MS$ATTRIBUTES REFERENCE :: ArrowheadAngle,Scope,NumSegments
          !MS$ATTRIBUTES REFERENCE :: NumSegPts,XGeomData,YGeomData
          !MS$ATTRIBUTES REFERENCE :: ZGeomData,mfc
          REAL*8        XPos
          REAL*8        YPos
          REAL*8        ZPos
          INTEGER*4     PosCoordMode
          INTEGER*4     AttachToZone
          INTEGER*4     Zone
          INTEGER*4     Color
          INTEGER*4     FillColor
          INTEGER*4     IsFilled
          INTEGER*4     GeomType
          INTEGER*4     LinePattern
          REAL*8        PatternLength
          REAL*8        LineThickness
          INTEGER*4     NumEllipsePts
          INTEGER*4     ArrowheadStyle
          INTEGER*4     ArrowheadAttachment
          REAL*8        ArrowheadSize
          REAL*8        ArrowheadAngle
          INTEGER*4     Scope
          INTEGER*4     NumSegments
          INTEGER*4     NumSegPts
          REAL*4        XGeomData(*)
          REAL*4        YGeomData(*)
          REAL*4        ZGeomData(*)
          character*(*) mfc
        END FUNCTION tecgeo

        INTEGER*4 FUNCTION tectxt &
         (XPos, &
          YPos, &
          PosCoordMode, &
          AttachToZone, &
          Zone, &
          BFont, &
          FontHeightUnits, &
          FontHeight, &
          BoxType, &
          BoxMargin, &
          BoxLineThickness, &
          BoxColor, &
          BoxFillColor, &
          Angle, &
          Anchor, &
          LineSpacing, &
          TextColor, &
          Scope, &
          Text, &
          mfc)
          !MS$ATTRIBUTES STDCALL :: tectxt
          !MS$ATTRIBUTES REFERENCE :: XPos,YPos,PosCoordMode,AttachToZone
          !MS$ATTRIBUTES REFERENCE :: Zone,BFont,FontHeightUnits
          !MS$ATTRIBUTES REFERENCE :: FontHeight,BoxType,BoxMargin
          !MS$ATTRIBUTES REFERENCE :: BoxLineThickness,BoxColor
          !MS$ATTRIBUTES REFERENCE :: BoxFillColor,Angle,Anchor
          !MS$ATTRIBUTES REFERENCE :: LineSpacing,TextColor,Scope
          !MS$ATTRIBUTES REFERENCE :: Text,mfc
          REAL*8     XPos
          REAL*8     YPos
          INTEGER*4  PosCoordMode
          INTEGER*4  AttachToZone
          INTEGER*4  Zone
          INTEGER*4  BFont
          INTEGER*4  FontHeightUnits
          REAL*8     FontHeight
          INTEGER*4  BoxType
          REAL*8     BoxMargin
          REAL*8     BoxLineThickness
          INTEGER*4  BoxColor
          INTEGER*4  BoxFillColor
          REAL*8     Angle
          INTEGER*4  Anchor
          REAL*8     LineSpacing
          INTEGER*4  TextColor
          INTEGER*4  Scope
          character*(*) Text
          character*(*) mfc
        END FUNCTION tectxt

        INTEGER*4 FUNCTION tecfil &
         (F)
          !MS$ATTRIBUTES STDCALL :: tecfil
          !MS$ATTRIBUTES REFERENCE :: F
          INTEGER*4  F
        END FUNCTION tecfil

      END INTERFACE

    end module TecplotIOInterface
