!     SUBROUTINE        P  T  R  S  E  T       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO SET UP POINTER ARRAYS NEEDED TO SPECIFY THE MATRIX STRUCTURE. 
!
! ***  SUBROUTINE ALSO CALLS TRIAD2COLUMN FUNCTIONS TO GENERATE
! ***  THE MAPPING ARRAY FOR CONVERTING A TRIAD STORAGE METHOD
! ***  TO A COLUMN STORAGE METHOD.
!                                                                       

      SUBROUTINE PTRSET (NBI27, ITRI, JTRI, MIOFF) 

      USE CONTRL 
      USE SOLVI 
      USE DIMS
      USE DIMX
      USE FUNITS
      use SutraMSPrecision
      USE MSErrorHandler
      USE Triad2Column

      IMPLICIT NONE
      
      LOGICAL :: lrc
      INTEGER (I4B) :: ival
      INTEGER (I4B) :: ISYM
      integer (I4B) :: &
        NBI27 (NBIX), MIOFF (27), ITRI (NELT), JTRI (NELT) 

      !LOCAL VARIABLES
      INTEGER (I4B) :: &
        IS, IBEG, IEND, ILEN, ILOC, INDOFF, INDLOC, &
        JB, JS, JLOC, &
        KLOC, KS, &
        MBEG, &
        NN12, NB1, NBM, NST, NBH1, NBMC, NBMJ, NBMK, NSTK, NSTJ, NNNBH
!                                                                       
!.....CREATE THE POINTER ARRAY "NBI27", WHICH GIVES THE CORRESPONDENCE  
!        BETWEEN COLUMN INDICES IN SUTRA BANDED FORMAT AND "27-WIDE"    
!        FORMAT.  (FOR 2-D PROBLEMS, IT'S "9-WIDE".)                    
!                                                                       
      IF (IABS (KTYPE) .EQ.3) THEN 
!.....3-D PROBLEM.                                                      
!                                                                       
         DO 200 ILOC = - 1, 1 
            DO 200 JLOC = - 1, 1 
               DO 200 KLOC = - 1, 1 
!                                                                       
!........COMPUTE COLUMN INDEX OFFSET IN SUTRA BANDED FORMAT.            
!           THIS IS EQUIVALENT TO THE FOLLOWING:                        
!           I1 = I + ILOC ;  J1 = J + JLOC ;  fINP = K + KLOC ;           
!           N = NN1*(NN2*(K - 1) + J - 1) + I ;                         
!           N1 = NN1*(NN2*(fINP - 1) + J1 - 1) + I1 ;                     
!           INDOFF = N - N1 .                                           
!                                                                       
                  INDOFF = NN1 * (NN2 * KLOC + JLOC) + ILOC 
!                                                                       
!........COMPUTE LOCAL NODE NUMBER IN 27-NODE "MOLECULE".               
!           THIS IS EQUIVALENT TO THE FOLLOWING:                        
!           I = ILOC + 2 ;  J = JLOC + 2 ;  K = KLOC + 2 ;              
!           INDLOC = NI*(NJ*(K - 1) + J - 1) + I  (NI=NJ=3).            
!                                                                       
                  INDLOC = 9 * KLOC + 3 * JLOC + ILOC + 14 
!                                                                       
!........COMPUTE AND STORE CORRESPONDENCE BETWEEN COLUMN INDEX IN SUTRA 
!           BANDED FORMAT AND LOCAL NODE NUMBER IN 27-NODE "MOLECULE".  
!                                                                       
                  JB = INDOFF + NBHALF 
                  NBI27 (JB) = INDLOC 
!                                                                       
  200    CONTINUE 
!                                                                       
!.....DEFINE CERTAIN QUANTITIES FOR CONVENIENCE AND EFFICIENCY.         
!                                                                       
         NN12 = NN1 * NN2 
         NNNBH = NN - NBHALF 
         NBH1 = NBHALF + 1 
         NB1 = NB + 1 
!                                                                       
!.....CREATE THE POINTER ARRAY "MIOFF", WHICH IS USED AS FOLLOWS TO     
!        COMPUTE THE POSITION, M, OF A MATRIX COEFFICIENT IN THE ARRAY: 
!        M = MIOFF(J27) + I, WHERE "I" IS THE ROW INDEX AND "J27" IS    
!        THE COLUMN INDEX IN "27-WIDE" FORMAT (I.E., THE LOCAL NODE     
!        NUMBER IN THE 27-NODE "MOLECULE").  THIS IS USED IN THE GLOBAL 
!        MATRIX ASSEMBLY ROUTINE "GLOTRI".                              
!                                                                       
         MBEG = 1 
         DO 400 KS = 0, 2 
            NBMK = KS * NN12 
            NSTK = KS * 9 
            DO 400 JS = 0, 2 
               NBMJ = NBMK + JS * NN1 
               NSTJ = NSTK + JS * 3 
               DO 400 IS = 1, 3 
                  NBM = NBMJ + IS 
                  NBMC = NB1 - NBM 
                  NST = NSTJ + IS 
                  IF (NST.LT.14) THEN 
                     IBEG = NBH1 - NBM 
                     IEND = NN 
                  ELSE 
                     IBEG = 1 
                     IEND = NNNBH + NBMC 
                  ENDIF 
                  MIOFF (NST) = MBEG - IBEG 
                  ILEN = IEND-IBEG + 1 
                  MBEG = MBEG + ILEN 
  400    CONTINUE 
!                                                                       
      ELSE 
!.....2-D PROBLEM.                                                      
!                                                                       
         DO 1200 ILOC = - 1, 1 
            DO 1200 JLOC = - 1, 1 
!                                                                       
!........COMPUTE COLUMN INDEX OFFSET IN SUTRA BANDED FORMAT.            
!           THIS IS EQUIVALENT TO THE FOLLOWING:                        
!           I1 = I + ILOC ;  J1 = J + JLOC ;                            
!           N = NN1*(J - 1) + I ;                                       
!           N1 = NN1*(J1 - 1) + I1 ;                                    
!           INDOFF = N - N1 .                                           
!                                                                       
               INDOFF = NN1 * JLOC + ILOC 
!                                                                       
!........COMPUTE LOCAL NODE NUMBER IN 9-NODE "MOLECULE".                
!           THIS IS EQUIVALENT TO THE FOLLOWING:                        
!           I = ILOC + 2 ;  J = JLOC + 2 ;                              
!           INDLOC = NI*(J - 1) + I  (NI=3).                            
!                                                                       
               INDLOC = 3 * JLOC + ILOC + 5 
!                                                                       
!........COMPUTE AND STORE CORRESPONDENCE BETWEEN COLUMN INDEX IN SUTRA 
!           BANDED FORMAT AND LOCAL NODE NUMBER IN 27-NODE "MOLECULE".  
!                                                                       
               JB = INDOFF + NBHALF 
               NBI27 (JB) = INDLOC 
!                                                                       
 1200    CONTINUE 
!                                                                       
!.....DEFINE CERTAIN QUANTITIES FOR CONVENIENCE AND EFFICIENCY.         
!                                                                       
         NNNBH = NN - NBHALF 
         NBH1 = NBHALF + 1 
         NB1 = NB + 1 
!                                                                       
!.....CREATE THE POINTER ARRAY "MIOFF", WHICH IS USED AS FOLLOWS TO     
!        COMPUTE THE POSITION, M, OF A MATRIX COEFFICIENT IN THE ARRAY: 
!        M = MIOFF(J9) + I, WHERE "I" IS THE ROW INDEX AND "J9" IS      
!        THE COLUMN INDEX IN "9-WIDE" FORMAT (I.E., THE LOCAL NODE      
!        NUMBER IN THE 9-NODE "MOLECULE").  THIS IS USED IN THE GLOBAL  
!        MATRIX ASSEMBLY ROUTINE "GLOTRI".                              
!                                                                       
         MBEG = 1 
         DO 1400 JS = 0, 2 
            NBMJ = JS * NN1 
            NSTJ = JS * 3 
            DO 1400 IS = 1, 3 
               NBM = NBMJ + IS 
               NBMC = NB1 - NBM 
               NST = NSTJ + IS 
               IF (NST.LT.5) THEN 
                  IBEG = NBH1 - NBM 
                  IEND = NN 
               ELSE 
                  IBEG = 1 
                  IEND = NNNBH + NBMC 
               ENDIF 
               MIOFF (NST) = MBEG - IBEG 
               ILEN = IEND-IBEG + 1 
               MBEG = MBEG + ILEN 
 1400    CONTINUE 
!                                                                       
      ENDIF 
!                                                                       
!.....CREATE THE POINTER ARRAYS "ITRI" AND "JTRI".                      
!                                                                       
      CALL TRISET (ITRI, JTRI) 
!
!.....PUSH ITRI AND JTRI INTO nIA AND nJA SO IT DOES NOT HAVE TO BE
!     REGENERATED EACH ITERATION FOR SLAP SOLVER
!     SOLVI
      nIA=ITRI
      nJA=JTRI

      do ival=1,NELT
        nAMAP(ival)=ival
      end do
      rAMAP=dble(nAMAP)

!
!.....CONVERT SLAP TRIAD FORMAT TO SLAP COLUMN FORMAT
      ISYM=0
      MSErrorValue%cDataSet='PTR'
      lrc = DS2Yws(NN, NELT, nIA, nJA, rAMAP, ISYM)
      if (.not. lrc) &
        call ErrorIO('PTRSET:: Error converting from SLAP Triad format to SLAP Column format')
      !SUCCESSFUL CONVERSION TO COLUMN FORMAT
      write(*,*) 'successful conversion to SLAP column format'
      !
      !convert rAMAP to nAMAP
      nAMAP=int(rAMAP)
!!
      if (lDebugData) then
        write (fLST,'(6a)') &
          '     Entry','     TriIA','     TriJA','      MapA','     ColIA','     ColJA'
        do ival=1,NELT
          write(fLST,'(6(i10))') &
            ival,ITRI(ival),JTRI(ival),nAMAP(ival),nIA(ival),nJA(ival)
        end do
      end if

!                                                                       
      RETURN 
      END SUBROUTINE PTRSET                         
