
!     SUBROUTINE     S  U  T  R  A  _  M  S     SUTRA-MS VERSION 2004.1
! ***
! *** PURPOSE :                                                         
! ***  MAIN PROGRAM FILE FOR SUTRA-MS.
! ***
!  
!
!|_____________________________________________________________________|
!|                                                                     |
!|                                                                     |
!|                  UNITED STATES GEOLOGICAL SURVEY                    |
!|   GROUNDWATER FLOW AND ENERGY OR SOLUTE TRANSPORT SIMULATION MODEL  |
!|                                                                     |
!|                                                                     |
!|                                                                     |
!|                                                                     |
!|                       _______________________                       |
!|                      |                       |                      |
!|                      |   S   U   T   R   A   |                      |
!|                      |         M   S         |                      |
!|                      |_______________________|                      |
!|                                                                     |
!|                                                                     |
!|                 Saturated    Unsaturated    TRAnsport               |
!|                 =            =              ===                     |
!|                                                                     |
!|                             Multi-Species                           |
!|                             =     =                                 |
!|                                                                     |
!|                                                                     |
!|     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     |
!|     *                                                         *     |
!|     * SUTRA: June 2003                                        *     |
!|     * ->saturated and/or unsaturated groundwater flow         *     |
!|     * ->multiple species reactive solute transport            *     |
!|     *    and/or thermal energy transport                      *     |
!|     * ->two-dimensional areal or cross-sectional simulation   *     |
!|     * ->fully three-dimensional simulation                    *     |
!|     * ->either (2-D, 3-D) cartesian or                        *     |
!|     *   (2-D) radial/cylindrical coordinates                  *     |
!|     * ->hybrid galerkin-finite-element method and             *     |
!|     *    integrated-finite-difference method                  *     |
!|     *    with two-dimensional quadrilateral or                *     |
!|     *    three-dimensional hexahedral finite elements         *     |
!|     * ->finite-difference time discretization                 *     |
!|     * ->non-linear iterative, sequential or steady-state      *     |
!|     *    solution modes                                       *     |
!|     * ->direct and iterative solvers                          *     |
!|     * ->optional fluid velocity calculation                   *     |
!|     * ->optional observation well output                      *     |
!|     * ->optional printer plots of output                      *     |
!|     * ->optional fluid mass and solute mass or energy budget  *     |
!|     * ->flexible, columnwise output                           *     |
!|     *                                                         *     |
!|     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     |
!|                                                                     |
!|                                                                     |
!|     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     |
!|     *                                                         *     |
!|     * SUTRA-MS: December 2004, Version 2004.1                 *     |
!|     * ->allows sequential solution of pressure, temperature,  *     |
!|     *    and concentration for more than one species          *     |
!|     * ->allows solution of multiple chemical species.         *     |
!|     *    Each species can affect the fluid density and        *     |
!|     *    viscosity.                                           *     |
!|     * ->allows the bulk thermal conductivity to be calculated *     |
!|     *    using a three phase (air, water, and solid matrix)   *     |
!|     *    geometric mean equation.                             *     |
!|     * ->allows spatially varying solid matrix thermal         *     |
!|     *    conductivities to be applied.                        *     |
!|     * ->allows a different dispersivity multiplier to be      *     |
!|     *    to each simulated species.                           *     |
!|     * ->Time-dependent boundary conditions can be specified   *     |
!|     *    in a separate input data file (sutra.tbc)            *     |
!|     * ->User specified output times can be specified for      *     |
!|     *    nodal and elemental output in a separate input data  *     |
!|     *    file (sutra.otm).                                    *     |
!|     * ->A simple Automatic Time-Stepping algorithm can be     *     |
!|     *    used which adjusts the simulation time-step based    *     |
!|     *    on user specified parameters and the number of       *     |
!|     *    iterations required to achieve convergence.  The     *     |
!|     *    algorithm also has the optional to rerun a time-     *     |
!|     *    step if the number of iterations exceeds user        *     |
!|     *    specified values.  Data for the Automatic Time-      *     |
!|     *    Stepping algorithm is specified in a separate input  *     |
!|     *    data file (sutra.ats).                               *     |
!|     * ->Observations can be specified using x, y, and z       *     |
!|     *    coordinates in a separate input data file            *     |
!|     *   (sutra.sob).                                          *     |
!|     *                                                         *     |
!|     * SUTRA-MS: November 2009, Version 1.1                    *     |
!|     * ->Nodal and Elemental aquifer parameters can be         *     |
!|     *    specified using zones defined in the main SUTRA-MS   *     |
!|     *    input data file (sutra.inp) using the NREG and LREG  *     |
!|     *    parameters in DataSets 14B and 15B and aquifer       *     |
!|     *    parametersx specified in a separate input data file  *     |
!|     *   (sutra.zon).                                          *     |
!|     * ->Optional output of nodal and elemental data to        *     |
!|     *    ASCII Tecplot files.  Specified using TPN and TPE    *     |
!|     *    Keywords in sutra.fil.  Output frequency based on    *     |
!|     *    output control in main file (sutra.inp) or in user   *     |
!|     *    specified output times file (sutra.otm).             *     |
!|     *                                                         *     |
!|     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     |
!|                                                                     |
!|                                                                     |
!|                                                                     |
!|       Complete explanation of function and use of SUTRA-MS          |
!|       is given in :                                                 |
!|                                                                     |
!|       Hughes, Joseph D. and Sanford, Ward E., 2004, SUTRA-MS        |
!|            A Version of SUTRA Modified to Simulate Heat and         |
!|            Multiple-Solute Transport, U.S. Geological Survey        |
!|            Open-File Report 2004-1207.                              |
!|                                                                     |
!|                                                                     |
!|       Complete explanation of function and use of SUTRA             |
!|       is given in :                                                 |
!|                                                                     |
!|       Voss, Clifford I. and Provost, Alden M., 2002, SUTRA          |
!|            A Model for Saturated-Unsaturated, Variable-Density      |
!|            Ground-Water Flow with Solute or Energy Transport,       |
!|            U.S. Geological Survey Water-Resources                   |
!|            Investigations Report 02-4231.                           |
!|                                                                     |
!|                       which supersedes                              |
!|                                                                     |
!|       Voss, Clifford I., 1984, SUTRA: A Finite-Element              |
!|            Simulation Model for Saturated-Unsaturated               |
!|            Fluid-Density-Dependent Ground-Water Flow                |
!|            with Energy Transport or Chemically-Reactive             |
!|            Single-Species Solute Transport, U.S. Geological         |
!|            Survey Water-Resources Investigations Report             |
!|            84-4369.                                                 |
!|                                                                     |
!|                                                                     |
!|                                                                     |
!|       Users who wish to be notified of updates of the SUTRA-MS      |
!|       code and documentation may be added to the mailing            |
!|       by sending a request to :                                     |
!|                                                                     |
!|                     Chief Hydrologist - SUTRA-MS                    |
!|                       U.S. Geological Survey                        |
!|                        431 National Center                          |
!|                       Reston, Virginia 20192                        |
!|                                USA                                  |
!|                                                                     |
!|     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     |
!|     * The SUTRA code and documentation were prepared under a  *     |
!|     * joint research project of the U.S. Geological Survey,   *     |
!|     * Department of the Interior, Reston, Virginia, and the   *     |
!|     * Engineering and Services Laboratory,  U.S. Air Force    *     |
!|     * Engineering and Services Center, Tyndall A.F.B.,        *     |
!|     * Florida.                                                *     |
!|     *                                                         *     |
!|     * The SUTRA-MS code and documentation are available for   *     |
!|     * unlimited distribution.                                 *     |
!|     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     |
!|                                                                     |
!|     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     |
!|     *                  S U T R A  -  M  S                     *     |
!|     *                   REVISION HISTORY                      *     |
!|     *                                                         *     |
!|     * First Revision: November 2009, Version 1.1              *     |
!|     * by: Joseph D. Hughes, U.S. Geological Survey            *     |
!|     *                                                         *     |
!|     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     |
!|                                                                     |
!|     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     |
!|     *                      S U T R A                          *     |
!|     *                   REVISION HISTORY                      *     |
!|     *                                                         *     |
!|     * First Revision: June 1990, Version V06902D              *     |
!|     * by: Clifford I. Voss, U.S. Geological Survey            *     |
!|     *                                                         *     |
!|     * Second Revision: September 1997, Version V09972D        *     |
!|     * by: C.I. Voss and David Boldt, U.S. Geological Survey   *     |
!|     *                                                         *     |
!|     * Third Revision: June 2003, Version 2D3D.1               *     |
!|     * by: Clifford I. Voss and Alden M. Provost,              *     |
!|     * U.S. Geological Survey                                  *     |
!|     *                                                         *     |
!|     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *     |
!|_____________________________________________________________________|
!                                                                       
!                                                                       
!                                                                       
      PROGRAM SUTRA_MS 
      USE ITERAT 
      USE PARAMS 
      USE CONTRL 
      USE SOLVI 
      USE ITSOLR 
      USE ITSOLI 
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE OBS
      USE SOLVN
      USE SOLVC
      USE PLT1
      USE SOURCEITEMS
      !MS specific modules
      use SutraMSPrecision
      USE ATSDATA 
      USE MSErrorHandler
      USE SutraZoneModule
      USE SutraStorage
      USE SobIO
      USE Triad2Column
      USE SolverStorage
      USE TECPLOT
      USE ColumnStorage

      IMPLICIT NONE

      CHARACTER(len=80):: TITLE1, TITLE2 
      CHARACTER(len=80):: SIMULA(4)
      CHARACTER(len=80):: MSHTYP(2)
      CHARACTER(len=80):: SIMSTR, MSHSTR 
      CHARACTER(len=80):: CUNSAT(1), CSSFLO(1), CSSTRA(1), CREAD (1)
      CHARACTER(len=80):: UNSSTR, SSFSTR, SSTSTR, RDSTR 
      INTEGER(I4B)     :: RMDIM, RVDIM, IMVDIM 
      INTEGER(I4B)     :: &
        NNSP, &
        NNIZ, &
        NEIZ, &
        NNZV, &
        NEZV, &
        NEZV2, &
        NEZV3, &
        NNV, &
        NEV, &
        NEV2, &
        NEV3, &
        NEVX, &
        NEVG, &
        NE8, &
        NE8X, &
        NMULTI, &
        IDUM0, &
        IDUM1, &
        IDUM2, &
        NOBS, &
        NTOBS, &
        NN123, &
        NE123, &
        NFILE, &
        NNU, &
        NNS, &
        ISMERR, &
        NLSTOT, &
        NWIP, &
        NWFP, &
        NWIU, &
        NWFU, &
        KMXSOLVU, &
        NMXSAVEU, &
        MATDIM, &
        MATOBS, &
        NSPEDIM, &
        MRMD, &
        MRVD, &
        MIMVD
      INTEGER(I4B)     :: &
        NLSKIP, &
        NWORDS
      INTEGER(I4B)     :: &
        I, &
        K, &
        I0, &
        I1, &
        I2, &
        NF, &
        M, &
        JT, &
        JTMAX, &
        KT

      REAL(DP)         :: &
        RT0, &
        RT1, &
        RTF, &
        TS, &
        DELTK

      REAL(DP), ALLOCATABLE :: &
        TEMPTIMOBS(:)

      
!....."NSLVRS" AND THE ARRAYS "SOLWRD" AND "SOLNAM" ARE INITIALIZED     
!.....IN MODULE.F90 IN THE SOLVC AND SOLVN MODULES
!                                                                       
!                                                                       
!_____________________________________________________________________  
!|* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *|  
!|* *************************************************************** *|  
!|* *                                                             * *|  
!|* *   Three arrays are dimensioned dynamically as follows:      * *|  
!|* *                                                             * *|  
!|* *   ALLOCATE(RM(1:RMDIM), RV(1:RVDIM), IMV(1:IMVDIM))        * *|  
!|* *                                                             * *|  
!|* *   RMDIM = 2*NN*NCBI                                         * *|  
!|* *                                                             * *|  
!|* *   RVDIM = NNV*NN +(NEV+NEVG)*NE + NEVX + NE8X + NBCN*3     * *|  
!|* *              +(NOBS+1)*(NTOBS+2)*2 + NTOBS + NDIMXR + 5    * *|  
!|* *              + NNSP*NN*NSPE                                 * *|  
!|* *                                                             * *|  
!|* *   IMVDIM = NE*9 + NN + NSOP + MSOU*NSPE                     * *|  
!|* *              + NBCN + NBCN*NSPE + NOBS + NTOBS +            * *|  
!|* *              + NDIMXI + 9                                   * *|  
!|* *                                                             * *|  
!|* *   where:                                                    * *|  
!|* *                                                             * *|  
!|* *    NNV   = 27 - NNSP                                        * *|  
!|* *    NEV   = 24 for 3-D; 12 for 2-D                           * *|  
!|* *    NEVG  = 16 for 3-D; 8 for 2-D                            * *|  
!|* *    NEVX  = 24 - NEV                                         * *|  
!|* *    NE8X  = NE*8 for 3-D; 4 for 2-D                          * *|  
!|* *    NBCN  = NPBC + MNUBC*NSPE                                * *|  
!|* *    NNSP  = 7                                                * *|  
!|* *                                                             * *|  
!|* *   and:                                                      * *|  
!|* *                                                             * *|  
!|* *    NN = number of nodes in finite element mesh              * *|  
!|* *    NE = number of elements in finite element mesh           * *|  
!|* *    NBI = estimated full bandwidth in finite element mesh    * *|  
!|* *    NCBI = NBI for Gaussian (direct) solver; =27 otherwise   * *|  
!|* *    NDIMXI = integer array storage associated with solver(s),* *|  
!|* *       exclusive of the matrix arrays (see NELT below);      * *|  
!|* *       if two different solvers are used, the storage        * *|  
!|* *       required is the greater of the two values             * *|  
!|* *       = 6 for Gaussian (direct) solver;                     * *|  
!|* *       = 5*NL-NN+NBI+52 for CG with non-symmetric storage;   * *|  
!|* *       = 3*NELT+3*NN+NBI+73 for GMRES;                       * *|  
!|* *       = 3*NELT+3*NN+NBI+53 for ORTHOMIN                     * *|  
!|* *    NDIMXR = real array storage associated with solver(s),   * *|  
!|* *       exclusive of the matrix arrays (see NELT below);      * *|  
!|* *       if two different solvers are used, the storage        * *|  
!|* *       required is the greater of the two values             * *|  
!|* *       = 4 for Gaussian (direct) solver;                     * *|  
!|* *       = NL+6*NN+41 for CG with non-symmetric storage;       * *|  
!|* *       = NELT+NN*(NSAVE+7)+NSAVE*(NSAVE+3)+42 for GMRES;     * *|  
!|* *       = NELT+NN*(3*NSAVE+10)+NSAVE+41 for ORTHOMIN          * *|  
!|* *    NELT = length of matrix arrays                           * *|  
!|* *       = NN for Gaussian (direct) solver;                    * *|  
!|* *       = 27*NN-6*NN1*(1+3*NN2)-2 for 3-D, using GMRES,       * *|  
!|* *         ORTHOMIN, or CG with non-symmetric storage;         * *|  
!|* *       = 9*NN-6*NN1-2 for 2-D, using GMRES,                  * *|  
!|* *         ORTHOMIN, or CG with non-symmetric storage          * *|  
!|* *    NL = (NELT + NN)/2 for CG with non-symmetric storage     * *|  
!|* *    NN1 = number of nodes in 1st node numbering direction,   * *|  
!|* *       for iterative solver                                  * *|  
!|* *    NN2 = number of nodes in 2nd node numbering direction,   * *|  
!|* *       for iterative solver                                  * *|  
!|* *    NN3 = number of nodes in 3rd node numbering direction,   * *|  
!|* *       for 3-D using iterative solver                        * *|  
!|* *    NSAVE = maximum number of direction vectors for GMRES    * *|  
!|* *       or ORTHOMIN                                           * *|  
!|* *    NOBS = number of observation nodes in mesh               * *|  
!|* *    NTOBS = number of time steps with observations           * *|  
!|* *    NSOP = number of fluid mass source nodes in mesh         * *|  
!|* *    NSOU = number of energy or solute mass source nodes      * *|  
!|* *    NPBC = number of specified pressure nodes in mesh        * *|  
!|* *    NUBC = number of specified concentration or temperature  * *|  
!|* *       nodes in mesh                                         * *|  
!|* *                                                             * *|  
!|* *                                                             * *|  
!|* *    NSPE  = NUMBER OF SPECIES                                * *|  
!|* *    MNSOU = MAXIMUM NUMBER OF NSOU FOR ALL SPECIES           * *|  
!|* *    MNUBC = MAXIMUM NUMBER OF NUBC FOR ALL SPECIES           * *|  
!|* *    NLBC  = NUMBER OF LINEAR TIME-DEPENDENT B.C.'s           * *|  
!|* *                                                             * *|  
!|* *                                                             * *|  
!|* *   NOTE :                                                    * *|  
!|* *    Two files must be permanently assigned just below for    * *|  
!|* *    your computer installation.  One file captures error     * *|  
!|* *    output written during subsequent file opening.  The      * *|  
!|* *    other file contains the unit numbers and file names      * *|  
!|* *    to be assigned as SUTRA input and output files           * *|  
!|* *    for each simulation.                                     * *|  
!|* *                                                             * *|  
!|* *    STANDARD ASSIGNMENTS TO BE MADE:                         * *|  
!|* *    for Error Output:                                        * *|  
!|* *        Filename is contained in ENAME                       * *|  
!|* *        Unit Number is contained in fSMY                     * *|  
!|* *    for Simulation Units and Files:                          * *|  
!|* *        Filename is contained in UNAME                       * *|  
!|* *        Unit Number is contained in fSutraFil                * *|  
!|* *                                                             * *|  
!|* *************************************************************** *|  
!|* *************************************************************** *|  
!                                                                       
!
!                                                                       
!.....SIMULATION TIME                                                   
      CALL CPU_TIME(RT0)
!.....Normal Termination Flag
      LNormal=.false.
!                                                                       
!                                                                       
!|* *****  S T A N D A R D   F I L E   A S S I G N M E N T S  ***** *|  
!|*   E R R O R   O U T P U T   &   L O G                         * *|  
      cSMY = 'SUTRA.SMY' 
      fSMY = 1 
!|*   S I M U L A T I O N   U N I T S   A N D   F I L E S         * *|  
      cSutraFil = 'SUTRA.FIL' 
      fSutraFil = 99 
!|*                                                               * *|  
!|*   -------> Suggested Format of Unit fSutraFil :               * *|  
!|*                                                               * *|  
!|*       V A R I A B L E S                     F O R M A T       * *|  
!|*                                                               * *|  
!|*       Key Word   Unit Number for fINP        (free format)    * *|  
!|*                   File Name for fINP         (A80)            * *|  
!|*       Key Word   Unit Number for fICS        (free format)    * *|  
!|*                   File Name for fICS         (A80)            * *|  
!|*       Key Word   Unit Number for fLST        (free format)    * *|  
!|*                   File Name for fLST         (A80)            * *|  
!|*       Key Word   Unit Number for fRST        (free format)    * *|  
!|*                   File Name for fRST         (A80)            * *|  
!|*       Key Word   Unit Number for fNOD        (free format)    * *|  
!|*                   File Name for fNOD         (A80)            * *|  
!|*       Key Word   Unit Number for fELE        (free format)    * *|  
!|*                   File Name for fELE         (A80)            * *|  
!|*       Key Word   Unit Number for fOBS        (free format)    * *|  
!|*                   File Name for fOBS         (A80)            * *|  
!|*       Key Word   Unit Number for fOTM       (free format)     * *|  
!|*                   File Name for fOTM        (A80)             * *|  
!|*       Key Word   Unit Number for fATS       (free format)     * *|  
!|*                   File Name for fATS        (A80)             * *|  
!|*       Key Word   Unit Number for fSOB       (free format)     * *|  
!|*                   File Name for fSOB        (A80)             * *|  
!|*       Key Word   Unit Number for fZON       (free format)     * *|  
!|*                   File Name for fZON        (A80)             * *|  
!|*                                                               * *|  
!|*        fINP IS THE SUTRA2D3DMS INPUT DATA FILE                * *|  
!|*        fICS IS THE SUTRA2D3DMS INITIAL CONDITIONS DATA FILE   * *|  
!|*        fLST IS THE SUTRA2D3DMS OUTPUT DATA FILE               * *|  
!|*        fRST IS THE SUTRA2D3DMS RESTART DATA FILE              * *|  
!|*        fNOD IS THE SUTRA2D3DMS NODAL OUTPUT DATA FILE         * *|  
!|*        fELE IS THE SUTRA2D3DMS VELOCITY OUTPUT DATA FILE      * *|  
!|*        fOBS IS THE SUTRA2D3DMS OBSERVATION DATA OUTPUT FILE   * *|  
!|*        fTBC IS THE TIME-DEPENDENT B.C. INPUT DATA FILE        * *|  
!|*        fOTM IS THE PRINT OUTPUT CONTROL DATA FILE FOR         * *|  
!|*           UNITS fNOD, AND fELE.                               * *|  
!|*        fATS IS AUTOMATIC TIMESTEP DATA FILE WHICH CONTROLS    * *|  
!|*           REDUCTIONS IN THE TIMESTEP BASED ON ITERATIVE       * *|  
!|*           SOLVER PERFORMANCE (NUMBER OF ITERATIONS).          * *|  
!|*        fSOB IS SPECIFIED OBSERVATION LOCATIONS THAT ARE       * *|  
!|*           INTERNALLY CONVERTED TO THE CLOSEST NODAL           * *|  
!|*           LOCATION.                                           * *|  
!|*        fZON IS HYDRAULIC PARAMETERS FOR NODE WISE AND CELL    * *|  
!|*           WISE PARAMETERS TO BE APPLIED TO ZONES SPECIFIED    * *|  
!|*           BY LREG AND NREG.                                   * *|  
!|*                                                               * *|  
!|*                                                               * *|  
!|*       THE Key Word IS A TWO TO THREE CHARACTER KEYWORD FOR    * *|  
!|*       WHICH INPUT OR OUTPUT DATA FILE IS TO BE ASSIGNED       * *|  
!|*       THE SPECIFIED UNIT NUMBER AND FILE NAME.  THE           * *|  
!|*       KEYWORD CORRESPONDS DIRECTLY TO THE FUNITS MODULE       * *|  
!|*       VARIABLE NAME FOR THE FILE.                             * *|  
!|*                                                               * *|  
!|*       THE Key Word VARIABLE CAN BE HAVE ANY OF THE FOLLOWING  * *|  
!|*       VALUES - K1,K2,K3,K4,K5,K6,K7,                          * *|     
!|*                TBC,OTM,ATS,                                   * *|  
!|*                SOb, or ZON                                    * *|  
!|*       EQUIVALENT Key Words ARE                                * *|  
!|*                INP, ICS, LST, RST, NOD, ELE, OBS,             * *|  
!|*                TBC, OTM, ATS,                                 * *|
!|*                SOB, or ZON                                    * *|  
!|*                                                               * *|  
!|*       FOR EXAMPLE:                                            * *|  
!|*       THE INP INPUT DATA FILE WOULD REQUIRE THE Key Word      * *|  
!|*       "INP" OR "K1"                                           * *|  
!|*       and                                                     * *|  
!|*       THE OBS OUTPUT DATA FILE WOULD REQUIRE THE Key Word     * *|  
!|*       "OBS" OR "K7"                                           * *|  
!|*       etc.                                                    * *|  
!|*                                                               * *|  
!|*   Some of the last nine lines need not be included if         * *|  
!|*   UNIT-RST,NOD,ELE,OBS,TBC,OTM,ATS,or ZON are not used.       * *|  
!|*     This file has between six and twelve lines.               * *|  
!|*                                                               * *|  
!|*   THE DATA BLOCKS FOR EACH INPUT AND OUTPUT DATA FILE CAN     * *|  
!|*   BE SPECIFIED IN ANY ORDER SINCE A KEYWORD (i.e., INP,       * *|  
!|*   ICS,LST,RST,NOD,ELE,OBS,TBC,OTM,ATS,SOB, and/or ZON         * *|  
!|*   IS REQUIRED.                                                * *|  
!|*                                                               * *|  
!|*                                                               * *|  
!|* *************************************************************** *|  
!|* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *|  
!|___________________________________________________________________|  
!                                                                       
!                                                                       
! ---> Programmers making code changes that affect dimensions must      
! --->  check and change the following assignments for NNV and NEV:     
!                                                                       
!      NNIZ IS THE NUMBER OF INTEGER VECTORS THAT ARE NN LONG AND DEFINED AS ZONES.        
      NNIZ = 1
!
!      NEIZ IS THE NUMBER OF INTEGER VECTORS THAT ARE NE LONG AND DEFINED AS ZONES.        
      NEIZ = 1

!      NNSP IS THE NUMBER OF REAL VECTORS THAT ARE NSPE*NN LONG.        
      NNSP = 7 
!
!      NNZV is the number of REAL VECTORS that are defined as zones
!      POR and SOP
      NNZV = 2
!                                                                       
!      NNV IS NUMBER OF REAL VECTORS THAT ARE NN LONG                   
      NNV = 27 - NNSP - NNZV

!      NEZV2 is the number of REAL VECTORS that are defined as zones(2D)
!      POR and SOP
      NEZV2 = 9
!      NEZV3 is the number of REAL VECTORS that are defined as zones(3D)
      NEZV3 = 21

!      NEV IS NUMBER OF REAL VECTORS THAT ARE NE LONG.                  
!         NEV = NEV2 for 2-D;                                           
!         NEV = NEV3 for 3-D                                            
      NEV2 = 12 - NEZV2
      NEV3 = 24 - NEZV3
!                                                                       
!                                                                       
!.....KEEP TRACK IF OUTPUT ROUTINES HAVE BEEN EXECUTED, TO PRINT        
!      HEADERS ONLY ONCE                                                
      ONCEK5 = .FALSE. 
      ONCEK6 = .FALSE. 
!                                                                       
!.....ASSIGN UNIT NUMBERS AND OPEN FILE UNITS FOR THIS SIMULATION       
      write (*,*) 'Opening input files...'
      CALL FOPEN()
!                                                                       
!.....OUTPUT BANNER                                                     
      WRITE (fLST, 110)
  110 FORMAT(1H1,131(1H*)////3(132(1H*)////)////                        &
     &   47X,' SSSS   UU  UU  TTTTTT  RRRRR     AA  '/                  &
     &   47X,'SS   S  UU  UU  T TT T  RR  RR   AAAA '/                  &
     &   47X,'SSSS    UU  UU    TT    RRRRR   AA  AA'/                  &
     &   47X,'    SS  UU  UU    TT    RR R    AAAAAA'/                  &
     &   47X,'SS  SS  UU  UU    TT    RR RR   AA  AA'/                  &
     &   47X,' SSSS    UUUU     TT    RR  RR  AA  AA'/                  &
     &   /,                                                             &
     &   47X,'           MMM  MM   SSSS             '/                  &
     &   47X,'           MMM MMM  SS   S            '/                  &
     &   47X,'           MMMMMMM  SSSS              '/                  &
     &   47X,'           MM M MM      SS            '/                  &
     &   47X,'           MM   MM  SS  SS            '/                  &
     &   47X,'           MM   MM   SSSS             '/                  &
     &   7(/),37X,'U N I T E D    S T A T E S   ',                      &
     &   'G E O L O G I C A L   S U R V E Y'////                        &
     &   45X,'SUBSURFACE FLOW AND TRANSPORT SIMULATION MODEL'/          &
     &   //57X,'-SUTRA-MS VERSION  1.1 11/2009-'///                     &
     &   32X,'*  SATURATED-UNSATURATED FLOW AND SOLUTE AND/OR ENERGY',  &
     &   ' TRANSPORT  *'////4(////132(1H*)))                           
!                                                                       
      write (*,*) 'Initializing input data...'
      write (*,*) '  Reading Data Sets 1-3...'
!.....INPUT DATASET 1:  OUTPUT HEADING                                  
      MSErrorValue%cDataSet='  1'
      CALL SKPCOM(fINP, NLSKIP)
      READ (fINP,'(a80,/,a80)',iostat=ios) TITLE1, TITLE2 
      IF (ios/=0) call ErrorIO('Error specifying title lines 1 and 2')
!                                                                       
!.....INPUT DATASET 2A:  INPUT DATA HEADING -- TYPE OF TRANSPORT        
!.....( SET ME=-1 FOR SOLUTE TRANSPORT, ME=+1 FOR ENERGY TRANSPORT,     
!.....  SET ME= 0 FOR SOLUTE AND ENERGY TRANSPORT                  )   
      MSErrorValue%cDataSet=' 2A'
      CALL SKPCOM(fINP, NLSKIP)
      READ (fINP,*,iostat=ios)SIMSTR 
      IF (ios/=0)call ErrorIO('Error specifying transport type')
      CALL PRSWDS(SIMSTR, 4, SIMULA, NWORDS)
!                                                                       
!                                                                       
!.....SEE IF MULTI SPECIES SIMULATION                                   
!     IF ANY SIMULA IS 'MULTI' THEN IT MODIFIES THE DATASET 3 READ      
!     STATEMENT                                                         
      NMULTI = 0 
      DO K = 2, 4 
      IF (INDEX(SIMULA(K), 'MULTI').GT.0) THEN 
         NMULTI = 1 
         EXIT 
      ENDIF 
      ENDDO 
!                                                                       
      IF (SIMULA(1).NE.'SUTRA     ')call ErrorIO('Error - first word must be "SUTRA"')
!.....TEST FOR VALID KEY WORDS                                          
      DO K = 2, 4 
      IF (SIMULA(K).EQ.'SOLUTE    ')GOTO 120 
      IF (SIMULA(K).EQ.'ENERGY    ')GOTO 140 
      ENDDO 
!                                                                       
  115 call ErrorIO('Error - "SOLUTE" and/or "ENERGY" transport must be specified')

  120 ME = - 1 
      WRITE (fLST, 130)
  130 FORMAT(1H1//132(1H*)///20X,'* * * * *   S U T R A   S O L U ',    &
     &   'T E   T R A N S P O R T   S I M U L A T I O N   * * * * *'//  &
     &   /132(1H*)/)                                                   
!                                                                       
!.....TEST IF ALSO SIMULATING ENERGY TRANSPORT                          
      DO K = 2, 4 
      IF (SIMULA(K).EQ.'ENERGY    ') THEN 
         ME = 0 
         NMULTI = 1 
         WRITE (fLST, 150)
         EXIT 
      ENDIF 
      ENDDO 
!                                                                       
      GOTO 160 
  140 ME = + 1 
      WRITE (fLST, 150)
  150 FORMAT(1H1//132(1H*) &
             ///20X,'* * * * *   S U T R A   E N E R G Y   T R A N S P O R T   S I M U L A T I O N   * * * * *'//  &
             /132(1H*)/)                                                   
!                                                                       
!.....TEST IF ALSO SIMULATING SOLUTE TRANSPORT                          
      DO K = 2, 4 
      IF (SIMULA(K).EQ.'SOLUTE    ') THEN 
         ME = 0 
         NMULTI = 1 
         WRITE (fLST, 130)
         EXIT 
      ENDIF 
      ENDDO 
!                                                                       
  160 CONTINUE 
!                                                                       
!.....INPUT DATASET 2B:  INPUT DATA HEADING -- MESH TYPE                
      MSErrorValue%cDataSet=' 2B'
      CALL SKPCOM(fINP, NLSKIP)
      READ (fINP,*)MSHSTR 
      IF (ios/=0)call ErrorIO('Error specifying mesh type')
      CALL PRSWDS(MSHSTR, 2, MSHTYP, NWORDS)
!.....KTYPE SET ACCORDING TO THE TYPE OF FINITE-ELEMENT MESH:           
!        3-D, IRREGULAR MESH   ==>   KTYPE = +3                         
!        3-D, REGULAR MESH     ==>   KTYPE = -3                         
!        2-D, IRREGULAR MESH   ==>   KTYPE = +2                         
!        2-D, REGULAR MESH     ==>   KTYPE = -2                         
      IF (MSHTYP(1).EQ.'2D        ') THEN 
         KTYPE = 2 
      ELSEIF (MSHTYP(1).EQ.'3D        ') THEN 
         KTYPE = 3 
      ELSE 
         call ErrorIO('Error: mesh type must be 2D or 3D')
      ENDIF 
      IF (MSHTYP(2).EQ.'IRREGULAR ') THEN 
!         IF (KTYPE.EQ.3) THEN 
!           call ErrorIO('Error: 3D IRREGULAR MESH is currently not allowed')
!         ENDIF 
      ELSEIF ((MSHTYP(2).EQ.'REGULAR   ').OR. &
              (MSHTYP(2).EQ.'BLOCKWISE ')) THEN                                                    
         KTYPE = - KTYPE 
         BACKSPACE(fINP)
         IF (KTYPE.EQ. - 2) THEN 
            READ (fINP,*,iostat=ios)MSHSTR, NN1, NN2 
            IF (ios/=0)call ErrorIO('Error specifying NN1 and/or NN2')
            NN3 = 1 
         ELSE 
            READ (fINP,*,iostat=ios)MSHSTR, NN1, NN2, NN3 
            IF (ios/=0)call ErrorIO('Error specifying NN1, NN2, and/or NN3')
         ENDIF 
         IF ((NN1.LE.0).OR.(NN2.LE.0).OR.(NN3.LE.0)) THEN 
            call ErrorIO('Error NN1, NN2, and NN3(for 3D)must be > 0')
         ENDIF 
         IF (MSHTYP(2).EQ.'BLOCKWISE ') THEN 
            BACKSPACE(fINP)
            READ (fINP,*,iostat=ios)MSHSTR,(IDUM0, I0 = 1, - KTYPE),(IDUM1,     &
           (IDUM2, I2 = 1, IDUM1), I1 = 1, - KTYPE)                   
            IF (ios/=0)call ErrorIO('Error specifying BLOCKWISE mesh data')
         ENDIF 
      ELSE 
         call ErrorIO('Error specifying mesh data.  Must be either "REGULAR", "BLOCKWISE" or "IRREGULAR"')
      ENDIF 
!                                                                       
!.....OUTPUT DATASET 1                                                  
      WRITE (fLST, 180) TITLE1, TITLE2 
  180 FORMAT(////1X,131(1H-)//26X,80A1//26X,80A1//1X,131(1H-))
!                                                                       
!.....OUTPUT FILE UNIT ASSIGNMENTS                                      
      WRITE (fLST, 202) fINP, cINP, fICS, cICS, fLST, cLST
  202 FORMAT(/////11X,'F I L E   U N I T   A S S I G N M E N T S'//     &
     &   13X,'INPUT UNITS:'/                                            &
     &   13X,' SIMULATION DATA          ',I3,4X,'ASSIGNED TO ',A80/     &
     &   13X,' INITIAL CONDITIONS       ',I3,4X,'ASSIGNED TO ',A80//    &
     &   13X,'OUTPUT UNITS:'/                                           &
     &   13X,' SIMULATION RESULTS       ',I3,4X,'ASSIGNED TO ',A80)        
      IF (fRST.GT.0)WRITE (fLST, 9203)   fRST, cRST
 9203 FORMAT(13X,' RESTART DATA             ',I3,4X,'ASSIGNED TO ',A80)
      IF (fNOD.GT.0)WRITE (fLST, 9204)   fNOD, cNOD
 9204 FORMAT(13X,' PLOT DATA                ',I3,4X,'ASSIGNED TO ',A80)
      IF (fELE.GT.0)WRITE (fLST, 9206)   fELE, cELE
 9206 FORMAT(13X,' VECTOR PLOT DATA         ',I3,4X,'ASSIGNED TO ',A80)
      IF (fOBS.GT.0)WRITE (fLST, 9207)   fOBS, cOBS
 9207 FORMAT(13X,' OBSERVATION DATA         ',I3,4X,'ASSIGNED TO ',A80)
!.....VERSION 1.1
      IF(fTPN.GT.0)WRITE(fLST, 9208)   fTPN, cTPN
 9208 FORMAT(13X,' TECPLOT NODE DATA        ',I3,4X,'ASSIGNED TO ',A80)
      IF(fTPE.GT.0)WRITE(fLST, 9209)   fTPE, cTPE
 9209 FORMAT(13X,' TECPLOT ELEMENT DATA     ',I3,4X,'ASSIGNED TO ',A80)

      WRITE (fLST,203)
  203 FORMAT(//13X,'ADDITIONAL SUTRA-MS INPUT UNITS(OPTIONAL DATASETS):')
      IF (fTBC.GT.0)WRITE (fLST, 9211)   fTBC, cTBC
 9211 FORMAT(13X,' TRANSIENT BC DATA        ',I3,4X,'ASSIGNED TO ',A80)
      IF (fOTM.GT.0)WRITE (fLST, 9212) fOTM, cOTM
 9212 FORMAT(13X,' OUTPUT TIME DATA         ',I3,4X,'ASSIGNED TO ',A80)
      IF (fATS.GT.0)WRITE (fLST, 9213) fATS, cATS
 9213 FORMAT(13X,' AUTOMATIC TIME STEP      ',I3,4X,'ASSIGNED TO ',A80)
      IF (fZON.GT.0)WRITE (fLST, 9214) fZON, cZON
 9214 FORMAT(13X,' NODE AND ELEM. ZONES     ',I3,4X,'ASSIGNED TO ',A80)
      IF (fSOB.GT.0)WRITE (fLST, 9215) fSOB, cSOB
 9215 FORMAT(13X,' SPECIFIED OBS. PNTS.     ',I3,4X,'ASSIGNED TO ',A80)
!                                                                       
!.....INPUT DATASET 3A:  SIMULATION CONTROL NUMBERS                     
      MSErrorValue%cDataSet=' 3A'
      CALL SKPCOM(fINP, NLSKIP)
!.....SINGLE SPECIES SIMULATION                                         
      IF (NMULTI.EQ.0) THEN 
         NSPE = 1 
         READ (fINP,*,iostat=ios) NN, NE, NPBC, MNUBC, NSOP, MNSOU, NOBS 
         IF (ios/=0)call ErrorIO('Error specifying single species simulation control numbers')
!.....MULTI SPECIES SIMULATION                                          
      ELSE 
         READ (fINP,*,iostat=ios) NN, NE, NPBC, MNUBC, NSOP, MNSOU, NOBS, NSPE                                                           
         IF (ios/=0)call ErrorIO('Error specifying multi-species simulation control numbers')
      ENDIF 
!                                                                       
!.....TEST FOR ERRORS IN NUMBER OF SPECIFIED SPECIES                    
      IF (NSPE.LT.1) NSPE = 1 
      IF (ME.NE.0) GOTO 2010 
      IF (NSPE.GT.1) GOTO 2010 
!.....ERROR DETECTED IN INPUT FILE                                      
      call ErrorIO('Data Set 2A indicates multi-species transport but Data Set 3A indicates single species transport')
!                                                                       
!.....TEST TO SEE IF ENERGY TRANSPORT IS TO BE SIMULATED(ME=+1)       
!.....AND NS > 1                                                        
 2010 IF (ME.LT.1) GOTO 2020 
      IF (NSPE.LT.2) GOTO 2020 
!.....ERROR DETECTED IN INPUT FILE                                      
      call ErrorIO('Data Set 2A indicates energy transport but Data Set 3A indicates multi-species transport')
!                                                                       
 2020 IF (KTYPE.LT.0) THEN 
         NN123 = NN1 * NN2 * NN3 
         IF (NN123.NE.NN) THEN 
            call ErrorIO('Error: Number of NODES do not match rectangular dimensions(NN /= NN1*NN2*NN3)')
         ENDIF 
         IF (IABS(KTYPE).EQ.3) THEN 
            NE123 =(NN1 - 1)*(NN2 - 1)*(NN3 - 1)
         ELSE 
            NE123 =(NN1 - 1)*(NN2 - 1)
         ENDIF 
         IF (NE123.NE.NE) THEN 
            call ErrorIO('Error: Number of ELEMENTS do not match rectangular dimensions(NE /=(NN1-1)*(NN2-1)*(NN3-1))')
         ENDIF 
      ENDIF

      !check for consistency between NOBS and presence of fSOB(specified obs. locations)
      if (fSOB>0 .and. NOBS>0) &
        call ErrorIO('Error: Observation points cannot be specified in main input file if Specified Observation &
Locations are being used')
      !allocate storage and initialize SOB arrays
      if (fSOB>0) then
        write (*,*) 'Reading Specified Observation Location data...'
        IF (.not.Allo_Rd_SobData(NOBS,NOBCYC)) call ErrorIO('Error processing specified observation location data') 
      end if
      write (*,*) 'Initializing initial SUTRA storage...'
      !Allocate initial SUTRA storage for 
      !ITERAT, PARAMS, CONTRL, SOLVI, ITSOLI, ITSOLR, JCOLS, and MODSPR modules
      if (.not. AllocateInitialStorage()) call ErrorIO('Main2D3D: Error allocating inital SUTRA module storage')
!
!.....Read Zone Data
!     If Zone Parameter File is not specified(fZON==0)this routine is still called
!     since all spatial data for nodes and elements is contained in structures defined
!     and allocated in the zone module.  If the Zone Parameter File is not specified(fZON==0)
!     then NodeData and ElemData are dimensioned to NN and NE, respectively.
      write (*,*) 'Allocating hydraulic parameter zones and reading data if zone file is specified...'
      if (.not.RdZoneData()) call ErrorIO('Error in SutraZoneModule::RdZoneData')
!                                                                       
!.....INPUT AND OUTPUT DATASET 3B IF NSPE > 1                           
!.....READ 10 CHARACTER NAME FOR EACH SPECIES AND DETERMINE WHICH,      
!.....IF ANY ARE THE ENERY TRANSPORT SPECIES                            
      MSErrorValue%cDataSet=' 3B'
      NESP = 0 
!.....SOLUTE TRANSPORT ONLY                                             
      IF (ME.EQ. - 1) THEN 
         IF (NSPE.LT.2) THEN 
            SPNAME(1) = 'Concentration'
            NSOU(1)= MNSOU 
            NUBC(1)= MNUBC 
!.........SET NBCN BASED ON NUBC(1)FOR EFFICIENT STORAGE               
            NBCN = NPBC + MNUBC + 1 
!.........SET MNSOU BASED ON NSOU(1)FOR EFFICIENT STORAGE              
            MNSOU = MNSOU 
         ELSE 
            SPNAME = '    SOLUTE' 
            NSOU = MNSOU 
            NUBC = MNUBC 
!.........SKIP COMMENTS                                                 
            CALL SKPCOM(fINP, NLSKIP)
            WRITE (fLST, 2025)
!                                                                       
            NNU = 0 
            NNS = 0 
 2050       READ (fINP,*,iostat=ios) K 
            IF (ios/=0) call ErrorIO('Error: specifying species number')

!.........TERMINATE DATASET 3B READ IF K EQUALS ZERO                    
            IF (K.EQ.0)GOTO 2055 
!.........CONTINUE READING TILL K EQUALS ZERO IF                        
!         K EXCEEDS THE NUMBER OF SPECIES                               
            IF (K.GT.NSPE)GOTO 2050 
!.........BACKSPACE AND RE-READ LINE                                    
            BACKSPACE(fINP)
            READ (fINP,*,iostat=ios)K, SPNAME(K), NUBC(K), NSOU(K)
            IF (ios/=0) call ErrorIO('Error specifying data for species '//trim(adjustl(Val2Char(K))))

            WRITE (fLST, 2035)K, SPNAME(K), NUBC(K), NSOU(K)
!.........COUNTERS FOR MAXIMUM NUMBER OF U BOUNDARY CONDITIONS          
            NNU = MAX(NNU, NUBC(K))
            NNS = MAX(NNS, NSOU(K))
!.........CONTINUE READING UNTIL K EQUALS ZERO                          
            GOTO 2050 
!.........SET MNUBC BASED ON NUBC(K)FOR EFFICIENT STORAGE              
 2055       MNUBC = NNU 
!.........SET NBCN BASED ON NUBC(K)FOR EFFICIENT STORAGE               
            NBCN = NPBC + NNU + 1 
!.........SET MNSOU BASED ON NSOU(K)FOR EFFICIENT STORAGE              
            MNSOU = NNS 
         ENDIF 
!                                                                       
!.....ENERGY TRANSPORT ONLY                                             
      ELSEIF (ME.EQ.1) THEN 
         NESP = 1 
         SPNAME(1) = 'Temperature' 
         NSOU(1)= MNSOU 
         NUBC(1)= MNUBC 
!.......SET NBCN BASED ON NUBC(1)FOR EFFICIENT STORAGE                 
         NBCN = NPBC + MNUBC + 1 
!.......SET MNSOU BASED ON NSOU(1)FOR EFFICIENT STORAGE                
         MNSOU = MNSOU 
!                                                                       
!.....SOLUTE AND ENERGY TRANSPORT SPECIFIED(ME=+0)                    
      ELSE 
         SPNAME = '    SOLUTE' 
         NSOU = MNSOU 
         NUBC = MNUBC 
!.......SKIP COMMENTS                                                   
         CALL SKPCOM(fINP, NLSKIP)
         WRITE (fLST, 2030)
         NNU = 0 
         NNS = 0 
!.......READ SPECIES NUMBER                                             
 2060    READ (fINP,*,iostat=ios)K 
         IF (ios/=0)call ErrorIO('Error: specifying species number')
!.......TERMINATE DATASET 3B READ IF K EQUALS ZERO                      
         IF (K.EQ.0)GOTO 2065 
!.......CONTINUE READING TILL K EQUALS ZERO IF                          
!       K EXCEEDS THE NUMBER OF SPECIES                                 
         IF (K.GT.NSPE)GOTO 2060 
!.......BACKSPACE AND RE-READ LINE                                      
         BACKSPACE(fINP)
!                                                                       
         READ (fINP,*,iostat=ios)K, SPNAME(K), NUBC(K), NSOU(K)
         IF (ios/=0)call ErrorIO('Error specifying data for species '//trim(adjustl(Val2Char(K))))
         WRITE (fLST, 2035)K, SPNAME(K), NUBC(K), NSOU(K)
!.......COUNTERS FOR MAXIMUM NUMBER OF U BOUNDARY CONDITIONS            
         NNU = MAX(NNU, NUBC(K))
         NNS = MAX(NNS, NSOU(K))
!.......DETERMINE WHICH SPECIES IS ENERGY FOR SOLUTE AND ENERGY         
!.......TRANSPORT ARE SPECIFIED(ME=+0)                                
         IF (INDEX(SPNAME(K), 'ENERGY')>0) THEN 
!.........ERROR IF MORE THAN ONE ENERGY TRANSPORT SPECIES               
!.........IS SPECIFIED                                                  
            IF (NESP.GT.0) THEN 
               call ErrorIO('Error: more than one energy transport species specified')
            ENDIF 
            NESP = K 
            SPNAME(K)= 'Temperature'
         ENDIF 
!.......CONTINUE READING UNTIL K EQUALS ZERO                            
         GOTO 2060 
!.......SET NBCN BASED ON NUBC(K)FOR EFFICIENT STORAGE                 
 2065    NBCN = NPBC + NNU + 1 
!.......SET MNUBC BASED ON NUBC(K)FOR EFFICIENT STORAGE                
         MNUBC = NNU 
!.......SET MNSOU BASED ON NSOU(K)FOR EFFICIENT STORAGE                
         MNSOU = NNS 
!.......CONFIRM THAN AT LEAST ONE SPECIES IS ENERGY TRANSPORT SPECIES   
         IF (NESP.LT.1) THEN 
            call ErrorIO('Error: Energy and solute transport is specified but no energy species is specified')
         ENDIF 
      ENDIF 
!                                                                       
 2025 FORMAT(//,19X,                                                    &
     &  'DATASET 3B SOLUTE TRANSPORT SPECIES',/11X,                     &
     &   '   SPECIES',10X,'   NAME',19X,'NUBC',16X,'NSOU',/,3X,         &
     &   66('-'),/)                                                    
 2030 FORMAT(//,14X,                                                    &
     &  'DATASET 3B SOLUTE AND ENERGY TRANSPORT SPECIES',/11X,          &
     &  '   SPECIES',10X,'   NAME',19X,'NUBC',16X,'NSOU',/,3X,          &
     &   66('-'),/)                                                    
 2035 FORMAT(11X,I10,10X,A10,10X,I10,10X,I10)
 2040 FORMAT(//,10X,                                                    &
     &'MORE THAN ONE ENERGY SPECIES INDICATED IN DATASET 3B.',/,11X,    &
     &'ONLY ONE ENERGY TRANSPORT SPECIES ALLOWED.',/,12X,               &
     &'PLEASE CORRECT THIS IN THE INPUT DATA AND RERUN.',/,9X,          &
     &'S I M U L A T I O N   H A L T E D   DUE TO INPUT ERROR')        
 2045 FORMAT(//,15X,                                                    &
     &'SOLUTE AND ENERGY TRANSPORT IS SPECIFIED BUT',/,13X,             &
     &'NO ENERGY SPECIES IS SPECIFIED IN DATASET 3B.',/,12X,            &
     &'PLEASE CORRECT THIS IN THE INPUT DATA AND RERUN.',/,9X,          &
     &'S I M U L A T I O N   H A L T E D   DUE TO INPUT ERROR')        
!                                                                       
!.....INPUT AND OUTPUT DATASET 4:  SIMULATION MODE OPTIONS              
      write (*,*) 'Reading Data Set 4...'
      MSErrorValue%cDataSet='  4'
      CALL SKPCOM(fINP, NLSKIP)
      READ (fINP,*,iostat=ios) UNSSTR, SSFSTR, SSTSTR, RDSTR, ISTORE 
      IF (ios/=0) call ErrorIO('Error specifying simulation mode options')
!                                                                       
      CALL PRSWDS(UNSSTR, 1, CUNSAT, NWORDS)
      CALL PRSWDS(SSFSTR, 1, CSSFLO, NWORDS)
      CALL PRSWDS(SSTSTR, 1, CSSTRA, NWORDS)
      CALL PRSWDS(RDSTR,  1,  CREAD, NWORDS)
      ISMERR = 0 
      IF (CUNSAT(1).EQ.'UNSATURATED') THEN 
         IUNSAT = + 1 
      ELSEIF (CUNSAT(1).EQ.'SATURATED') THEN 
         IUNSAT = 0 
      ELSE 
         ISMERR = 1 
      ENDIF 
      IF (CSSFLO(1).EQ.'TRANSIENT') THEN 
         ISSFLO = 0 
      ELSEIF (CSSFLO(1).EQ.'STEADY') THEN 
         ISSFLO = + 1 
      ELSE 
         ISMERR = 1 
      ENDIF 
      IF (CSSTRA(1).EQ.'TRANSIENT') THEN 
         ISSTRA = 0 
      ELSEIF (CSSTRA(1).EQ.'STEADY') THEN 
         ISSTRA = + 1 
      ELSE 
         ISMERR = 1 
      ENDIF 
      IF (CREAD (1).EQ.'COLD') THEN 
         IREAD = + 1 
      ELSEIF (CREAD (1).EQ.'WARM') THEN 
         IREAD = - 1 
      ELSE 
         ISMERR = 1 
      ENDIF 
      IF (ISMERR.EQ.1) THEN 
         call ErrorIO('Error: unrecognized simulation mode')
      ENDIF 
      WRITE (fLST, 205)
  205 FORMAT(////11X,'S I M U L A T I O N   M O D E   ',                &
     &   'O P T I O N S'/)                                             
      IF (ISSTRA.EQ.1.AND.ISSFLO.NE.1) THEN 
         call ErrorIO('Error: steady state transport(CSSTRA="STEADY" requires CSSFLO="STEADY"')
      ENDIF 
      IF (IUNSAT.EQ. + 1)WRITE (fLST, 215)
      IF (IUNSAT.EQ.0)WRITE (fLST, 216)
  215 FORMAT(11X,'- ALLOW UNSATURATED AND SATURATED FLOW:  UNSATURATED',&
     &   ' PROPERTIES ARE USER-PROGRAMMED IN SUBROUTINE   U N S A T')  
  216 FORMAT(11X,'- ASSUME SATURATED FLOW ONLY')
      IF (ISSFLO.EQ. + 1.AND.ME.EQ. - 1)WRITE (fLST, 219)
      IF (ISSFLO.EQ. + 1.AND.ME.EQ. + 1)WRITE (fLST, 220)
      IF (ISSFLO.EQ.0)WRITE (fLST, 221)
  219 FORMAT(11X,'- ASSUME STEADY-STATE FLOW FIELD CONSISTENT WITH ',   &
     &   'INITIAL CONCENTRATION CONDITIONS')                           
  220 FORMAT(11X,'- ASSUME STEADY-STATE FLOW FIELD CONSISTENT WITH ',   &
     &   'INITIAL TEMPERATURE CONDITIONS')                             
  221 FORMAT(11X,'- ALLOW TIME-DEPENDENT FLOW FIELD')
      IF (ISSTRA.EQ. + 1)WRITE (fLST, 225)
      IF (ISSTRA.EQ.0)WRITE (fLST, 226)
  225 FORMAT(11X,'- ASSUME STEADY-STATE TRANSPORT')
  226 FORMAT(11X,'- ALLOW TIME-DEPENDENT TRANSPORT')
      IF (IREAD.EQ. - 1)WRITE (fLST, 230)
      IF (IREAD.EQ. + 1)WRITE (fLST, 231)
  230 FORMAT(11X,'- WARM START - SIMULATION IS TO BE ',                 &
     &   'CONTINUED FROM PREVIOUSLY-STORED DATA')                      
  231 FORMAT(11X,'- COLD START - BEGIN NEW SIMULATION')
      IF (ISTORE.GT.0)WRITE (fLST, 240)ISTORE, fRST 
      IF (ISTORE.EQ.0)WRITE (fLST, 241)
  240 FORMAT(11X,'- STORE RESULTS AFTER EVERY',I5,' TIME STEPS ON UNIT-'&
     &   I3,' AS BACK-UP AND FOR USE IN A SIMULATION RE-START')        
  241 FORMAT(11X,'- DO NOT STORE RESULTS FOR USE IN A ',                &
     &   'RE-START OF SIMULATION')                                     
!.....OUTPUT DATASET 3                                                  
!.....SOLUTE TRANSPORT - MUTIPLE SPECIES ALLOWED                        
      IF (ME.EQ. - 1) THEN 
         WRITE (fLST, 245) NN, NE, NBI, NPBC, MNUBC, NSOP, MNSOU, NOBS 
         WRITE (fLST, 246) NSPE 
      ENDIF 
  245 FORMAT(////11X,'S I M U L A T I O N   C O N T R O L   ',          &
     &   'N U M B E R S'// 8X,I9,5X,'NUMBER OF NODES IN FINITE-',       &
     &   'ELEMENT MESH'/ 8X,I9,5X,'NUMBER OF ELEMENTS IN MESH'/         &
     &    8X,I9,5X,'ESTIMATED MAXIMUM FULL BAND WIDTH FOR MESH'//       &
     &    8X,I9,5X,'EXACT NUMBER OF NODES IN MESH AT WHICH ',           &
     &   'PRESSURE IS A SPECIFIED CONSTANT OR FUNCTION OF TIME'/        &
     &    8X,I9,5X,'MAXIMUM NUMBER OF NODES IN MESH AT WHICH ',         &
     &   'SOLUTE CONCENTRATION IS A SPECIFIED CONSTANT OR ',            &
     &   'FUNCTION OF TIME'// 8X,I9,5X,'MAXIMUM NUMBER OF NODES AT',    &
     &   ' WHICH FLUID INFLOW OR OUTFLOW IS A SPECIFIED CONSTANT',      &
     &   ' OR FUNCTION OF TIME'/ 8X,I9,5X,'MAXIMUM NUMBER OF NODES AT', &
     &   ' WHICH A SOURCE OR SINK OF SOLUTE MASS IS A SPECIFIED ',      &
     &   'CONSTANT OR FUNCTION OF TIME'// 8X,I9,5X,'EXACT NUMBER OF ',  &
     &   'NODES AT WHICH PRESSURE AND CONCENTRATION WILL BE OBSERVED') 
!                                                                       
  246 FORMAT(                                                           &
     &   8X,I9,5X,'NUMBER OF SPECIES FOR WHICH TRANSPORT WILL ',        &
     &            'BE SIMULATED')                                      
!                                                                       
!.....ENERGY TRANSPORT                                                  
      IF (ME.EQ. + 1) THEN 
         WRITE (fLST, 247)NN, NE, NBI, NPBC, MNUBC, NSOP, MNSOU, NOBS 
         WRITE (fLST, 246)NSPE 
      ENDIF 
  247 FORMAT(////11X,'S I M U L A T I O N   C O N T R O L   ',          &
     &   'N U M B E R S'// 8X,I9,5X,'NUMBER OF NODES IN FINITE-',       &
     &   'ELEMENT MESH'/ 8X,I9,5X,'NUMBER OF ELEMENTS IN MESH'/         &
     &    8X,I9,5X,'ESTIMATED MAXIMUM FULL BAND WIDTH FOR MESH'//       &
     &    8X,I9,5X,'EXACT NUMBER OF NODES IN MESH AT WHICH ',           &
     &   'PRESSURE IS A SPECIFIED CONSTANT OR FUNCTION OF TIME'/        &
     &    8X,I9,5X,'EXACT NUMBER OF NODES IN MESH AT WHICH ',           &
     &   'TEMPERATURE IS A SPECIFIED CONSTANT OR ',                     &
     &   'FUNCTION OF TIME'// 8X,I9,5X,'EXACT NUMBER OF NODES AT',      &
     &   ' WHICH FLUID INFLOW OR OUTFLOW IS A SPECIFIED CONSTANT',      &
     &   ' OR FUNCTION OF TIME'/ 8X,I9,5X,'EXACT NUMBER OF NODES AT',   &
     &   ' WHICH A SOURCE OR SINK OF ENERGY IS A SPECIFIED CONSTANT',   &
     &   ' OR FUNCTION OF TIME'// 8X,I9,5X,'EXACT NUMBER OF NODES ',    &
     &   'AT WHICH PRESSURE AND TEMPERATURE WILL BE OBSERVED')         
!                                                                       
!.....SOLUTE AND ENERGY TRANSPORT                                       
      IF (ME.EQ.0) THEN 
         WRITE (fLST, 248)NN, NE, NBI, NPBC, MNUBC, NSOP, MNSOU, NOBS 
         WRITE (fLST, 246)NSPE 
      ENDIF 
  248 FORMAT(////11X,'S I M U L A T I O N   C O N T R O L   ',          &
     &   'N U M B E R S'// 8X,I9,5X,'NUMBER OF NODES IN FINITE-',       &
     &   'ELEMENT MESH'/ 8X,I9,5X,'NUMBER OF ELEMENTS IN MESH'/         &
     &    8X,I9,5X,'ESTIMATED MAXIMUM FULL BAND WIDTH FOR MESH'//       &
     &    8X,I9,5X,'EXACT NUMBER OF NODES IN MESH AT WHICH ',           &
     &   'PRESSURE IS A SPECIFIED CONSTANT OR FUNCTION OF TIME'/        &
     &    8X,I9,5X,'MAXIMUM NUMBER OF NODES IN MESH AT WHICH ',         &
     &   'SOLUTE CONCENTRATION OR TEMPERATURE IS A SPECIFIED ',         &
     &   'CONSTANT OR FUNCTION OF TIME'// 8X,I9,5X,'MAXIMUM NUMBER '    &
     &   'OF NODES AT WHICH FLUID INFLOW OR OUTFLOW IS A SPECIFIED ',   &
     &   'CONSTANT OR FUNCTION OF TIME'/ 8X,I9,5X,'MAXIMUM NUMBER ',    &
     &   'OF NODES AT WHICH A SOURCE OR SINK OF SOLUTE MASS OR '        &
     &   'ENERGY IS A SPECIFIED CONSTANT OR FUNCTION OF TIME'//         &
     &    8X,I9,5X,'EXACT NUMBER OF NODES AT WHICH PRESSURE, ',         &
     &   'CONCENTRATION, AND TEMPERATURE WILL BE OBSERVED')            
!                                                                       
!.....INPUT DATASETS 5 - 7                                              
      write (*,*) 'Reading Data Sets 5-7...'
      CALL INDAT0()
!.....KSOLVP AND KSOLVU HAVE BEEN SET ACCORDING TO THE SOLVERS SELECTED:
!        INVALID SELECTION                      ==>  -1                 
!        BANDED GAUSSIAN ELIMINATION(DIRECT)    ==>   0                 
!        IC-PRECONDITIONED CONJUGATE GRADIENT   ==>   1                 
!        ILU-PRECONDITIONED GMRES               ==>   2                 
!        ILU-PRECONDITIONED ORTHOMIN            ==>   3                 
!        Additional solvers                     ==>   4                 
!                                                                       
!.....OUTPUT DATASETS 7B AND 7C                                         
      WRITE (fLST, 261)
  261 FORMAT(////11X,'S O L V E R - R E L A T E D   ',                  &
     &   'P A R A M E T E R S')                                        
!                                                                       
      DO 271 K = 1, NSPE 
         IF ((KSOLVP.LT.0).OR.(KSOLVU(K).LT.0).OR. &
           ((KSOLVP * KSOLVU(K).EQ.0).AND.(KSOLVP + KSOLVU(K).NE.0))) THEN   
            WRITE (fLST, '(////11X,a)') 'VALID SOLVER OPTIONS ARE: '
            WRITE (fLST, 263)(SOLWRD(M), SOLNAM(M), M = 0, NSLVRS - 1)
            WRITE (fLST, 264)
            call ErrorIO('Error: INVALID SOLVER SELECTION')
  263 FORMAT(11X,A10,' --> ',A40)
  264 FORMAT(/11X,'SOLVER SELECTIONS FOR P AND U MUST BE BOTH',        &
     &    ' DIRECT OR BOTH ITERATIVE.'                                  &
     &    //11X,'PLEASE CORRECT THE INPUT DATA AND RERUN.'////////      &
     &    45X,'S I M U L A T I O N   H A L T E D   DUE TO INPUT ERROR')
         ENDIF 
  271 END DO 
!                                                                       
      IF ((KTYPE.GT.0).AND.(KSOLVP.NE.0) .and. (.not.lColumnStorage)) THEN 
         call ErrorIO('Error: Currently the direct solver must be used for the irregular meshes if not using AP column storage')
      ENDIF 

!.....OUTPUT DATASETS 3B & 3C                                           
  266 IF (KSOLVP.NE.0) THEN 
         WRITE (fLST, 268) NN1, NN2, NN3, SOLNAM(KSOLVP), ITRMXP, ITOLP, &
         TOLP, NSAVEP                                                   
         DO K = 1, NSPE 
           WRITE (fLST, 269) TRIM(ADJUSTL(SPNAME(K))), SOLNAM(KSOLVU(K)), &
                             ITRMXU(K), ITOLU(K), TOLU(K), StartUTime(K), NSAVEU(K)
         ENDDO 
  268 FORMAT(                                                           &
     &   /8X,I9,5X,'NUMBER OF NODES IN 1ST NUMBERING DIRECTION'         &
     &   /8X,I9,5X,'NUMBER OF NODES IN 2ND NUMBERING DIRECTION'         &
     &   /8X,I9,5X,'NUMBER OF NODES IN 3RD NUMBERING DIRECTION'         &
     &   //13X,'SOLVER FOR P: ',A40                                     &
     &   //20X,I6,5X,'MAXIMUM NUMBER OF MATRIX SOLVER ITERATIONS',      &
     &        ' DURING P SOLUTION'                                      &
     &   /20X,I6,5X,'TYPE OF CONVERGENCE CRITERION FOR MATRIX',         &
     &        ' SOLVER ITERATIONS DURING P SOLUTION',                   &
     &   /11X,1PD15.4,5X,'CONVERGENCE TOLERANCE FOR MATRIX',            &
     &        ' SOLVER ITERATIONS DURING P SOLUTION'                    &
     &   /20X,I6,5X,'NUMBER OF DIRECTION VECTORS DURING P SOLUTION',    &
     &        '(USED BY GMRES AND ORTHOMIN)')                         
  269 FORMAT(                                                           &
     &   //13X,'SOLVER FOR U [',A,']: ',A40                             &
     &   //20X,I6,5X,'MAXIMUM NUMBER OF MATRIX SOLVER ITERATIONS',      &
     &        ' DURING U SOLUTION'                                      &
     &   /20X,I6,5X,'TYPE OF CONVERGENCE CRITERION FOR MATRIX',         &
     &        ' SOLVER ITERATIONS DURING U SOLUTION'                    &
     &   /11X,1PD15.4,5X,'CONVERGENCE TOLERANCE FOR MATRIX',            &
     &        ' SOLVER ITERATIONS DURING U SOLUTION'                    &
     &   /11X,1PD15.4,5X,'START TIME FOR U SOLUTION (-1e6 = TSTART)',   &
     &   /20X,I6,5X,'NUMBER OF DIRECTION VECTORS DURING U SOLUTION',    &
     &        '(USED BY GMRES AND ORTHOMIN)' )                        
      ELSE 
         WRITE (fLST, 270)SOLNAM(KSOLVP)
  270 FORMAT(/13X,'SOLVER FOR P AND U: ',A40)
      ENDIF 
!                                                                       
!                                                                       
!.....FOR CG OR ORTHOMIN SOLVER, IF ITOL=0, DEFAULT TO ITOL=1.          
      IF (((KSOLVP.EQ.1).OR.(KSOLVP.EQ.3)).AND.(ITOLP.EQ.0)) ITOLP = 1
      DO K = 1, NSPE 
        IF (((KSOLVU(K).EQ.1).OR.(KSOLVU(K).EQ.3)).AND.(ITOLU(K).EQ.0)) ITOLU(K) = 1
      ENDDO 
!                                                                       
!                                                                       
!.....CALCULATE THE NUMBER OF TIME STEPS ON WHICH OBSERVATIONS WILL     
!.....BE MADE, NTOBS.  THIS REQUIRES LOOKING AHEAD TO DATASET 8 OF      
!.....UNIT fINP AND DATASET 1 OF UNIT fICS.  THE UNITS ARE THEN "BACKSPACED"
!.....SO THAT THESE DATASETS CAN BE READ AGAIN LATER BY SUBROUTINES     
!.....INDAT1 AND INDAT2.                                                
!                                                                       
      MSErrorValue%cDataSet='Ad8'
      NTOBS = 0

      IF (NOBS.EQ.0)GOTO 11 
      IF (fSOB==0) then
          NLSTOT = 0 
          DO 7 I = 1, 4 
             CALL SKPCOM(fINP, NLSKIP)
             NLSTOT = NLSTOT + NLSKIP 
             READ (fINP,*,iostat=ios) NOBCYC 
             IF (ios/=0) call ErrorIO('Error specifying NOBCYC - in look ahead to Data Set 8')
    7   END DO 
        DO 9 I = 1, NLSTOT + 4 
    9   BACKSPACE(fINP)
      end if
      MSErrorValue%cDataSet='ICS'
      CALL SKPCOM(fICS, NLSKIP)
      READ (fICS,*,iostat=ios) TSTART 
      IF (ios/=0) call ErrorIO('Error specifying TSTART - in look ahead to starting time in ICS file')
      BACKSPACE(fICS)

   11 CONTINUE 
!                                                                       
!.....CALCULATE PROBLEM DIMENSIONS, EXCLUDING SOLVER DIMENSIONS
!     PARAMETER ARRAY LENGTHS                                           
      NSOP = NSOP + 1 
      MNSOU = MNSOU + 1 
      NIN = NE * 8 
      NOBSN = NOBS + 1 
      NSPEDIM = NN * NSPE 
!.....3-D MESH                                                          
      IF (IABS(KTYPE).EQ.3) THEN 
         NEV = NEV3 
         N48 = 8 
         NEX = NE
         NEZV = NEZV3
!.....2-D MESH                                                          
      ELSE 
         NEV = NEV2 
         N48 = 4 
         NEX = 1 
         NEZV = NEZV2
      ENDIF 
!.....NEXV IS THE NUMBER OF VECTORS THAT ARE OF LENGTH NE IN 3-D AND    
!.....ARE TO BE DIMENSIONED TO LENGTH 1 BECAUSE THEY ARE NOT NEEDED.    
!.....THUS, IN 3-D, NEXV=0; IN 2-D, NEXV=NEV3-NEV2.                     
      NEVX = NEV3 - NEV 
      NEVG = 2 * N48 
      NE8 = NE * N48 
      NE8X = NEX * N48 

      !Set Observation times
      if (NOBS>0) then
        if (.not. AllocateStandardSutraObservations() ) &
          call ErrorIO('OBS::AllocateStandardSutraObservations Error could not allocate and initialize OBS')
      end if
!                                                                       
!.....GET AUTOMATIC TIME STEPPING DATA 
      IF (fATS.GT.0) THEN 
         IF (.not. MAutomaticTimeStep()) call ErrorIO('Error reading ATS package')
      ENDIF 

!!     Dimension Main Arrays
      IF (.NOT. AllocateMainStorage()) call ErrorIO('Main2D3D: Error allocating primary SUTRA storage')

!.....INPUT SIMULATION DATA FROM UNIT-fINP (DATASETS 8 THROUGH 15)        
      write (*,*) 'Reading Data Sets 8-15...'
      CALL INDAT1()
!                                                                       
!.....INPUT FLUID MASS, AND ENERGY OR SOLUTE MASS SOURCES               
!        (DATASETS 17 AND 18)                                           
      CALL ZERO (QIN, NN, 0.0D0) 
      CALL ZERO (UIN, NN * NSPE, 0.0D0) 
      CALL ZERO (QUIN, NN * NSPE, 0.0D0) 
      IF (NSOP - 1.GT.0.OR.MNSOU - 1.GT.0) &
        CALL SOURCE (QIN, UIN, IQSOP, QUIN, IQSOU, IQSOPT, IQSOUT)                                      
!                                                                       
!.....INITIALIZE PBC                                                    
      PBC = 0.0 
!                                                                       
!.....INPUT SPECIFIED P AND U BOUNDARY CONDITIONS (DATASETS 19 AND 20)  
      IF (NBCN - 1.GT.0) CALL BOUND (IPBC, PBC, IUBC, UBC, IPBCT, IUBCT) 
!                                                                       
!.....SET FLAG FOR TIME-DEPENDENT SOURCES OR BOUNDARY CONDITIONS.       
!     WHEN IBCT=+4, THERE ARE NO TIME-DEPENDENT SPECIFICATIONS.         
      IBCT = IQSOPT + IQSOUT + IPBCT + IUBCT 
!                                                                       
!.....INPUT MESH CONNECTION DATA (DATASET 22)                           
      write (*,*) 'Reading Data Set 22 (Mesh Connection Data)...'
      CALL CONNEC (IN)
!                                                                       
!.....CALCULATE AND CHECK BAND WIDTH                                    
      write (*,*) 'Calculating Problem Bandwidth...'
      CALL BANWID (IN) 

!                                                                       
!.....CALCULATE SOLVER DIMENSIONS
      IF (KSOLVP.EQ.0) THEN 
!........SET DIMENSIONS FOR DIRECT SOLVER                               
         NCBI = NBI 
         NBIX = 1 
         NNNX = 1 
         NIPARM = 1 
         NFPARM = 1 
         NELT = NN 
         NWI = 1 
         NWF = 1 
      ELSE 
!........SET DIMENSIONS FOR ITERATIVE SOLVER(S)                        
         NCBI = 1 
         NBIX = NBI 
         NNNX = NN 
         NIPARM = 1 
         NFPARM = 1 
!........2-D MESH                                                       
         IF (IABS(KTYPE).EQ.2) THEN 
            NELT = 9 * NN - 6 * NN1 - 2 
!........3-D MESH                                                       
         ELSE 
            NELT = 27 * NN - 6 * NN1 *(1 + 3 * NN2)- 2 
         ENDIF 
         CALL DIMWRK(KSOLVP, NSAVEP, NN, NELT, NWIP, NWFP)
!........MAXIMUM DIMENSION OF U ITERATIVE SOLUTION                      
         KMXSOLVU = MAXVAL(KSOLVU)
         NMXSAVEU = MAXVAL(NSAVEU)
         CALL DIMWRK(KMXSOLVU, NMXSAVEU, NN, NELT, NWIU, NWFU)
         NWI = MAX(NWIP, NWIU)
         NWF = MAX(NWFP, NWFU)
      ENDIF 
      MATDIM=NELT*NCBI
!.....END SOLVER POINTERS
!

!
!.....CALCULATE IA AND JA IN COLUMN STORAGE FORMAT IF REQUIRED
      IF ( KSOLVP.EQ.0 ) then
!.......ALLOCATE STORAGE FOR SOLVER NOW USING SOLVER SIZES ESTIMATED
!         WITH ORIGINAL SUTRA2D3D METHODS IF NOT USING NEW COLUMSTORAGE INDEX ROUTINES
         WRITE (*,*) 'Allocating storage for solvers...'
         MSErrorValue%cDataSet='TRI'
         IF(.NOT.AllocateSolverStorage()) &
           call ErrorIO('AllocateSolverStorage: Could not allocate storage for iterative solvers (triad storage)')
!         IF (.not.AllocateTriad2Column()) &
!           call ErrorIO('AllocateTriad2Column(): Could not allocate storage for conversion from triad to column format')
      ELSE
        IF ( .not.lColumnStorage ) THEN
           write (*,*) 'Allocating storage for solvers...'
           MSErrorValue%cDataSet='TRI'
           IF (.not.AllocateSolverStorage()) &
             call ErrorIO('AllocateSolverStorage: Could not allocate storage for iterative solvers (triad storage)')
           write (*,*) 'Allocating storage for Triad to Column conversion subroutines...'
           IF (.not.AllocateTriad2Column()) &
             call ErrorIO('AllocateTriad2Column(): Could not allocate storage for conversion from triad to column format')
        ELSE
          write (*,*) 'Allocating solver storage and making AP column storage format...'
          if(.not.MakeITRIJTRI()) &
            call ErrorIO('MakeITRIJTRI:: Could not create column storage format')
        END IF
      END IF
!                                                                       
!.....SET UP ARRAY NBI27, WHICH GIVES THE TRANSFORMATION BETWEEN        
!     SUTRA'S BANDED MATRIX FORMAT AND A COMPRESSED FORMAT BASED        
!     ON THE 27-NODE "MOLECULE".  ALSO, SET UP POINTERS RELATED         
!     TO THE "SLAP TRIAD" MATRIX STRUCTURE.
      IF ( .not.lColumnStorage .and. KSOLVP.NE.0 ) then
        write (*,*) 'Setting up SLAP Triad Pointers...'
        CALL PTRSET (NBI27, ITRI, JTRI, MIOFF)
      END IF
!.....MOVED FROM SutraMSSubroutine - VERSION 1.1                                                                       
!.....INPUT INITIAL OR RESTART CONDITIONS AND INITIALIZE PARAMETERS     
!        (READ UNIT-fICS DATA)                                            
      write (*,*) 'Reading ICS data...'
      CALL INDAT2 (PVEC, UVEC, PM1, UM1, UM2, CS1, CS2, CS3, SL, SR,    &
                   RCIT, SW, DSWDP, PBC, IPBC, IPBCT)                          
!.....VERSION 1.1
!.....Allocate storage for Tecplot, if required
      IF (fTPN>0 .or. fTPE>0) THEN
        write (*,*) 'Allocating storage for Tecplot output...'
        MSErrorValue%cDataSet='TPL'
        if (.not. AllocateTecplot() ) call ErrorIO('Tecplot::AllocateTecplot:: Could not allocate required Tecplot arrays')
      END IF
!
!.....Calculate total storage used
      write (*,*) 'Calculating and output of total storage allocated...'
      IF (.not. CalculateTotalStorage() ) &
          call ErrorIO('TotalStorage::CalculateTotalStorage Could not write summary of memory allocated at runtime')

!                                                                       
!.....MAIN CONTROL ROUTINE, SUTRA                                       
      write (*,*) 'Calling main SUTRA subroutine...'
      CALL SUTRA(TITLE1,TITLE2)
!                                                                      
!
!.....TOTAL SIMULATION TIME
      CALL CPU_TIME(RT1)
      RTF=RT1-RT0
!
      WRITE (*,*) 'TOTAL SIMULATION TIME: ',RTF,' SEC.'
      WRITE (fLST,*)'TOTAL SIMULATION TIME: ',RTF,' SEC.'
      RTF=RTF/60.0D0
      WRITE (*,*) 'TOTAL SIMULATION TIME: ',RTF,' MIN.'
      WRITE (fLST,*)'TOTAL SIMULATION TIME: ',RTF,' MIN.'
      RTF=RTF/60.0D0
      WRITE (*,*) 'TOTAL SIMULATION TIME: ',RTF,'  HR.'
      WRITE (fLST,*)'TOTAL SIMULATION TIME: ',RTF,'  HR.'
      RTF=RTF/24.0D0
      WRITE (*,*) 'TOTAL SIMULATION TIME: ',RTF,' DAYS'
      WRITE (fLST,*)'TOTAL SIMULATION TIME: ',RTF,' DAYS'
!                                                                       
      LNormal=.true.
      call SutraEnd()
      END PROGRAM SUTRA_MS
