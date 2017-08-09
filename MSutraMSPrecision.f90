
!     MODULE            SUTRAMSPRECISION       SUTRA-MS VERSION 2004.1
! ***
! *** PURPOSE :                                                         
! ***  MODULE DEFINES THE PRECISION FOR 2 AND 4 BYTE INTEGERS AND
! ***  SINGLE (4 BYTE) AND DOUBLE (8 BYTE) PRECISION REAL VARIABLES.
! ***  HARDWARE INDEPENDENT IMPLEMENTATION OF TYPE PRECISION.
! ***

  module SutraMSPrecision

    integer, parameter :: I2B = selected_int_kind(4)
    integer, parameter :: I4B = selected_int_kind(9)

    integer, parameter :: SP  = selected_real_kind(6,37)
    integer, parameter :: DP  = selected_real_kind(15,307)

  end module SutraMSPrecision
