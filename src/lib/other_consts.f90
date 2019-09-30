!> \file
!> \author Merryn Tawhai
!> \brief This module contains definition of constants
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>
!>\Description
!> This module contains definition of constants (note that in the future this should be merged into a 'types' module
MODULE other_consts
  USE arrays, ONLY: dp
  IMPLICIT NONE

  INTEGER, PARAMETER :: MAX_FILENAME_LEN = 255, MAX_STRING_LEN = 100, MAX_SUBNAME_LEN = 60

  REAL(dp), PARAMETER :: PI = 3.14159265358979_dp
  REAL(dp), PARAMETER :: TOLERANCE =EPSILON (EPSILON(1.0_dp))

  PRIVATE
  PUBLIC MAX_SUBNAME_LEN, MAX_STRING_LEN, MAX_FILENAME_LEN, PI,TOLERANCE
END MODULE other_consts
