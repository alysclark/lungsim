!> \file
!> \author Merryn Tawhai
!> \brief This module handles diagnostics.
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>
!>\Description
!> This module handles diagnostics
MODULE diagnostics

  IMPLICIT NONE
  LOGICAL :: diagnostics_on

  PRIVATE
  PUBLIC enter_exit, get_diagnostics_on, set_diagnostics_on

CONTAINS

!!!######################################################################

  SUBROUTINE enter_exit(sub_name, state)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ENTER_EXIT" :: ENTER_EXIT
    USE other_consts, ONLY: MAX_SUBNAME_LEN
    IMPLICIT NONE

    INTEGER,INTENT(in) :: state
    CHARACTER(len=MAX_SUBNAME_LEN), INTENT(in) :: sub_name

    IF(diagnostics_on)THEN
       IF(state.EQ.1)THEN
          WRITE(*,'('' Entering subroutine '',60A,'':'')') sub_name(1:MAX_SUBNAME_LEN)
       ELSE
          WRITE(*,'('' Exiting subroutine '',60A,'':'')') sub_name(1:MAX_SUBNAME_LEN)
       ENDIF
    ENDIF

  END SUBROUTINE enter_exit

  SUBROUTINE set_diagnostics_on(state)
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"SO_SET_DIAGNOSTICS_ON":: SET_DIAGNOSTICS_ON
    IMPLICIT NONE

    LOGICAL, INTENT(in) :: state

    diagnostics_on = state

  END SUBROUTINE set_diagnostics_on

  SUBROUTINE get_diagnostics_on(state)
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"SO_GET_DIAGNOSTICS_ON":: GET_DIAGNOSTICS_ON
    IMPLICIT NONE

    LOGICAL :: state

    state = diagnostics_on

  END SUBROUTINE get_diagnostics_on

END MODULE diagnostics
