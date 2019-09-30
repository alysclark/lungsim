MODULE field_utilities
  !*Brief Description:* This module contains all the subroutines that perform general operations
  !on fields
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  !More info on what the module does if necessary
  !
  USE other_consts
  IMPLICIT NONE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  PRIVATE
  PUBLIC scale_flow_to_inlet

CONTAINS

  !
  !*scale_flow_field* Scales a flow field to an 'inlet flow' value (real units).
  SUBROUTINE scale_flow_to_inlet(inlet_flow,VorQ)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SCALE_FLOW_TO_INLET" :: SCALE_FLOW_TO_INLET
    USE arrays,ONLY: dp,elem_field,num_elems,num_units,unit_field,zero_tol
    USE indices,ONLY: ne_dvdt,ne_Vdot,nu_Vdot0,ne_Qdot,nu_perf
    USE diagnostics, ONLY: enter_exit

    REAL(dp),INTENT(in) :: inlet_flow
    CHARACTER(len=1), INTENT(in) :: VorQ
    REAL(dp) :: ratio
    CHARACTER(len=60) :: sub_name

    sub_name = 'scale_flow_to_inlet'
    CALL enter_exit(sub_name,1)

    IF(VorQ.EQ.'V')THEN
       IF(ABS(elem_field(ne_Vdot,1)).GT.zero_tol)THEN
          ratio = inlet_flow/elem_field(ne_Vdot,1)
          unit_field(nu_Vdot0,1:num_units) = unit_field(nu_Vdot0,1:num_units)*ratio
          elem_field(ne_Vdot,1:num_elems) = elem_field(ne_Vdot,1:num_elems)*ratio
          elem_field(ne_dvdt,1:num_elems) = elem_field(ne_dvdt,1:num_elems)*ratio
       ELSE
          WRITE(*,'('' Cannot scale to zero flow'')')
       ENDIF
    ELSEIF(VorQ.EQ.'Q')THEN
       IF(ABS(elem_field(ne_Qdot,1)).GT.zero_tol)THEN
          ratio = inlet_flow/elem_field(ne_Qdot,1)
          unit_field(nu_perf,1:num_units) = unit_field(nu_perf,1:num_units)*ratio
          elem_field(ne_Qdot,1:num_elems) = elem_field(ne_Qdot,1:num_elems)*ratio
       ELSE
          WRITE(*,'('' Cannot scale to zero flow'')')
       ENDIF
    ELSE
       PRINT*, ' WARNING: Not a ventilation or perfusion scaling, no calculations done'
    ENDIF

    CALL enter_exit(sub_name,2)
  END SUBROUTINE scale_flow_to_inlet

END MODULE field_utilities
