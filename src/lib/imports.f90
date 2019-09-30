MODULE imports
  !*Brief Description:* This module contains all the subroutines required to
  !import fields, previous model results, etc.
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  !
  !
  IMPLICIT NONE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  PRIVATE
  PUBLIC import_ventilation
  PUBLIC import_perfusion

CONTAINS
  !
  !##############################################################################
  !
  !>*import_ventilation:* This subroutine reads in the results of a ventilation model that
  ! has been saved in an exelem format as a single flow field (elements listed with
  ! ventilation as field values).
  SUBROUTINE import_ventilation(FLOWFILE)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_VENTILATION" :: IMPORT_VENTILATION
    USE arrays,ONLY: dp,elem_field,num_elems,num_units,units,unit_field,zero_tol,elem_cnct
    USE geometry,ONLY: get_final_real
    USE indices
    USE ventilation,ONLY: sum_elem_field_from_periphery
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    USE diagnostics, ONLY: enter_exit

    CHARACTER(len=MAX_FILENAME_LEN),INTENT(in) :: FLOWFILE
    !local variables
    INTEGER :: ierror,ne,nunit
    CHARACTER(LEN=132) :: ctemp1,exfile
    REAL(dp) :: flow,flow_unit,maxflow

    CHARACTER(len=60) :: sub_name

    sub_name = 'import_ventilation'
    CALL enter_exit(sub_name,1)

    PRINT *, 'Reading in ventilation results'
    CALL import_exelemfield(FLOWFILE,ne_Vdot)
    DO nunit = 1,num_units
       ne = units(nunit)
       IF(elem_field(ne_Vdot,ne).LT.0.0_dp) elem_field(ne_Vdot,ne) = zero_tol
       unit_field(nu_Vdot0,nunit) = elem_field(ne_Vdot,ne)
    ENDDO

!!! sum the fields up the tree
    CALL sum_elem_field_from_periphery(ne_Vdot) !sum the air flows recursively UP the tree
    maxflow = elem_field(ne_Vdot,1)


    CALL enter_exit(sub_name,2)
  END SUBROUTINE import_ventilation

  !
  !###########################################################################################
  !
  !>*import_perfusion:* This subroutine reads in the results of a ventilation model that
  ! has been saved in an exelem format as a single flow field (elements listed with
  ! ventilation as field values).
  SUBROUTINE import_perfusion(FLOWFILE)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_PERFUSION" :: IMPORT_PERFUSION
    USE arrays,ONLY: dp,elem_field,num_elems,num_units,units,unit_field,zero_tol,elem_cnct
    USE geometry,ONLY: get_final_real
    USE indices
    USE ventilation,ONLY: sum_elem_field_from_periphery
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    USE diagnostics, ONLY: enter_exit

    CHARACTER(len=MAX_FILENAME_LEN),INTENT(in) :: FLOWFILE
    !local variables
    INTEGER :: ierror,ne,nunit
    CHARACTER(LEN=132) :: ctemp1,exfile
    REAL(dp) :: flow,flow_unit,maxflow

    CHARACTER(len=60) :: sub_name

    sub_name = 'import_perfusion'
    CALL enter_exit(sub_name,1)

    PRINT *, 'Reading in perfusion results'
    CALL import_exelemfield(FLOWFILE,ne_Qdot)
    DO nunit = 1,num_units
       ne = units(nunit)
       IF(elem_field(ne_Qdot,ne).LT.0.0_dp) elem_field(ne_Qdot,ne) = zero_tol
       unit_field(nu_perf,nunit) = elem_field(ne_Qdot,ne)
    ENDDO

!!! sum the fields up the tree
    CALL sum_elem_field_from_periphery(ne_Qdot) !sum the air flows recursively UP the tree
    maxflow = elem_field(ne_Qdot,1)

    CALL enter_exit(sub_name,2)
  END SUBROUTINE import_perfusion

  !
  !##############################################################################
  !
  !>*import_exelemfield:* This subroutine reads in the content of an exelem field file (1 field)
  SUBROUTINE import_exelemfield(FLOWFILE,field_no)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_EXELEMFIELD" :: IMPORT_EXELEMFIELD
    USE arrays,ONLY: dp,elem_field,num_elems,num_units,units,unit_field,zero_tol,elem_cnct
    USE geometry,ONLY: get_final_real
    USE indices
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    USE diagnostics, ONLY: enter_exit

    CHARACTER(len=MAX_FILENAME_LEN),INTENT(in) :: FLOWFILE
    INTEGER, INTENT(in) :: field_no
    !local variables
    INTEGER :: ierror,ne,nunit
    CHARACTER(LEN=132) :: ctemp1,exfile
    REAL(dp) :: flow,flow_unit,maxflow

    CHARACTER(len=60) :: sub_name

    sub_name = 'import_exelemfield'
    CALL enter_exit(sub_name,1)

    OPEN(10, file=FLOWFILE, status='old')
    ne = 0
    read_elem_flow : DO !define a do loop name
       !.......read element flow
       READ(unit=10, fmt="(a)", iostat=ierror) ctemp1
       IF(INDEX(ctemp1, "Values:")> 0) THEN
          ne = ne+1
          READ(unit=10, fmt="(a)", iostat=ierror) ctemp1
          CALL get_final_real(ctemp1,flow)
          IF(flow.LT.0.0_dp) flow = zero_tol
          elem_field(field_no,ne) = flow! read it in
       END IF
       IF(ne.GE.num_elems) EXIT read_elem_flow
    END DO read_elem_flow

    CLOSE(10)

    CALL enter_exit(sub_name,2)
  END SUBROUTINE import_exelemfield

END MODULE imports
