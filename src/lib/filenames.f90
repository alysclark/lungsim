!> \file
!> \author Merryn Tawhai, Alys Clark
!> \brief This module handles all read and write filenames.
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>

!> This module handles all read and write filenames.
MODULE filenames
  USE other_consts, ONLY : MAX_FILENAME_LEN
  IMPLICIT NONE
  !Module parameters

  CHARACTER(LEN=MAX_FILENAME_LEN) :: AIRWAY_ELEMFILE,AIRWAY_NODEFILE,AIRWAY_FIELDFILE, &
       AIRWAY_EXNODEFILE,AIRWAY_EXELEMFILE,AIRWAYFIELD_EXELEMFILE,&
       AIRWAY_MESHFILE, &
       ARTERY_ELEMFILE,ARTERY_NODEFILE,ARTERY_FIELDFILE, &
       ARTERY_EXNODEFILE,ARTERY_EXELEMFILE,ARTERYFIELD_EXELEMFILE,&
       ARTERY_MESHFILE, &
       TERMINAL_EXNODEFILE, MAIN_GEOMETRY_FILE, &
       FLOW_GEOMETRY_FILE,MAIN_PARAMETER_FILE,FLOW_PARAMETER_FILE,&
       FLOW_EXELEMFILE, FLOW_EXNODEFILE, FLOW_RADIUS_EXELEM, EMPTY_FILENAME

  PUBLIC AIRWAY_ELEMFILE,AIRWAY_NODEFILE,AIRWAY_FIELDFILE, &
       AIRWAY_EXNODEFILE,AIRWAY_EXELEMFILE,AIRWAYFIELD_EXELEMFILE,&
       AIRWAY_MESHFILE, &
       ARTERY_ELEMFILE,ARTERY_NODEFILE,ARTERY_FIELDFILE, &
       ARTERY_EXNODEFILE,ARTERY_EXELEMFILE,ARTERYFIELD_EXELEMFILE,&
       ARTERY_MESHFILE, &
       TERMINAL_EXNODEFILE, MAIN_GEOMETRY_FILE, &
       FLOW_GEOMETRY_FILE,MAIN_PARAMETER_FILE,FLOW_PARAMETER_FILE,&
       FLOW_EXELEMFILE, FLOW_EXNODEFILE, FLOW_RADIUS_EXELEM
  !Module types

  !Module variables

  !Interfaces

  PRIVATE
  PUBLIC read_geometry_main
  PUBLIC read_geometry_evaluate_flow
  PUBLIC get_filename

CONTAINS
  !
  !###################################################################################
  !
  !> reads in output filenames typically used to analyse and visualise ventilation model results
  SUBROUTINE read_geometry_evaluate_flow()
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_READ_GEOMETRY_EVALUATE_FLOW" :: READ_GEOMETRY_EVALUATE_FLOW

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    ! Input related variables
    CHARACTER(len=255) :: buffer, label
    INTEGER :: pos
    INTEGER, PARAMETER :: fh = 15
    INTEGER :: ios
    INTEGER :: line
    CHARACTER(len=60) :: sub_name

    sub_name = 'read_geometry_evaluate_flow'
    CALL enter_exit(sub_name,1)

    ios = 0
    line = 0

    OPEN(fh, file='Parameters/geometry_evaluate_flow.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    DO WHILE (ios == 0)
       READ(fh, '(A)', iostat=ios) buffer
       IF (ios == 0) THEN
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = SCAN(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          SELECT CASE (label)
          CASE ('exelem')
             READ(buffer, *, iostat=ios) AIRWAYFIELD_EXELEMFILE
          CASE ('exnode')
             READ(buffer, *, iostat=ios) TERMINAL_EXNODEFILE
          CASE ('flowexelem')
             READ(buffer, *, iostat=ios) FLOW_EXELEMFILE
          CASE ('flowexnode')
             READ(buffer, *, iostat=ios) FLOW_EXNODEFILE
          CASE ('flowradiusexelem')
             READ(buffer, *, iostat=ios) FLOW_RADIUS_EXELEM
          CASE default
             !print *, 'Skipping invalid label at line', line
          END SELECT
       END IF
    END DO

    CLOSE(fh)
    CALL enter_exit(sub_name,2)

  END SUBROUTINE read_geometry_evaluate_flow

  !###################################################################################

  SUBROUTINE read_geometry_main()
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_READ_GEOMETRY_MAIN" :: READ_GEOMETRY_MAIN

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    ! Input related variables
    CHARACTER(len=255) :: buffer, label
    INTEGER :: pos
    INTEGER, PARAMETER :: fh = 15
    INTEGER :: ios
    INTEGER :: line
    CHARACTER(len=60) :: sub_name
    !
    ! ###########################################################################
    !
    !> Reads in input and output filenames required to generate and export a geometry
    sub_name = 'read_geometry_main'
    CALL enter_exit(sub_name,1)

    ios = 0
    line = 0
    OPEN(fh, file='Parameters/geometry_main.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    DO WHILE (ios == 0)
       READ(fh, '(A)', iostat=ios) buffer
       IF (ios == 0) THEN
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = SCAN(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          SELECT CASE (label)
          CASE ('airway_ipnode')
             READ(buffer, *, iostat=ios) AIRWAY_NODEFILE
          CASE ('airway_ipelem')
             READ(buffer, *, iostat=ios) AIRWAY_ELEMFILE
          CASE ('airway_ipfiel')
             READ(buffer, *, iostat=ios) AIRWAY_FIELDFILE
          CASE ('airway_ipmesh')
             READ(buffer, *, iostat=ios) AIRWAY_MESHFILE
          CASE ('airway_exnode')
             READ(buffer, *, iostat=ios) AIRWAY_EXNODEFILE
          CASE ('airway_exelem')
             READ(buffer, *, iostat=ios) AIRWAY_EXELEMFILE
          CASE ('artery_ipnode')
             READ(buffer, *, iostat=ios) ARTERY_NODEFILE
          CASE ('artery_ipelem')
             READ(buffer, *, iostat=ios) ARTERY_ELEMFILE
          CASE ('artery_ipfiel')
             READ(buffer, *, iostat=ios) ARTERY_FIELDFILE
          CASE ('artery_ipmesh')
             READ(buffer, *, iostat=ios) ARTERY_MESHFILE
          CASE ('artery_exnode')
             READ(buffer, *, iostat=ios) ARTERY_EXNODEFILE
          CASE ('artery_exelem')
             READ(buffer, *, iostat=ios) ARTERY_EXELEMFILE
          CASE default
             !print *, 'Skipping invalid label at line', line
          END SELECT
       END IF
    END DO

    CLOSE(fh)
    CALL enter_exit(sub_name,2)

  END SUBROUTINE read_geometry_main
  !
  !#####################################################################################################
  !
  FUNCTION get_filename(label) RESULT(str)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GET_FILENAME" :: GET_FILENAME
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    IMPLICIT NONE

    CHARACTER(len=MAX_STRING_LEN), INTENT(in) :: label
    CHARACTER(len=MAX_FILENAME_LEN) :: str

    SELECT CASE (label)
    CASE ('airway_ipnode')
       str = AIRWAY_NODEFILE
    CASE ('airway_ipelem')
       str = AIRWAY_ELEMFILE
    CASE ('airway_ipfiel')
       str = AIRWAY_FIELDFILE
    CASE ('airway_ipmesh')
       str = AIRWAY_MESHFILE
    CASE ('airway_exnode')
       str = AIRWAY_EXNODEFILE
    CASE ('airway_exelem')
       str = AIRWAY_EXELEMFILE
    CASE ('artery_ipnode')
       str = ARTERY_NODEFILE
    CASE ('artery_ipelem')
       str = ARTERY_ELEMFILE
    CASE ('artery_ipfiel')
       str = ARTERY_FIELDFILE
    CASE ('artery_ipmesh')
       str = ARTERY_MESHFILE
    CASE ('artery_exnode')
       str = ARTERY_EXNODEFILE
    CASE ('artery_exelem')
       str = ARTERY_EXELEMFILE
    CASE ('exelem')
       str =  AIRWAYFIELD_EXELEMFILE
    CASE ('exnode')
       str =  TERMINAL_EXNODEFILE
    CASE ('flowexelem')
       str =  FLOW_EXELEMFILE
    CASE ('flowexnode')
       str =  FLOW_EXNODEFILE
    CASE ('flowradiusexelem')
       str =  FLOW_RADIUS_EXELEM
    CASE default
       PRINT *, 'Umm I dont know this label, sorry'
       str = EMPTY_FILENAME
    END SELECT

  END FUNCTION get_filename
END MODULE filenames
