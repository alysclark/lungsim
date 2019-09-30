!> \file
!> \author Merryn Tawhai
!> \brief This module handles all export functions
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>
!>\Description
!> This module handles all export functions
MODULE exports
  IMPLICIT NONE

  PRIVATE
  PUBLIC export_1d_elem_geometry,export_elem_geometry_2d,export_node_geometry,export_node_geometry_2d,&
       export_node_field,export_elem_field,export_terminal_solution,export_terminal_perfusion,&
       export_terminal_ssgexch,export_1d_elem_field,export_data_geometry

CONTAINS
!!!################################################################

  SUBROUTINE export_1d_elem_field(ne_field, EXELEMFILE, group_name, field_name )
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_1D_ELEM_FIELD" :: EXPORT_1D_ELEM_FIELD
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    USE arrays,ONLY: elem_field,num_elems
    IMPLICIT NONE

!!! Parameters
    INTEGER, INTENT(in) :: ne_field
    CHARACTER(len=MAX_FILENAME_LEN), INTENT(in) :: EXELEMFILE
    CHARACTER(len=MAX_STRING_LEN), INTENT(in) :: field_name
    CHARACTER(len=MAX_STRING_LEN), INTENT(in) :: group_name

!!! Local Variables
    INTEGER :: len_end,ne
    LOGICAL :: CHANGED

    OPEN(10, file=EXELEMFILE, status='replace')

    len_end=len_TRIM(group_name)
    !**     write the group name
    WRITE(10,'( '' Group name: '',A)') group_name(:len_end)
    !**         write the elements
    WRITE(10,'( '' Shape.  Dimension=1'' )')
    CHANGED=.TRUE. !initialise to force output of element information
    len_end=len_TRIM(field_name)
    DO ne=1,num_elems
       IF(ne>1) THEN
          CHANGED=.FALSE.
       ENDIF
       IF(CHANGED)THEN
          WRITE(10,'( '' #Scale factor sets=0'' )')
          WRITE(10,'( '' #Nodes= 0'' )')
          WRITE(10,'( '' #Fields= 1'' )')
          WRITE(10,'( '' 1)'',A,'', field, rectangular cartesian, #Components=1'')')&
               field_name(:len_end)
          WRITE(10,'( ''  '',A,''.  l.Lagrange, no modify, grid based.'')') &
               field_name(:len_end)
          WRITE(10,'( ''  #xi1=1'')')
       ENDIF

       WRITE(10,'(1X,''Element: '',I12,'' 0 0'' )') ne
       WRITE(10,'(3X,''Values:'' )')
       WRITE(10,'(4X,2(1X,E12.5))') elem_field(ne_field,ne),elem_field(ne_field,ne)
    ENDDO !no_nelist (ne)
    CLOSE(10)

  END SUBROUTINE export_1d_elem_field

!!!########################################################################

  SUBROUTINE export_1d_elem_geometry(EXELEMFILE, name)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_1D_ELEM_GEOMETRY" :: EXPORT_1D_ELEM_GEOMETRY

    USE arrays,ONLY: elem_nodes,num_elems
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    IMPLICIT NONE

!!! Parameters
    CHARACTER(len=MAX_FILENAME_LEN), INTENT(in) :: EXELEMFILE
    CHARACTER(len=MAX_STRING_LEN), INTENT(in) :: name

!!! Local Variables
    INTEGER :: len_end,ne,nj,nn
    CHARACTER(len=1) :: char1
    LOGICAL :: CHANGED

    OPEN(10, file=EXELEMFILE, status='replace')
    len_end=len_TRIM(name)
    !**     write the group name
    WRITE(10,'( '' Group name: '',A)') name(:len_end)
    !**         write the elements
    WRITE(10,'( '' Shape.  Dimension=1'' )')
    CHANGED=.TRUE. !initialise to force output of element information
    DO ne=1,num_elems
       IF(ne>1) THEN
          CHANGED=.FALSE.
       ENDIF
       IF(CHANGED)THEN
          WRITE(10,'( '' #Scale factor sets=1'' )')
          WRITE(10,'( ''   l.Lagrange, #Scale factors= 2'' )')
          WRITE(10,'( '' #Nodes= 2'' )')
          WRITE(10,'( '' #Fields= 1'' )')
          WRITE(10,'( '' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
          DO nj=1,3
             IF(nj==1) char1='x'; IF(nj==2) char1='y'; IF(nj==3) char1='z';
             WRITE(10,'(''  '',A2,''.  l.Lagrange, no modify, standard node based.'')') char1
             WRITE(10,'( ''     #Nodes= 2'')')
             DO nn=1,2
                WRITE(10,'(''      '',I1,''.  #Values=1'')') nn
                WRITE(10,'(''       Value indices:      1 '')')
                WRITE(10,'(''       Scale factor indices:'',I4)') nn
             ENDDO !nn
          ENDDO !nj
       ENDIF
       WRITE(10,'(1X,''Element: '',I12,'' 0 0'' )') ne
       !**               write the nodes
       WRITE(10,'(3X,''Nodes:'' )')
       WRITE(10,'(4X,2(1X,I12))') elem_nodes(1,ne),elem_nodes(2,ne)
       !**                 write the scale factors
       WRITE(10,'(3X,''Scale factors:'' )')
       WRITE(10,'(4X,2(1X,E12.5))') 1.d0,1.d0
    ENDDO !no_nelist (ne)
    CLOSE(10)

  END SUBROUTINE export_1d_elem_geometry

!!!############################################################################

  SUBROUTINE export_elem_geometry_2d(EXELEMFILE, name, offset_elem, offset_node)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_ELEM_GEOMETRY_2D" :: EXPORT_ELEM_GEOMETRY_2D

    USE arrays,ONLY: num_lines_2d,lines_2d,elem_versn_2d,elem_nodes_2d,nodes_2d,scale_factors_2d, &
         elem_lines_2d,num_elems_2d
    USE diagnostics,ONLY: enter_exit

!!! Parameters
    INTEGER :: offset_elem,offset_node
    CHARACTER(len=*) :: EXELEMFILE
    CHARACTER(len=*) :: name

!!! Local Variables
    INTEGER :: ne,nj,nk,nl,nn,nn_index(4),np_index(4),numnodes_ex,nvv(4)
    CHARACTER(len=1) :: char1
    CHARACTER(len=200) :: exfile
    LOGICAL :: CHANGED
    CHARACTER(len=60) :: sub_name = 'export_elem_geometry_2d'

    CALL enter_exit(sub_name,1)

    exfile = TRIM(exelemfile)//'.exelem'
    OPEN(10, file=exfile, status='replace')

    !**     write the group name
    WRITE(10,'( '' Group name: '',a)') TRIM(name)

    !**         write the lines
    IF(num_lines_2d.GT.0) THEN
       WRITE(10,'( '' Shape.  Dimension=1'' )')
       DO nl=1,num_lines_2d
          WRITE(10,'( '' Element: 0 0 '',I5)') lines_2d(nl)
       ENDDO !nl
    ENDIF

    !**         write the elements
    WRITE(10,'( '' Shape.  Dimension=2'' )')

    CHANGED=.TRUE. !initialise to force output of element information

    nvv=0

    DO ne=1,num_elems_2d
       IF(nvv(1)==elem_versn_2d(1,ne) .AND. nvv(2)==elem_versn_2d(2,ne) .AND. &
            nvv(3)==elem_versn_2d(3,ne) .AND. nvv(4)==elem_versn_2d(4,ne)) THEN
          CHANGED=.FALSE.
       ELSE
          CHANGED=.TRUE.
       ENDIF

       FORALL (nn=1:4) nvv(nn)=elem_versn_2d(nn,ne)
       numnodes_ex = 4
       FORALL (nn=1:4) nn_index(nn) = nn
       np_index(1:4) = elem_nodes_2d(1:4,ne)
       !       if(elem_nodes_2d(1,ne).eq.elem_nodes_2d(2,ne))then
       !          numnodes_ex = 3
       !          forall (nn=2:4) nn_index(nn) = nn-1
       !          np_index(2) = np_index(3)
       !          np_index(3) = np_index(4)
       !       elseif(elem_nodes_2d(1,ne).eq.elem_nodes_2d(3,ne))then
       !          numnodes_ex = 3
       !          nn_index(3) = 1
       !          nn_index(4) = 3
       !          np_index(3) = np_index(4)
       !       elseif(elem_nodes_2d(2,ne).eq.elem_nodes_2d(4,ne))then
       !          numnodes_ex = 3
       !          nn_index(4) = 2
       !       elseif(elem_nodes_2d(3,ne).eq.elem_nodes_2d(4,ne))then
       !          numnodes_ex = 3
       !          nn_index(4) = 3
       !       endif

       IF(CHANGED)THEN
          WRITE(10,'( '' #Scale factor sets=1'' )')
          WRITE(10,'( ''   c.Hermite*c.Hermite, #Scale factors=16'' )')
          WRITE(10,'( '' #Nodes= '',I2 )') numnodes_ex
          WRITE(10,'( '' #Fields= 1'' )')
          WRITE(10,'( '' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')

          DO nj=1,3
             IF(nj==1) char1='x'; IF(nj==2) char1='y'; IF(nj==3) char1='z';
             WRITE(10,'(''  '',A2,''.  c.Hermite*c.Hermite, no modify, standard node based.'')') char1
             WRITE(10,'( ''     #Nodes= 4'')')
             DO nn=1,4
                WRITE(10,'(''      '',I1,''.  #Values=4'')') nn_index(nn)
                WRITE(10,'(''       Value indices:       '',4I4)') 4*(nvv(nn)-1)+1, &
                     4*(nvv(nn)-1)+2,4*(nvv(nn)-1)+3,4*(nvv(nn)-1)+4
                WRITE(10,'(''       Scale factor indices:'',4I4)') 4*nn-3,4*nn-2,4*nn-1,4*nn
             ENDDO !nn
          ENDDO !nj
       ENDIF
       !**               write the element
       WRITE(10,'(1X,''Element: '',I12,'' 0 0'' )') ne+offset_elem
       !**                 write the faces
       WRITE(10,'(3X,''Faces: '' )')

       WRITE(10,'(5X,''0 0'',I6)')  elem_lines_2d(3,ne)
       WRITE(10,'(5X,''0 0'',I6)')  elem_lines_2d(4,ne)
       WRITE(10,'(5X,''0 0'',I6)')  elem_lines_2d(1,ne)
       WRITE(10,'(5X,''0 0'',I6)')  elem_lines_2d(2,ne)

       !**               write the nodes
       WRITE(10,'(3X,''Nodes:'' )')
       !       WRITE(10,'(4X,16(1X,I12))') (elem_nodes_2d(nn,ne)+offset_node,nn=1,4)
       WRITE(10,'(4X,16(1X,I12))') (nodes_2d(np_index(nn))+offset_node,nn=1,numnodes_ex)
       !**                 write the scale factors
       WRITE(10,'(3X,''Scale factors:'' )')
       WRITE(10,'(4X,4(1X,E12.5))') (scale_factors_2d(nk,ne),nk=1,4) !node 1
       WRITE(10,'(4X,4(1X,E12.5))') (scale_factors_2d(nk,ne),nk=5,8) !node 2
       WRITE(10,'(4X,4(1X,E12.5))') (scale_factors_2d(nk,ne),nk=9,12) !node 3
       WRITE(10,'(4X,4(1X,E12.5))') (scale_factors_2d(nk,ne),nk=13,16) !node 4
    ENDDO
    CLOSE(10)
    CALL enter_exit(sub_name,2)

  END SUBROUTINE export_elem_geometry_2d


!!!##########################################################################

  SUBROUTINE export_node_geometry(EXNODEFILE, name)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_NODE_GEOMETRY" :: EXPORT_NODE_GEOMETRY

    USE arrays,ONLY: node_xyz,num_nodes
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    IMPLICIT NONE

!!! Parameters
    CHARACTER(len=MAX_FILENAME_LEN),INTENT(in) :: EXNODEFILE
    CHARACTER(len=MAX_STRING_LEN),INTENT(in) :: name

!!! Local Variables
    INTEGER :: len_end,nj,np,np_last,VALUE_INDEX
    LOGICAL :: FIRST_NODE

    len_end=len_TRIM(name)
    IF(num_nodes.GT.0) THEN
       OPEN(10, file=EXNODEFILE, status='replace')
       !**     write the group name
       WRITE(10,'( '' Group name: '',A)') name(:len_end)
       FIRST_NODE=.TRUE.
       np_last=1
       !*** Exporting Geometry
       DO np=1,num_nodes
          IF(np.GT.1) np_last = np
          !*** Write the field information
          VALUE_INDEX=1
          IF(FIRST_NODE)THEN
             WRITE(10,'( '' #Fields=1'' )')
             WRITE(10,'('' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
             DO nj=1,3
                IF(nj.EQ.1) WRITE(10,'(2X,''x.  '')',advance="no")
                IF(nj.EQ.2) WRITE(10,'(2X,''y.  '')',advance="no")
                IF(nj.EQ.3) WRITE(10,'(2X,''z.  '')',advance="no")
                WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="no") 1,0
                WRITE(10,'()')
             ENDDO
          ENDIF !FIRST_NODE
          !***      write the node
          WRITE(10,'(1X,''Node: '',I12)') np
          DO nj=1,3
             WRITE(10,'(2X,4(1X,F12.6))') (node_xyz(nj,np))
          ENDDO !njj2
          FIRST_NODE=.FALSE.
          np_last=np
       ENDDO !nolist (np)
    ENDIF !num_nodes
    CLOSE(10)

  END SUBROUTINE export_node_geometry

!!!########################################################################

  SUBROUTINE export_node_geometry_2d(EXNODEFILE, name, offset)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_NODE_GEOMETRY_2D" :: EXPORT_NODE_GEOMETRY_2D

    USE arrays!,only: nodes_2d,node_xyz_2d,num_nodes_2d,node_versn_2d
    USE diagnostics, ONLY: enter_exit
    INTEGER :: offset
    CHARACTER(len=*) :: EXNODEFILE
    CHARACTER(len=*) :: name
    CHARACTER(len=60) :: sub_name = 'export_node_geometry_2d'


    !     Local Variables
    INTEGER :: nderiv,nversions,nj,nk,np,np_last,nv,VALUE_INDEX
    LOGICAL :: FIRST_NODE
    CHARACTER(len=200) :: exfile

    CALL enter_exit(sub_name,1)

    IF(num_nodes_2d.GT.0)THEN
       exfile = TRIM(exnodefile)//'.exnode'
       OPEN(10, file=exfile, status='replace')
       !**     write the group name
       WRITE(10,'( '' Group name: '',A)') TRIM(name)

       FIRST_NODE=.TRUE.
       np_last=1
       !*** Exporting Geometry
       DO np=1,num_nodes_2d
          IF(np.GT.1) np_last = np-1
          nderiv = 3
          nversions=node_versn_2d(np)
          !*** Write the field information
          VALUE_INDEX=1
          IF(FIRST_NODE.OR.node_versn_2d(np).NE.node_versn_2d(np_last))THEN
             WRITE(10,'( '' #Fields=1'' )')
             WRITE(10,'('' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
             DO nj=1,3
                IF(nj.EQ.1) WRITE(10,'(2X,''x.  '')',advance="no")
                IF(nj.EQ.2) WRITE(10,'(2X,''y.  '')',advance="no")
                IF(nj.EQ.3) WRITE(10,'(2X,''z.  '')',advance="no")
                IF(VALUE_INDEX<10)THEN
                   WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="no") VALUE_INDEX,NDERIV
                ELSE
                   WRITE(10,'(''Value index='',I2,'', #Derivatives='',I1)',advance="no") VALUE_INDEX,NDERIV
                ENDIF
                IF(nderiv.GE.1) WRITE(10,'('' (d/ds1'')',advance="no")
                IF(nderiv.GE.2) WRITE(10,'('',d/ds2'')',advance="no")
                IF(nderiv.GE.3) WRITE(10,'('',d2/ds1ds2'')',advance="no")

                IF(NVERSIONS.GT.1)THEN
                   WRITE(10,'(''),#Versions='',I2)') NVERSIONS
                ELSE IF(NDERIV.GT.0)THEN
                   WRITE(10,'('')'')')
                ELSE
                   WRITE(10,'()')
                ENDIF

                VALUE_INDEX=VALUE_INDEX+MAX(4*node_versn_2d(np),1)
             ENDDO

          ENDIF !FIRST_NODE
          !***      write the node
          WRITE(10,'(1X,''Node: '',I12)') nodes_2d(NP)+OFFSET
          DO nj=1,3
             IF(node_versn_2d(NP).GT.0) THEN
                DO nv=1,node_versn_2d(np)
                   WRITE(10,'(2X,4(1X,F12.6))') &
                        (node_xyz_2d(nk,nv,nj,NP),nk=1,4)
                ENDDO
             ELSE
                WRITE(10,'(3X,I1)') 0
             ENDIF
          ENDDO !njj2
          FIRST_NODE=.FALSE.
          np_last=np
       ENDDO !nolist (np)
    ENDIF
    CLOSE(10)

    CALL enter_exit(sub_name,2)

  END SUBROUTINE export_node_geometry_2d

!!!####################################################################

  SUBROUTINE export_data_geometry(EXDATAFILE, name, offset)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_DATA_GEOMETRY" :: EXPORT_DATA_GEOMETRY

    USE arrays, ONLY: num_data,data_xyz
    USE diagnostics,ONLY: enter_exit
!!! dummy arguments
    INTEGER :: offset
    CHARACTER(len=*) :: EXDATAFILE
    CHARACTER(len=*) :: name
!!! local variables
    INTEGER,PARAMETER :: num_coords = 3
    INTEGER nd,nj
    CHARACTER(len=200) :: exfile
    CHARACTER(len=60) :: sub_name = 'export_data_geometry'

    CALL enter_exit(sub_name,1)

    exfile = TRIM(exdatafile)//'.exdata'
    OPEN(10, file = exfile, status = 'replace')
    !**   write the group name
    WRITE(10,'( '' Group name: '',A)') TRIM(name)
    WRITE(10,'(1X,''#Fields=1'')')
    WRITE(10,'(1X,''1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
    WRITE(10,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
    WRITE(10,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
    WRITE(10,'(1X,''  z.  Value index= 3, #Derivatives=0'')')

    DO nd = 1,num_data
       WRITE(10,'(1X,''Node: '',I9)') nd + offset
       WRITE(10,'(1X,3E13.5)')  (data_xyz(nj,nd),nj=1,num_coords)
    ENDDO !NOLIST
    CLOSE(10)
    CALL enter_exit(sub_name,2)

  END SUBROUTINE export_data_geometry

!!!########################################################################

  SUBROUTINE export_terminal_solution(EXNODEFILE, name)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_TERMINAL_SOLUTION" :: EXPORT_TERMINAL_SOLUTION

    USE arrays,ONLY: elem_nodes,&
         node_xyz,num_units,units,unit_field
    USE indices,ONLY: nu_comp,nu_pe,nu_vt,nu_vent
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    IMPLICIT NONE

!!! Parameters
    CHARACTER(len=MAX_FILENAME_LEN),INTENT(in) :: EXNODEFILE
    CHARACTER(len=MAX_STRING_LEN),INTENT(in) :: name

!!! Local Variables
    INTEGER :: len_end,ne,nj,NOLIST,np,np_last,VALUE_INDEX
    LOGICAL :: FIRST_NODE

    len_end=len_TRIM(name)
    IF(num_units.GT.0) THEN
       OPEN(10, file=EXNODEFILE, status='replace')
       !**     write the group name
       WRITE(10,'( '' Group name: '',A)') name(:len_end)
       FIRST_NODE=.TRUE.
       np_last=1
       !*** Exporting Terminal Solution
       DO nolist=1,num_units
          IF(nolist.GT.1) np_last = np
          ne=units(nolist)
          np=elem_nodes(2,ne)
          !*** Write the field information
          VALUE_INDEX=1
          IF(FIRST_NODE)THEN
             WRITE(10,'( '' #Fields=5'' )')
             WRITE(10,'('' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
             DO nj=1,3
                IF(nj.EQ.1) WRITE(10,'(2X,''x.  '')',advance="no")
                IF(nj.EQ.2) WRITE(10,'(2X,''y.  '')',advance="no")
                IF(nj.EQ.3) WRITE(10,'(2X,''z.  '')',advance="no")
                WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
                VALUE_INDEX=VALUE_INDEX+1
             ENDDO
             !Ventilation (tidal volume/insp time)
             WRITE(10,'('' 2) flow, field, rectangular cartesian, #Components=1'')')
             WRITE(10,'(2X,''1.  '')',advance="no")
             WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
             !VALUE_INDEX=VALUE_INDEX+1
             !Volume
             !write(10,'('' 3) volume, field, rectangular cartesian, #Components=1'')')
             !write(10,'(2X,''1.  '')',advance="no")
             !write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
             !VALUE_INDEX=VALUE_INDEX+1
             !!Pressure
             !write(10,'('' 4) pressure, field, rectangular cartesian, #Components=1'')')
             !write(10,'(2X,''1.  '')',advance="no")
             !write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
             !Compliance
             WRITE(10,'('' 5) compliance, field, rectangular cartesian, #Components=1'')')
             WRITE(10,'(2X,''1.  '')',advance="no")
             WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
             VALUE_INDEX=VALUE_INDEX+1
             !Pleural pressure
             WRITE(10,'('' 6) pleural pressure, field, rectangular cartesian, #Components=1'')')
             WRITE(10,'(2X,''1.  '')',advance="no")
             WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
             VALUE_INDEX=VALUE_INDEX+1
             !Tidal volume
             WRITE(10,'('' 7) tidal volume, field, rectangular cartesian, #Components=1'')')
             WRITE(10,'(2X,''1.  '')',advance="no")
             WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
          ENDIF !FIRST_NODE
          !***      write the node
          WRITE(10,'(1X,''Node: '',I12)') np
          DO nj=1,3
             WRITE(10,'(2X,4(1X,F12.6))') (node_xyz(nj,np))      !Coordinates
          ENDDO !njj2
          WRITE(10,'(2X,4(1X,F12.6))') (unit_field(nu_vent,NOLIST)) !Ventilation
          !write(10,'(2X,4(1X,F12.6))') (unit_field(nu_vol,nolist))   !Volume (end expiration)
          !write(10,'(2X,4(1X,F12.6))') (unit_field(nu_press,nolist)) !Pressure
          WRITE(10,'(2X,4(1X,F12.6))') (unit_field(nu_comp,nolist))  !Compliance (end exp)
          WRITE(10,'(2X,4(1X,F12.6))') (unit_field(nu_pe,nolist))    !Recoil pressure
          WRITE(10,'(2X,4(1X,F12.6))') (unit_field(nu_vt,nolist))    !Tidal volume
          FIRST_NODE=.FALSE.
          np_last=np
       ENDDO !nolist (np)
    ENDIF !num_nodes
    CLOSE(10)

  END SUBROUTINE export_terminal_solution
!!! ##########################################################
  SUBROUTINE export_terminal_perfusion(EXNODEFILE, name)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_TERMINAL_PERFUSION" :: EXPORT_TERMINAL_PERFUSION

    USE arrays,ONLY: elem_nodes,&
         node_xyz,num_units,units,unit_field
    USE indices
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    IMPLICIT NONE

!!! Parameters
    CHARACTER(len=MAX_FILENAME_LEN),INTENT(in) :: EXNODEFILE
    CHARACTER(len=MAX_STRING_LEN),INTENT(in) :: name

!!! Local Variables
    INTEGER :: len_end,ne,nj,NOLIST,np,np_last,VALUE_INDEX
    LOGICAL :: FIRST_NODE

    len_end=len_TRIM(name)
    IF(num_units.GT.0) THEN
       OPEN(10, file=EXNODEFILE, status='replace')
       !**     write the group name
       WRITE(10,'( '' Group name: '',A)') name(:len_end)
       FIRST_NODE=.TRUE.
       np_last=1
       !*** Exporting Terminal Solution
       DO nolist=1,num_units
          IF(nolist.GT.1) np_last = np
          ne=units(nolist)
          np=elem_nodes(2,ne)
          !*** Write the field information
          VALUE_INDEX=1
          IF(FIRST_NODE)THEN
             WRITE(10,'( '' #Fields=3'' )')
             WRITE(10,'('' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
             DO nj=1,3
                IF(nj.EQ.1) WRITE(10,'(2X,''x.  '')',advance="no")
                IF(nj.EQ.2) WRITE(10,'(2X,''y.  '')',advance="no")
                IF(nj.EQ.3) WRITE(10,'(2X,''z.  '')',advance="no")
                WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
                VALUE_INDEX=VALUE_INDEX+1
             ENDDO
             !perfusion
             WRITE(10,'('' 2) flow, field, rectangular cartesian, #Components=1'')')
             WRITE(10,'(2X,''1.  '')',advance="no")
             WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
             !Pressure
             VALUE_INDEX=VALUE_INDEX+1
             WRITE(10,'('' 3) pressure, field, rectangular cartesian, #Components=1'')')
             WRITE(10,'(2X,''1.  '')',advance="no")
             WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
          ENDIF !FIRST_NODE
          !***      write the node
          WRITE(10,'(1X,''Node: '',I12)') np
          DO nj=1,3
             WRITE(10,'(2X,4(1X,F12.6))') (node_xyz(nj,np))      !Coordinates
          ENDDO !njj2
          WRITE(10,'(2X,4(1X,F12.6))') (unit_field(nu_perf,NOLIST)) !flow
          WRITE(10,'(2X,4(1X,F12.6))') (unit_field(nu_blood_press,NOLIST)) !pressure
          FIRST_NODE=.FALSE.
          np_last=np
       ENDDO !nolist (np)
    ENDIF !num_nodes
    CLOSE(10)

  END SUBROUTINE export_terminal_perfusion
!!!################################################
  SUBROUTINE export_terminal_ssgexch(EXNODEFILE, name)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_TERMINAL_SSGEXCH" :: EXPORT_TERMINAL_SSGEXCH

    USE arrays,ONLY: elem_nodes,&
         node_xyz,num_units,units,unit_field,gasex_field
    USE indices
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    IMPLICIT NONE

!!! Parameters
    CHARACTER(len=MAX_FILENAME_LEN),INTENT(in) :: EXNODEFILE
    CHARACTER(len=MAX_STRING_LEN),INTENT(in) :: name

!!! Local Variables
    INTEGER :: len_end,ne,nj,NOLIST,np,np_last,VALUE_INDEX
    LOGICAL :: FIRST_NODE

    len_end=len_TRIM(name)
    IF(num_units.GT.0) THEN
       OPEN(10, file=EXNODEFILE, status='replace')
       !**     write the group name
       WRITE(10,'( '' Group name: '',A)') name(:len_end)
       FIRST_NODE=.TRUE.
       np_last=1
       !*** Exporting Terminal Solution
       DO nolist=1,num_units
          IF(nolist.GT.1) np_last = np
          ne=units(nolist)
          np=elem_nodes(2,ne)
          !*** Write the field information
          VALUE_INDEX=1
          IF(FIRST_NODE)THEN
             WRITE(10,'( '' #Fields=5'' )')
             WRITE(10,'('' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
             DO nj=1,3
                IF(nj.EQ.1) WRITE(10,'(2X,''x.  '')',advance="no")
                IF(nj.EQ.2) WRITE(10,'(2X,''y.  '')',advance="no")
                IF(nj.EQ.3) WRITE(10,'(2X,''z.  '')',advance="no")
                WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
                VALUE_INDEX=VALUE_INDEX+1
             ENDDO
             !ventilation
             WRITE(10,'('' 2) alv_ventilation, field, rectangular cartesian, #Components=1'')')
             WRITE(10,'(2X,''1.  '')',advance="no")
             WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
             !perfusion
             VALUE_INDEX=VALUE_INDEX+1
             WRITE(10,'('' 3) cap_perfusion, field, rectangular cartesian, #Components=1'')')
             WRITE(10,'(2X,''1.  '')',advance="no")
             WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
             !pc_o2
             VALUE_INDEX=VALUE_INDEX+1
             WRITE(10,'('' 4) p_c_o2, field, rectangular cartesian, #Components=1'')')
             WRITE(10,'(2X,''1.  '')',advance="no")
             WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
             !pc_co2
             VALUE_INDEX=VALUE_INDEX+1
             WRITE(10,'('' 5) p_c_co2, field, rectangular cartesian, #Components=1'')')
             WRITE(10,'(2X,''1.  '')',advance="no")
             WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
          ENDIF !FIRST_NODE
          !***      write the node
          WRITE(10,'(1X,''Node: '',I12)') np
          DO nj=1,3
             WRITE(10,'(2X,4(1X,F12.6))') (node_xyz(nj,np))      !Coordinates
          ENDDO !njj2
          WRITE(10,'(2X,4(1X,F12.6))') (unit_field(nu_Vdot0,NOLIST)) !ventilation
          WRITE(10,'(2X,4(1X,F12.6))') (unit_field(nu_perf,NOLIST)) !perfusion
          WRITE(10,'(2X,4(1X,F12.6))') (gasex_field(ng_p_cap_o2,NOLIST)) !end capillary o2
          WRITE(10,'(2X,4(1X,F12.6))') (gasex_field(ng_p_cap_co2,NOLIST)) !end capillary co2
          FIRST_NODE=.FALSE.
          np_last=np
       ENDDO !nolist (np)
    ENDIF !num_nodes
    CLOSE(10)

  END SUBROUTINE export_terminal_ssgexch



!!! #################################################################

  SUBROUTINE export_node_field(nj_field, EXNODEFIELD, name, field_name)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_NODE_FIELD" :: EXPORT_NODE_FIELD

    USE arrays,ONLY: node_field,num_nodes
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    IMPLICIT NONE

!!! Parameters
    INTEGER,INTENT(in) :: nj_field
    CHARACTER(len=MAX_FILENAME_LEN),INTENT(in) :: EXNODEFIELD
    CHARACTER(len=MAX_STRING_LEN),INTENT(in) :: field_name
    CHARACTER(len=MAX_STRING_LEN),INTENT(in) :: name

!!! Local Variables
    INTEGER :: len_end,np
    LOGICAL :: FIRST_NODE

    OPEN(10, file=EXNODEFIELD, status='replace')
    !**     write the group name
    len_end=len_TRIM(name)
    WRITE(10,'( '' Group name: '',A)') name(:len_end)
    len_end=len_TRIM(field_name)
    FIRST_NODE=.TRUE.
    !*** the field as specified by user
    DO np=1,num_nodes
       !*** Write the field information
       IF(FIRST_NODE)THEN
          WRITE(10,'( '' #Fields=1'' )')
          WRITE(10,'('' 1) '',A,'', field, rectangular cartesian, #Components=1'')') &
               field_name(:len_end)
          WRITE(10,'(2X,''1.  '')',advance="no")
          WRITE(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") 1,0
       ENDIF !FIRST_NODE
       !***      write the node
       WRITE(10,'(1X,''Node: '',I12)') np
       WRITE(10,'(2X,2(1X,F12.6))') (node_field(nj_field,np))
       FIRST_NODE=.FALSE.
    ENDDO !num_nodes
    CLOSE(10)

  END SUBROUTINE export_node_field


!!! ###########################################################

  SUBROUTINE export_elem_field(EXELEMFIELD, name, field_name)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_ELEM_FIELD" :: EXPORT_ELEM_FIELD

    USE arrays,ONLY: elem_nodes,num_elems
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    IMPLICIT NONE

!!! Parameters
    CHARACTER(len=MAX_FILENAME_LEN), INTENT(in) :: EXELEMFIELD
    CHARACTER(len=MAX_STRING_LEN), INTENT(in) :: field_name
    CHARACTER(len=MAX_STRING_LEN), INTENT(in) :: name

!!! Local Variables
    INTEGER :: len_end,ne,nn
    LOGICAL :: CHANGED

    OPEN(10, file=EXELEMFIELD, status='replace')
    len_end=len_TRIM(name)
    !**     write the group name
    WRITE(10,'( '' Group name: '',A)') name(:len_end)
    !**         write the elements
    WRITE(10,'( '' Shape.  Dimension=1'' )')
    CHANGED=.TRUE. !initialise to force output of element information
    len_end=len_TRIM(field_name)
    DO ne=1,num_elems
       IF(ne>1) THEN
          CHANGED=.FALSE.
       ENDIF
       IF(CHANGED)THEN
          WRITE(10,'( '' #Scale factor sets=1'' )')
          WRITE(10,'( ''   l.Lagrange, #Scale factors= 2'' )')
          WRITE(10,'( '' #Nodes= 2'' )')
          WRITE(10,'( '' #Fields= 1'' )')
          WRITE(10,'( '' 1) '',A,'', field, rectangular cartesian, #Components=1'')') &
               field_name(:len_end)
          WRITE(10,'(''   1.  l.Lagrange, no modify, standard node based.'')')
          WRITE(10,'( ''     #Nodes= 2'')')
          DO nn=1,2
             WRITE(10,'(''      '',I1,''.  #Values=1'')') nn
             WRITE(10,'(''       Value indices:      1 '')')
             WRITE(10,'(''       Scale factor indices:'',I4)') nn
          ENDDO !nn
       ENDIF
       !**               write the element
       WRITE(10,'(1X,''Element: '',I12,'' 0 0'' )') ne
       !**               write the nodes
       WRITE(10,'(3X,''Nodes:'' )')
       WRITE(10,'(4X,2(1X,I12))') elem_nodes(1,ne),elem_nodes(2,ne)
       !**                 write the scale factors
       WRITE(10,'(3X,''Scale factors:'' )')
       WRITE(10,'(4X,2(1X,E12.5))') 1.d0,1.d0
    ENDDO !no_nelist (ne)
    CLOSE(10)

  END SUBROUTINE export_elem_field

END MODULE exports
