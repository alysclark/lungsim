MODULE surface_fitting

  USE arrays
  USE other_consts
  USE mesh_functions
  !use precision
  USE solve

  IMPLICIT NONE

  PUBLIC fit_surface_geometry,pxi

  INTEGER,PARAMETER :: nmax_data_elem = 4000     ! max # data points on an element
  INTEGER,PARAMETER :: nmax_versn = 6            ! max # versions of node (derivative)
  INTEGER,PARAMETER :: num_coords = 3            ! x,y,z coordinates
  INTEGER,PARAMETER :: num_deriv = 4             ! coordinate + 3 derivatives
  INTEGER,PARAMETER :: num_deriv_elem = 16       ! nodes * num_deriv = 16
  INTEGER,PARAMETER :: num_elem_nodes = 4        ! code is unashamedly for 4-noded surface elements!
  INTEGER,PARAMETER :: num_fit = 3               ! number of parameters in fit (x,y,z)
  INTEGER,PARAMETER :: num_gauss = 9             ! number of gauss points

  INTEGER,PARAMETER :: nsize_gkk = 10000000        ! size of A matrix, until allocation sorted out

  !real(dp),parameter :: loose_tol = 1.0e-6_dp

CONTAINS

!!! ##########################################################################

  SUBROUTINE fit_surface_geometry(niterations,fitting_file)

!!! completes 'niterations' of geometry fitting to a surface, via minimising
!!! the least squares distance between a list of data points (3D RC coordinates)
!!! and a surface mesh (assumed bi-cubic Hermite only). 'fitting_file' lists
!!! the nodes/derivatives that are fixed, and any mapping of nodes and/or
!!! derivatives

    USE geometry,ONLY: get_local_node_f
!!! dummy arguments
    INTEGER,INTENT(in) :: niterations             ! user-specified number of fitting iterations
    CHARACTER(len=255),INTENT(in) :: fitting_file ! file that lists versions/mapping/BCs
!!! local variables
    INTEGER  :: nfit,nk,NOT_1,NOT_2,np,num_depvar,nv,ny_max
    LOGICAL :: first = .TRUE.
!!! local allocatable arrays
    INTEGER,ALLOCATABLE :: data_elem(:)
    INTEGER,ALLOCATABLE :: data_on_elem(:,:)     ! list of data closest to elements
    INTEGER,ALLOCATABLE :: elem_list(:)          ! list of elements in fit (all)
    INTEGER,ALLOCATABLE :: ndata_on_elem(:)      ! # of data points closest to element
    INTEGER,ALLOCATABLE :: npny(:,:)             ! maps deriv, vers, node, etc for a dep. variable
    INTEGER,ALLOCATABLE :: nynp(:,:,:,:)         ! dep. variable # for node, deriv, version etc.
    INTEGER,ALLOCATABLE :: nynr(:)               ! list of all dep. variables
    INTEGER,ALLOCATABLE :: nyny(:,:)             ! maps dep. variable to another dep. variable
    REAL(dp),ALLOCATABLE :: cyny(:,:)            ! weighting for mapped dep. variables
    REAL(dp),ALLOCATABLE :: data_xi(:,:)         ! xi location of data points
    REAL(dp),ALLOCATABLE :: fit_soln(:,:,:,:)    ! current fit geometry solution
    REAL(dp),ALLOCATABLE :: sobelov_wts(:,:)     ! Scaling factor for, and Sobelov weights
    LOGICAL,ALLOCATABLE :: fix_bcs(:)            ! logical for boundary conditions


!!! allocate element-sized arrays
    ALLOCATE(elem_list(0:num_elems_2d))
    ALLOCATE(sobelov_wts(0:6,num_elems_2d))

!!! allocate data arrays
    ALLOCATE(data_elem(num_data))
    ALLOCATE(data_on_elem(num_elems_2d,nmax_data_elem))
    ALLOCATE(ndata_on_elem(num_elems_2d))
    ALLOCATE(data_xi(2,num_data))
    data_elem = 0

!!! allocate dependent variable arrays
    ny_max = 0
    DO np = 1,num_nodes_2d
       DO nv = 1,node_versn_2d(np)
          ny_max = ny_max+1
       ENDDO
    ENDDO
    ny_max = ny_max*num_nodes_2d*num_fit*num_deriv  ! nodes * coordinates * #derivatives+1
    ALLOCATE(nynr(0:ny_max))
    ALLOCATE(npny(0:6,ny_max))
    ALLOCATE(nynp(num_deriv,nmax_versn,num_fit,num_nodes_2d))
    ALLOCATE(fix_bcs(ny_max))

!!! find the closest surface to each data point, and calculate the Xi
!!! coordinates of the data point to the surface
    WRITE(*,'('' Calculating normal projections: slow first time '')')
    CALL define_xi_closest(data_elem,data_on_elem,ndata_on_elem,data_xi,first)
    first = .FALSE.

!!! list the total error between the data points and the surface
    WRITE(*,'(/'' TOTAL RMS ERROR PRIOR TO FITTING:'')')
    CALL list_data_error(data_on_elem,ndata_on_elem,data_xi)

!!! read 'fitting_file' to define the fitting constraints. set up mapping
!!! arrays, dependent variable arrays, Sobelov smoothing weights (hardcoded)
    WRITE(*,'('' Define fitting problem '')')
    !    call define_geometry_fit(elem_list,npny,num_depvar,nynp,nynr,nyny,&
    !         cyny,sobelov_wts,fit_soln,fitting_file,fix_bcs)

    DO nfit = 1,niterations ! user-defined number of iterations
       CALL define_geometry_fit(elem_list,npny,num_depvar,nynp,nynr,nyny,&
            cyny,sobelov_wts,fit_soln,fitting_file,fix_bcs)
       WRITE(*,'(/'' FITTING ITERATION'',I3)') nfit
!!!    solve for new nodal coordinates and derivatives
       WRITE(*,'('' Solve fitting problem '')')
       CALL solve_geometry_fit(data_on_elem,ndata_on_elem,num_depvar,&
            elem_list,not_1,not_2,npny,nynp,nynr,&
            nyny,data_xi,cyny,sobelov_wts,fit_soln,fix_bcs)
!!!    update the scale factors for new geometry if NOT unit scale factors
       !       write(*,'('' Update scale factors '')')
       !       call update_scale_factor_norm
!!!    update the data point projections and their Xi coordinates
       WRITE(*,'('' Calculating normal projections '')')
       CALL define_xi_closest(data_elem,data_on_elem,ndata_on_elem,data_xi,first)
!!!    calculated the updated error between data and surface
       WRITE(*,'(/'' CURRENT RMS ERROR FOR ALL DATA:'')')
       CALL list_data_error(data_on_elem,ndata_on_elem,data_xi)
    ENDDO

    DEALLOCATE(elem_list)
    DEALLOCATE(sobelov_wts)
    DEALLOCATE(data_on_elem)
    DEALLOCATE(ndata_on_elem)
    DEALLOCATE(data_xi)
    DEALLOCATE(cyny)
    DEALLOCATE(nynr)
    DEALLOCATE(npny)
    DEALLOCATE(nynp)
    DEALLOCATE(nyny)
    DEALLOCATE(fix_bcs)

  END SUBROUTINE fit_surface_geometry

!!! ##########################################################################

  SUBROUTINE define_geometry_fit(elem_list,npny,num_depvar,nynp,nynr,nyny,&
       cyny,sobelov_wts,fit_soln,fitting_file,fix_bcs)

!!! read information from 'fitting_file' to determine the boundary conditions
!!! and mapping of dependent variables for geometric fitting. set up the
!!! dependent variable-to-mapping arrays

    USE arrays,ONLY: elems_2d,node_versn_2d,node_xyz_2d,num_elems_2d,num_nodes_2d
    USE geometry,ONLY: get_final_integer,get_local_node_f
!!! dummy arguments
    INTEGER :: elem_list(0:),npny(0:,:),num_depvar,nynp(:,:,:,:),nynr(0:)
    INTEGER,ALLOCATABLE :: nyny(:,:)
    REAL(dp) :: sobelov_wts(0:,:)
    REAL(dp),ALLOCATABLE :: cyny(:,:),fit_soln(:,:,:,:)
    CHARACTER(len=*),INTENT(in) :: fitting_file
    LOGICAL :: fix_bcs(:)
!!! local variables
    INTEGER :: i,ibeg,iend,ierror,IPFILE=10,i_ss_end,L,ne,nh,nj,nk,node,np,&
         np_global,number_of_fixed,nv,nv_fix,ny
    CHARACTER(len=132) :: readfile,string


    IF(.NOT.ALLOCATED(fit_soln)) ALLOCATE(fit_soln(4,10,16,num_nodes_2d))

    ! linear fitting for 3 geometric variables. solution stored in fields 1,2,3
    ! includes Sobelov smoothing on the geometry field

    !***Set up dependent variable interpolation information
    fit_soln = node_xyz_2d

!!! the following not correct because it refers to the global element #s
!!! use elem_list because we might want to fit only some of the elements
    !    elem_list(1:num_elems_2d) = elems_2d(1:num_elems_2d)
    FORALL (i=1:num_elems_2d) elem_list(i) = i
    elem_list(0) = num_elems_2d

    ! *** Specify smoothing constraints on each element
    DO L=1,elem_list(0)
       ne=elem_list(L)
       sobelov_wts(0,ne) = 1.0_dp
       sobelov_wts(1,ne) = 1.0_dp !the scaling factor for the Sobolev weights
       !  The 5 weights on derivs wrt Xi_1/_11/_2/_22/'_12 are:
       sobelov_wts(2,ne) = 1.0e-4_dp !weight for deriv wrt Xi_1
       sobelov_wts(3,ne) = 2.0e-3_dp
       sobelov_wts(4,ne) = 1.0e-4_dp
       sobelov_wts(5,ne) = 2.0e-3_dp
       sobelov_wts(6,ne) = 5.0e-3_dp
    ENDDO !L

    !*** Calculate ny maps
    CALL calculate_ny_maps(npny,num_depvar,nynp,nynr)

    fix_bcs = .FALSE. !initialise, default

    readfile = TRIM(fitting_file)//'.ipmap'
    OPEN(IPFILE, file = readfile, status='old')
    read_number_of_fixed : DO
       READ(unit=IPFILE, fmt="(a)", iostat=ierror) string
       IF(INDEX(string, "fixed")> 0) THEN
          CALL get_final_integer(string,number_of_fixed)
          EXIT read_number_of_fixed
       ENDIF
    END DO read_number_of_fixed

    DO node = 1,number_of_fixed
       READ(unit=IPFILE, fmt="(a)", iostat=ierror) string
       ibeg = 1
       i_ss_end = LEN(string) !get the end location of the sub-string
       iend=INDEX(string," ") !get location of next blank in sub-string
       READ (string(ibeg:iend-1), '(i6)' ) np_global
       np = get_local_node_f(2,np_global)

       string = ADJUSTL(string(iend:i_ss_end)) ! get chars beyond " " and remove the leading blanks
       iend=INDEX(string," ") !get location of next blank in sub-string
       READ (string(ibeg:iend-1), '(i6)' ) nv_fix

       string = ADJUSTL(string(iend:i_ss_end)) ! get chars beyond " " and remove the leading blanks
       READ (string(ibeg:i_ss_end), '(i6)' ) nk

       !       string = adjustl(string(iend:i_ss_end)) ! get chars beyond " " and remove the leading blanks
       !       read (string(ibeg:i_ss_end), '(i6)' ) nk

       nk=nk+1 !read in 0 for coordinate, 1 for 1st deriv, 2 for 2nd deriv
       IF(nv_fix.EQ.0)THEN ! do for all versions
          DO nv = 1,node_versn_2d(np)
             DO nh = 1,3
                ny = nynp(nk,nv,nh,np)
                fix_bcs(ny) = .TRUE.
             ENDDO !nh
          ENDDO !nv
       ELSE
          nv = nv_fix
          DO nh = 1,3
             ny = nynp(nk,nv,nh,np)
             fix_bcs(ny) = .TRUE.
             fit_soln(nk,nv,nh,np) = 0.0_dp
          ENDDO !nh
       ENDIF
    ENDDO !node

    CALL map_versions(IPFILE,num_depvar,nynp,nyny,cyny,fit_soln,fix_bcs)
    CLOSE(IPFILE)

    ! fix ALL of the cross derivatives, and set to zero
    DO np = 1,num_nodes_2d
       nk = 4 !index for 1-2 cross-derivative
       DO nv = 1,node_versn_2d(np)
          DO nj = 1,num_coords
             ny = nynp(nk,nv,nj,np)
             fix_bcs(ny) = .TRUE.
             node_xyz_2d(nk,nv,nj,np) = 0.0_dp
          ENDDO !nj
       ENDDO !nv
    ENDDO !np

  END SUBROUTINE define_geometry_fit

!!! ##########################################################################

  SUBROUTINE gauss1(PG)

!!! dummy arguments
    REAL(dp) :: PG(:,:,:)
!!! local variables
    INTEGER :: I,J,ng,nk,nn,ns,nu
    REAL(dp) :: D(3),XI(3),XIGG(3,3,2)

    D = [-0.3872983346207410_dp, 0.0_dp, 0.3872983346207410_dp]

    DO j=1,3
       DO i=1,3
          XIGG(i,j,1) = 0.50_dp+D(i)
          XIGG(i,j,2) = 0.50_dp+D(j)
          ng = I+(J-1)*3
          XI(1:2) = XIGG(i,j,1:2)
          ns=0
          DO nn=1,num_elem_nodes !number of local nodes
             DO nk=1,num_deriv !number of derivatives at node
                ns=ns+1
                DO nu=1,6 !number of xi coord derivs at node
                   PG(ns,nu,ng) = PSI1(nu,nk,nn,XI)
                ENDDO !nu
             ENDDO !nk
          ENDDO !nn
       ENDDO !i
    ENDDO !j

  END SUBROUTINE gauss1

!!! ##########################################################################

  FUNCTION getnyr(npny,ny,nynp)

!!! returns the dependent variable number
!!! dummy arguments
    INTEGER :: npny(0:,:),ny,nynp(:,:,:,:)
!!! local variables
    INTEGER :: nh,nk,np,nv
    INTEGER :: getnyr

    getnyr = 0
    nk = npny(1,ny)
    nv = npny(2,ny)
    nh = npny(3,ny)
    np = npny(4,ny)
    getnyr = nynp(nk,nv,nh,np)

  END FUNCTION getnyr

!!! ##########################################################################

  SUBROUTINE globalf(nony,not_1,not_2,npny,nyno,nynp,nyny,cony,cyno,cyny,fix_bcs)

!!! calculates the mapping arrays nyno/nony/cyno/cony

    USE arrays,ONLY: node_versn_2d,num_nodes_2d

!!! dummy arguments
    INTEGER :: nony(0:,:,:),not_1,not_2,npny(0:,:),nyno(0:,:,:),nynp(:,:,:,:),nyny(0:,:)
    REAL(dp) :: cony(0:,:,:),cyno(0:,:,:),cyny(0:,:)
    LOGICAL :: fix_bcs(:)
!!! local variables
    INTEGER :: nh,nv,nk,no,no_tot(2),np,nrc,ny,nyy(2),nyo,nyr,nyr2,nyy2(2),ny2
    REAL(dp) :: COY,RATIO
    LOGICAL :: done

!!!***  Initialise mapping arrays
    nony = 0
    cony = 0.0_dp
    nyno = 0
    cyno = 0.0_dp
    no_tot = 0

!!!*** Calculate mapping arrays
    DO np=1,num_nodes_2d
       DO nh=1,num_fit
          DO nv=1,node_versn_2d(np)
             DO nk=1,num_deriv
                ny=nynp(nk,nv,nh,np)
                IF(.NOT.fix_bcs(ny)) THEN
                   ! variable needs to be solved for
                   done = .FALSE.
                   nyy(1) = nynp(nk,nv,nh,np) !global row #
                   nyy(2) = ny !global variable #
                   ny2 = ny !the default
                   IF(nyny(0,ny).NE.0) THEN ! a special mapping
                      ny2 = nyny(1,ny)
                      RATIO = cyny(1,ny) ! the weighting of the mapping
                   ENDIF
                   IF(ny2.NE.ny) THEN ! for mapping, where ny exists
                      done = .TRUE.
                      IF(ny2.EQ.0) THEN !dof not used
                         fix_bcs(ny) = .TRUE.
                      ELSE IF(fix_bcs(ny2)) THEN
                         fix_bcs(ny) = .TRUE.
                      ELSE            ! no mapping
                         nyy2(1) = getnyr(npny,ny2,nynp) !row#
                         nyy2(2) = ny2 !global col#
                         DO nrc=1,2 !nrc=1,2 local row and local column
                            nyr = nyy(nrc)
                            nyr2 = nyy2(nrc)
                            nony(0,nyr,nrc) = 1
                            no = nony(1,nyr2,nrc)
                            nony(1,nyr,nrc) = no
                            COY = RATIO*cony(1,nyr2,nrc)
                            cony(1,nyr,nrc) = COY
                            !                            write(*,*) np,nh,nv,nk,ny,ny2,nrc
                            nyo = nyno(0,no,nrc)+1
                            nyno(0,no,nrc) = nyo
                            nyno(nyo,no,nrc) = nyr
                            cyno(nyo,no,nrc) = COY
                         ENDDO !nrc
                      ENDIF !ny2=0/fix_bcs
                   ENDIF !ny.NE.ny2

                   IF(.NOT.done) THEN
                      DO nrc=1,2 !rows and columns
                         no_tot(nrc) = no_tot(nrc)+1
                         nony(0,nyy(nrc),nrc) = 1
                         nony(1,nyy(nrc),nrc) = no_tot(nrc)
                         cony(0,nyy(nrc),nrc) = 0.0_dp
                         cony(1,nyy(nrc),nrc) = 1.0_dp
                         nyno(0,no_tot(nrc),nrc) = 1
                         nyno(1,no_tot(nrc),nrc) = nyy(nrc)
                         cyno(0,no_tot(nrc),nrc) = 0.0_dp
                         cyno(1,no_tot(nrc),nrc) = 1.0_dp
                      ENDDO !nrc
                   ENDIF !not done
                ENDIF !fix
             ENDDO !nk
          ENDDO !nv
       ENDDO !nh
    ENDDO !np

    NOT_1 = no_tot(1)
    NOT_2 = no_tot(2)

  END SUBROUTINE globalf

!!! ##########################################################################

  SUBROUTINE line_segments_for_2d_mesh

!!! sets up the line segment arrays for a 2d mesh

    USE arrays,ONLY: arclength,elem_cnct_2d,elem_nodes_2d,elem_versn_2d,elem_lines_2d,&
         lines_in_elem,line_versn_2d,&
         lines_2d,nodes_in_line,num_elems_2d,num_lines_2d,scale_factors_2d

!!! local variables
    INTEGER :: ne,ne_adjacent,ni1,nj,npn(2)
    LOGICAL :: MAKE

    num_lines_2d=0
    ! estimate number of lines, for allocating memory to arrays
    DO ne=1,num_elems_2d
       IF(elem_cnct_2d(-1,0,ne) == 0) num_lines_2d=num_lines_2d+1
       IF(elem_cnct_2d(-2,0,ne) == 0) num_lines_2d=num_lines_2d+1
       num_lines_2d=num_lines_2d+2 ! the minimum # of new lines for each element
    ENDDO

    IF(.NOT.ALLOCATED(lines_2d)) ALLOCATE (lines_2d(0:num_lines_2d))
    IF(.NOT.ALLOCATED(line_versn_2d)) ALLOCATE(line_versn_2d(2,3,num_lines_2d))
    IF(.NOT.ALLOCATED(elem_lines_2d)) ALLOCATE (elem_lines_2d(4,num_elems_2d))
    IF(.NOT.ALLOCATED(lines_in_elem)) ALLOCATE (lines_in_elem(0:4,num_lines_2d))
    IF(.NOT.ALLOCATED(nodes_in_line)) ALLOCATE (nodes_in_line(3,0:3,num_lines_2d))
    IF(.NOT.ALLOCATED(scale_factors_2d)) ALLOCATE(scale_factors_2d(16,num_elems_2d))
    IF(.NOT.ALLOCATED(arclength)) ALLOCATE(arclength(3,num_lines_2d))

    lines_in_elem=0
    lines_2d=0
    elem_lines_2d=0
    nodes_in_line=0
    line_versn_2d=0
    num_lines_2d=0

    DO ne=1,num_elems_2d
       !check whether to make a line
       MAKE=.FALSE.
       IF(elem_cnct_2d(-1,0,ne) == 0) MAKE=.TRUE. !exterior, make line
       ne_adjacent=elem_cnct_2d(-1,1,ne)
       IF(ne_adjacent.GT.0)THEN
          IF(elem_lines_2d(4,ne_adjacent) == 0) MAKE=.TRUE.
       ENDIF

       IF(MAKE)THEN
          num_lines_2d=num_lines_2d+1
          lines_2d(num_lines_2d)=num_lines_2d !record a new line number
          lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
          lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
          elem_lines_2d(3,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 3 of ne
          npn(1)=1
          npn(2)=3
          nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(1,ne) !records 1st node in line
          nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(3,ne) !records 2nd node in line
          nodes_in_line(1,0,num_lines_2d)=2 !Xi-direction of line segment num_lines_2d
          DO nj=1,3
             nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
             DO ni1=1,2
                line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
             ENDDO !n
          ENDDO !nj
       ELSE !get adjacent element line number
          !WARNING:: this only works if all Xi directions are consistent!!!!
          ne_adjacent=elem_cnct_2d(-1,1,ne)
          elem_lines_2d(3,ne)=elem_lines_2d(4,ne_adjacent)
       ENDIF

       !check whether to make a line
       MAKE=.FALSE.
       IF(elem_cnct_2d(-2,0,ne) == 0) MAKE=.TRUE. !exterior, make line
       ne_adjacent=elem_cnct_2d(-2,1,ne)
       IF(ne_adjacent.GT.0)THEN
          IF(elem_lines_2d(2,ne_adjacent) == 0) MAKE=.TRUE.
       ENDIF

       IF(MAKE)THEN
          num_lines_2d=num_lines_2d+1
          lines_2d(num_lines_2d)=num_lines_2d !record a new line number
          lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
          lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
          elem_lines_2d(1,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 1 of ne
          npn(1)=1
          npn(2)=2
          nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(1,ne) !records 1st node in line
          nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(2,ne) !records 2nd node in line
          !        write(*,*) 'line in -2 for ne',ne,num_lines_2d,' nodes',npne(1,ne),npne(2,ne)
          nodes_in_line(1,0,num_lines_2d)=1 !Xi-direction of line segment num_lines_2d
          DO nj=1,3
             nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
             DO ni1=1,2
                line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
             ENDDO !n
          ENDDO !nj
       ELSE !get adjacent element line number
          !WARNING:: this only works if all Xi directions are consistent!!!!
          ne_adjacent=elem_cnct_2d(-2,1,ne)
          elem_lines_2d(1,ne)=elem_lines_2d(2,ne_adjacent)
          !        write(*,*) 'adjacent in -2',ne,ne_adjacent,elem_lines_2d(1,ne)
       ENDIF

       num_lines_2d=num_lines_2d+1
       lines_2d(num_lines_2d)=num_lines_2d !record a new line number
       lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
       lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
       elem_lines_2d(4,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 4 of ne
       npn(1)=2
       npn(2)=4
       nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(2,ne) !records 1st node in line
       nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(4,ne) !records 2nd node in line
       !     write(*,*) 'line in +2 for ne',ne,num_lines_2d,' nodes',npne(2,ne),npne(4,ne)
       nodes_in_line(1,0,num_lines_2d)=2 !Xi-direction of line segment num_lines_2d
       DO nj=1,3
          nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
          DO ni1=1,2
             line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
          ENDDO !n
       ENDDO !nj

       num_lines_2d = num_lines_2d+1
       lines_2d(num_lines_2d) = num_lines_2d !record a new line number
       lines_in_elem(0,num_lines_2d) = lines_in_elem(0,num_lines_2d)+1
       lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d) = ne !line num_lines_2d is in element ne
       elem_lines_2d(2,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 2 of ne
       npn(1) = 3
       npn(2) = 4
       nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(3,ne) !records 1st node in line
       nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(4,ne) !records 2nd node in line
       nodes_in_line(1,0,num_lines_2d)=1 !Xi-direction of line segment num_lines_2d
       DO nj=1,3
          nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
          DO ni1=1,2
             line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
          ENDDO !n
       ENDDO !nj

    ENDDO !ne

    CALL calc_scale_factors_2d('arcl')

  END SUBROUTINE line_segments_for_2d_mesh

!!! ##########################################################################

  SUBROUTINE list_data_error(data_on_elem,ndata_on_elem,data_xi)

!!! calculate and write out the RMS error for distance between data points
!!! and 2d mesh surface

    USE arrays,ONLY: data_xyz,num_elems_2d
!!! dummy arguments
    INTEGER :: data_on_elem(:,:),ndata_on_elem(:)
    REAL(dp) :: data_xi(:,:)
!!! local variables
    INTEGER elem,nd,nde,num_data_infit,ne,nj
    REAL(dp) :: data_xi_local(2),EDD,SAED,SMED,SUM,SQED,X(6),&
         XE(num_deriv_elem,num_coords)


    SMED=0.0_dp
    SAED=0.0_dp
    SQED=0.0_dp
    num_data_infit=0

    DO ne=1,num_elems_2d
       CALL xpxe(ne,xe)
       elem=ne
       DO nde=1,ndata_on_elem(elem) !for each data point on element
          nd=data_on_elem(elem,nde) !the data point number
          data_xi_local(1:2) = data_xi(1:2,nd)
          DO nj=1,num_coords
             X(nj)=PXI(1,data_xi_local,XE(1,nj))
          ENDDO
          SUM=0.0_dp
          DO nj=1,num_coords
             SUM=SUM+(X(nj)-data_xyz(nj,nd))**2
          ENDDO !nj
          EDD=DSQRT(SUM)
          SMED=SMED+EDD
          SAED=SAED+DABS(EDD)
          SQED=SQED+EDD**2
          num_data_infit=num_data_infit+1
       ENDDO !nde
    ENDDO !list of elements

    IF(num_data_infit.GT.1) THEN
       WRITE(*,'('' Number of data points in fit ='',I8)') num_data_infit
       !     write(*,'('' Average error           : '',D12.6,'' +/- '',D12.6)') &
       !          SMED/DBLE(num_data_infit), &
       !          DSQRT((SQED-SMED**2/DBLE(num_data_infit))/DBLE(num_data_infit-1))

       WRITE(*,'('' Average absolute error  : '',D12.6,'' +/- '',D12.6)') &
            SAED/DBLE(num_data_infit),DSQRT((SQED-SAED**2/DBLE(num_data_infit))/ &
            DBLE(num_data_infit-1))
       WRITE(*,'('' Root mean squared error : '',D12.6)') &
            DSQRT(SQED/DBLE(num_data_infit))
    ELSE
       WRITE(*,'('' No data points in any elements'')')
       STOP
    ENDIF !ndtot>1

  END SUBROUTINE list_data_error

!!! ##########################################################################

  SUBROUTINE map_versions(IPFILE,num_depvar,nynp,nyny,cyny,fit_soln,fix_bcs)

    USE arrays,ONLY: node_versn_2d,node_xyz_2d,num_nodes_2d
    USE geometry,ONLY: get_final_integer,get_local_node_f
!!! dummy arguments
    INTEGER, INTENT(in) :: IPFILE,num_depvar
    INTEGER :: nynp(:,:,:,:)
    INTEGER,ALLOCATABLE :: nyny(:,:)
    REAL(dp),ALLOCATABLE :: cyny(:,:)
    REAL(dp) :: fit_soln(:,:,:,:)
    LOGICAL :: fix_bcs(:)
!!! local variables
    INTEGER :: i,ibeg,iend,ierror,i_ss_end,nj,node,np,number_of_maps,nv, &
         NV_MAX,ny,nk_t,nv_t,nj_t,np_t,nk_m,nv_m,nj_m,np_m, &
         nmap_info(100,7),ny_t
    REAL(dp) :: r_map_coef
    CHARACTER(len=132) :: string

!!! fix the boundary conditions for coordinates for nodes with versions, such that
!!! versions higher than 1 map to version 1
    DO np=1,num_nodes_2d
       NV_MAX=node_versn_2d(np)
       IF(NV_MAX>1)THEN
          DO nv=2,NV_MAX
             DO nj=1,num_coords
                node_xyz_2d(1,nv,nj,np) = node_xyz_2d(1,1,nj,np)
                fit_soln(1,nv,nj,np) = fit_soln(1,1,nj,np)
                ny = nynp(1,nv,nj,np)
                fix_bcs(ny) = .TRUE.
             ENDDO !nj
          ENDDO !nv
       ENDIF
    ENDDO !node

!!! read in the following for mapping:
!!!    node, version, derivative >> node, version, derivative, mapping coefficient
!!!    default is that all versions are independent

    read_number_of_mappings : DO
       READ(unit=IPFILE, fmt="(a)", iostat=ierror) string
       ! read line containing "Number of mappings"
       IF(INDEX(string, "mappings")> 0) THEN
          CALL get_final_integer(string,number_of_maps)
          EXIT read_number_of_mappings
       ENDIF
    END DO read_number_of_mappings

    ! allocate memory for dependent variable mapping arrays
    WRITE(*,*) 'Number of dependent variables =',num_depvar,'; squared =',num_depvar**2
    IF(.NOT.ALLOCATED(cyny)) ALLOCATE(cyny(0:number_of_maps,num_depvar))
    IF(.NOT.ALLOCATED(nyny)) ALLOCATE(nyny(0:number_of_maps,num_depvar))
    nyny = 0       ! initialise depvar to depvar mapping
    cyny = 0.0_dp  ! initialise weighting for mappings

!!! note that the global node numbers are used in the mapping file, whereas we need to use
!!! local numbering for the computation. Read in as global and then map to local below.
    DO node=1,number_of_maps ! for the number of nodes with mappings
       READ(unit=IPFILE, fmt="(a)", iostat=ierror) string
       ibeg=1
       i_ss_end=LEN(string) !get the end location of the sub-string
       DO i=1,6
          iend=INDEX(string," ") !get location of next blank in sub-string
          READ (string(ibeg:iend-1), '(i6)' ) nmap_info(node,i)
          string = ADJUSTL(string(iend:i_ss_end)) ! get the characters beyond " " and remove the leading blanks
       ENDDO !i
       READ (string(ibeg:i_ss_end), '(i6)' ) nmap_info(node,7)
    ENDDO !node

    DO node = 1,number_of_maps !for each mapping
       DO nj = 1,num_coords
          nk_m = nmap_info(node,3)+1 !derivative
          nv_m = nmap_info(node,2) !version
          nj_m = nj !coordinate
          np_m = get_local_node_f(2,nmap_info(node,1)) !global node mapped to local node
          ny = nynp(nk_m,nv_m,nj_m,np_m)
          nk_t = nmap_info(node,6)+1 !derivative
          nv_t = nmap_info(node,5) !version
          nj_t = nj !coordinate
          np_t = get_local_node_f(2,nmap_info(node,4)) !global node mapped to local node
          ny_t = nynp(nk_t,nv_t,nj_t,np_t)
          r_map_coef = REAL(nmap_info(node,7)) !mapping coefficient, +1 or -1

          IF(ny > 0) THEN
             nyny(0,ny) = nyny(0,ny)+1 ! increment array size
             nyny(nyny(0,ny),ny) = ny_t
             cyny(0,ny) = 0.0_dp
             cyny(nyny(0,ny),ny) = r_map_coef
             node_xyz_2d(nk_m,nv_m,nj_m,np_m) = node_xyz_2d(nk_t,nv_t,nj_t,np_t)*r_map_coef
             fit_soln(nk_m,nv_m,nj_m,np_m) = node_xyz_2d(nk_t,nv_t,nj_t,np_t)*r_map_coef
             !             write(*,*) 'mapping ny',ny_t,' to',ny,' with',r_map_coef
          ENDIF ! ny.GT.0
       ENDDO !nj
    ENDDO

!!! make all coordinates for different versions be the same
    DO np=1,num_nodes_2d
       DO nj=1,3
          ny_t=nynp(1,1,nj,np)
          DO nv=2,node_versn_2d(np) !for each version
             ny=nynp(1,nv,nj,np)
             IF(ny > 0) THEN
                nyny(0,ny)=nyny(0,ny)+1 ! increment array size
                nyny(nyny(0,ny),ny)=ny_t
                cyny(0,ny)=0.0_dp
                cyny(nyny(0,ny),ny)=1.0_dp
                node_xyz_2d(1,nv,nj,np)=node_xyz_2d(1,1,nj,np)
                fit_soln(1,nv,nj,np)=node_xyz_2d(1,1,nj,np)
             ENDIF ! ny.GT.0
          ENDDO !nv
       ENDDO !nj
    ENDDO

  END SUBROUTINE map_versions

!!! ##########################################################################

  SUBROUTINE melgef(LGE2,ne,NHST,nynp)

!!! calculates the row numbers (LGE(*,1)) and column numbers
!!! (LGE(*,2)) in the matrix for fitting for element variables nhs
!!! and fit variable njj in region nr.  It also returns the total
!!! number of element variables NHST(nrc).

    USE arrays,ONLY: elem_nodes_2d,elem_versn_2d
!!! dummy arguments
    INTEGER :: LGE2(num_fit*num_deriv_elem,2),ne,NHST(2),nynp(:,:,:,:)
!!! local variables
    INTEGER nh,nk,nn,np,nrc,nv

    DO nrc=1,2
       NHST(nrc)=0
       DO nh=1,num_fit
          DO nn=1,num_elem_nodes !nodal variables
             np=elem_nodes_2d(nn,ne)
             nv=elem_versn_2d(nn,ne)
             DO nk=1,num_deriv
                NHST(nrc)=NHST(nrc)+1
                LGE2(NHST(nrc),nrc)=nynp(nk,nv,nh,np)
             ENDDO !nk
          ENDDO !nn
       ENDDO !nhj
    ENDDO !nrc

  END SUBROUTINE melgef

!!! ##########################################################################

  FUNCTION psi1(nu,nk,nn,XI)

!!! dummy arguments
    INTEGER nk,nn,nu
    REAL(dp) :: XI(:)
!!! local variables
    INTEGER :: ido(num_deriv,2),inp(4,2),ipu(6,2),ni
    REAL(dp) :: psi1

    ipu = RESHAPE([1,2,3,1,1,2,1,1,1,2,3,2],SHAPE(ipu))
    inp = RESHAPE([1,2,1,2,1,1,2,2],SHAPE(inp)) ! the node index positions
    ido = RESHAPE([1,2,1,2,1,1,2,2],SHAPE(ido))

    psi1 = 1.0_dp
    DO ni=1,2
       psi1 = psi1*ph3(inp(nn,ni),ido(nk,ni),ipu(nu,ni),xi(ni))
    ENDDO

  END FUNCTION psi1

!!! ###############################################################

  FUNCTION pxi(nu,xi,x)

!!! dummy arguments
    INTEGER :: nu
    REAL(dp) :: xi(2),x(num_deriv_elem)
!!! local variables
    INTEGER :: nk,nn,ns
    REAL(dp) :: pxi

    pxi = 0.0_dp
    ns = 0
    DO nn=1,num_elem_nodes
       DO nk=1,num_deriv
          ns = ns+1
          pxi = pxi + psi1(nu,nk,nn,xi)*x(ns)
       ENDDO
    ENDDO

  END FUNCTION pxi

!!! ##########################################################################

  SUBROUTINE update_scale_factor_norm

    USE arrays,ONLY: node_versn_2d,node_xyz_2d,num_nodes_2d

!!! local variables
    INTEGER :: nj,nk,nk1,np,nv
    REAL(dp) :: SCALE,XD(3),ZERO_TOL=1.0e-12_dp

    DO np=1,num_nodes_2d
       DO nv=1,node_versn_2d(np)
          DO nk1=1,2 !loop over 1st derivs
             nk = nk1+1
             XD(1:3) = node_xyz_2d(nk,nv,1:3,np)
             SCALE = XD(1)**2+XD(2)**2+XD(3)**2
             IF(SCALE.GT.ZERO_TOL) THEN
                ! Normalise the nodal derivatives
                SCALE = 1.0_dp/DSQRT(SCALE)
                DO nj=1,3
                   node_xyz_2d(nk,nv,nj,np) = node_xyz_2d(nk,nv,nj,np)*SCALE
                   ! adjust cross-derivatives
                   node_xyz_2d(4,nv,nj,np) = node_xyz_2d(4,nv,nj,np)*SCALE
                ENDDO ! nj
             ENDIF !>0
          ENDDO ! nk1
       ENDDO !nv
    ENDDO !no_np

    CALL calc_scale_factors_2d('arcl')

  END SUBROUTINE update_scale_factor_norm

!!! ##########################################################################

  SUBROUTINE xpxe(ne,xe)

!!! copies geometry information from nodes into a local element array

    USE arrays,ONLY: elem_nodes_2d,elem_versn_2d,node_xyz_2d,scale_factors_2d
!!! dummy arguments
    INTEGER,INTENT(in) :: ne
    REAL(dp) :: xe(:,:)
!!! local variablesK     Local Variables
    INTEGER :: nj,nk,nn,np,ns,nv

    DO nj=1,3
       ns=0
       DO nn=1,num_elem_nodes
          np=elem_nodes_2d(nn,ne)
          nv=elem_versn_2d(nn,ne)
          DO nk=1,num_deriv
             ns=ns+1
             xe(ns,nj)=node_xyz_2d(nk,nv,nj,np)*scale_factors_2d(ns,ne)
          ENDDO
       ENDDO
    ENDDO !nj

  END SUBROUTINE xpxe

!!! ##########################################################################

  SUBROUTINE zder(data_on_elem,ndata_on_elem,ne,data_xi,ER,PG,WDL,WG,sobelov_wts,&
       XIDL,fit_soln_local)

!!!    Evaluates element rhs, ER(ns), in calculation of least squares
!!!    fit of linear field variables, defined by nodal values
!!!    node_xyz_2d(nk,nv,nj,np), to the set of data values data_xyz(nj,nd) with
!!!    weights data_weight(nj,nd) at local coordinate values data_xi(ni,nd).

    USE arrays,ONLY: data_weight,data_xyz,scale_factors_2d
!!! dummy arguments
    INTEGER :: data_on_elem(:,:),ndata_on_elem(:),ne
    REAL(dp) :: data_xi(:,:),ER(:),PG(:,:,:),WDL(:,:),WG(:),sobelov_wts(0:,:),&
         XIDL(:,:),fit_soln_local(:,:)
!!! local variables
    INTEGER nd,nde,ng,nh,nhs1,nk1,nn1,ns1,ns2,nu
    REAL(dp) :: SUM1,SUM2,SUM3,SUM4,X,ZDL(3,nmax_data_elem)

    DO nde = 1,ndata_on_elem(ne) ! for each data point on the element
       nd = data_on_elem(ne,nde) ! the data point number
       XIDL(1:2,nde) = data_xi(1:2,nd)
       ZDL(1:3,nde) = data_xyz(1:3,nd)
       WDL(1:3,nde) = data_weight(1:3,nd)
    ENDDO !nde

    nhs1=0
    DO nh = 1,num_fit
       DO nde = 1,ndata_on_elem(ne)
          X = PXI(1,XIDL(1:2,nde),fit_soln_local(1:num_deriv_elem,nh))
          ZDL(nh,nde) = ZDL(nh,nde)-X
       ENDDO !nde
       ns1 = 0
       DO nn1 = 1,num_elem_nodes
          DO nk1 = 1,num_deriv
             nhs1 = nhs1+1
             ns1 = ns1+1
             SUM1 = 0.0_dp
             DO nde = 1,ndata_on_elem(ne)
                SUM1 = SUM1+PSI1(1,nk1,nn1,XIDL(1:2,nde))*ZDL(nh,nde)*WDL(nh,nde)
             ENDDO !nde
             SUM2 = 0.0_dp
             DO ng = 1,num_gauss
                SUM3 = 0.0_dp
                DO nu = 2,6 !for 2d elements
                   SUM4 = 0.0_dp
                   DO ns2 = 1,num_deriv_elem
                      SUM4 = SUM4+fit_soln_local(ns2,nh)*PG(ns2,nu,ng)
                   ENDDO !ns2
                   SUM3 = SUM3+SUM4*PG(ns1,nu,ng)*sobelov_wts(nu,ne) !*sobelov_wts(1,ne)
                ENDDO !nu
                SUM2 = SUM2-SUM3*WG(ng) !*RG(ng)
             ENDDO !ng
             ER(nhs1) = ER(nhs1)+(SUM1+SUM2*sobelov_wts(0,ne))*scale_factors_2d(ns1,ne)
          ENDDO !nk1
       ENDDO !nn1
    ENDDO !nhj1

  END SUBROUTINE zder

!!! ##########################################################################

  SUBROUTINE zdes(ndata_on_elem,ne,ES,PG,WDL,WG,sobelov_wts,XIDL)

!!!    ZDES evaluates element stiffness matrix ES(ms,ns) in calculation
!!!    of least squares fit of linear field variables, defined by nodal
!!!    values node_xyz_2d(nk,nv,nj,np), to the set of data values XD(nj,nd) with
!!!    weights data_weight(nj,nd) at local coordinate values data_xi(ni,nd), where
!!!    nj=NJO.

    USE arrays,ONLY: scale_factors_2d
!!! dummy arguments
    INTEGER :: ndata_on_elem(:),ne
    REAL(dp) :: ES(:,:),PG(:,:,:),WDL(:,:),WG(:),sobelov_wts(0:,:),XIDL(:,:)
!!! local variables
    INTEGER nde,ng,nh1,nh2,nhj1,nhj2,nhs1,nhs1_for_nhj1,nhs2, &
         nk1,nk2,nn1,nn2,ns1,ns2,nu
    REAL(dp) :: PD(num_deriv_elem),SUM2,SUM3

    ES = 0.0_dp
    nhs1 = 0
    ! for each of the 3 dependent variables to be fitted (num_fit(1)=3)
    DO nhj1=1,num_fit !nhj are vars for the fit problem njj
       nh1=nhj1
       nhs1_for_nhj1 = nhs1
       DO nde=1,ndata_on_elem(ne)
          nhs1 = nhs1_for_nhj1
          ns1=0
          DO nn1=1,num_elem_nodes
             DO nk1=1,num_deriv
                nhs1=nhs1+1
                ns1=ns1+1
                PD(ns1)=PSI1(1,nk1,nn1,XIDL(1:2,nde))
             ENDDO !nk1
          ENDDO !nn1
          nhs1 = nhs1_for_nhj1
          DO ns1=1,num_deriv_elem
             nhs1=nhs1+1
             nhs2=0
             DO nhj2=1,num_fit !columns
                nh2=nhj2
                DO ns2=1,num_deriv_elem
                   nhs2=nhs2+1
                   IF(nhj2.EQ.nhj1) THEN !to avoid coupling for now
                      ES(nhs1,nhs2)=ES(nhs1,nhs2)+PD(ns1)*PD(ns2) &
                           *WDL(nh1,nde)*scale_factors_2d(ns1,ne)*scale_factors_2d(ns2,ne)
                   ENDIF !nhj2=nhj1
                ENDDO !ns2
             ENDDO !nhj2
          ENDDO !ns1
       ENDDO !nde

       ns1=0
       nhs1 = nhs1_for_nhj1
       DO nn1=1,num_elem_nodes
          DO nk1=1,num_deriv
             nhs1=nhs1+1
             ns1=ns1+1
             nhs2=0
             DO nhj2=1,num_fit !columns
                ns2=0
                DO nn2=1,num_elem_nodes
                   DO nk2=1,num_deriv
                      nhs2=nhs2+1
                      ns2=ns2+1
                      IF(nhj2.EQ.nhj1) THEN !to avoid coupling for now
                         SUM2=0.0_dp
                         DO ng=1,num_gauss
                            SUM3=0.0_dp
                            DO nu=2,6
                               SUM3=SUM3+ &
                                    PG(ns1,nu,ng)*PG(ns2,nu,ng)*sobelov_wts(nu,ne)
                            ENDDO !nu
                            SUM2=SUM2+SUM3*WG(ng)
                         ENDDO !ng
                         ES(nhs1,nhs2)=ES(nhs1,nhs2)+(SUM2*sobelov_wts(0,ne))* &
                              scale_factors_2d(ns1,ne)*scale_factors_2d(ns2,ne)
                      ENDIF !nhj2=nhj1
                   ENDDO !nk2
                ENDDO !nn2
             ENDDO !nhj2
          ENDDO !nk1
       ENDDO !nn1
    ENDDO !nhj1

  END SUBROUTINE zdes

!!! ##########################################################################

  SUBROUTINE zpze_fit(ne,fit_soln_local,fit_soln)

    USE arrays,ONLY: elem_nodes_2d,elem_versn_2d,&
         scale_factors_2d
!!! dummy arguments
    INTEGER,INTENT(in) :: ne
    REAL(dp) :: fit_soln_local(:,:)
    REAL(dp) :: fit_soln(:,:,:,:)
!!! local variables
    INTEGER :: nh,nk,nn,np,ns,nv

    DO nh=1,num_fit
       ns=0
       DO nn=1,num_elem_nodes
          np=elem_nodes_2d(nn,ne)
          nv=elem_versn_2d(nn,ne)
          DO nk=1,num_deriv
             ns=ns+1
             fit_soln_local(ns,nh)=fit_soln(nk,nv,nh,np)*scale_factors_2d(ns,ne)
          ENDDO !nk
       ENDDO !nn
    ENDDO !nhx

  END SUBROUTINE zpze_fit

!!! ##########################################################################

  SUBROUTINE calculate_ny_maps(npny,num_depvar,nynp,nynr)

    USE arrays,ONLY: node_versn_2d,num_nodes_2d
!!! dummy arguments
    INTEGER :: npny(0:,:),num_depvar,nynp(:,:,:,:),nynr(0:)
!!! local variables
    INTEGER nh,nk,np,nv,ny

    !***  Initialise mapping arrays
    nynp=0
    npny=0
    nynr=0

    !***  Set up mapping arrays
    ny = 0
    nynr(0)=0
    DO nh=1,num_fit
       DO np=1,num_nodes_2d
          DO nv=1,node_versn_2d(np)
             DO nk=1,num_deriv
                ny=ny+1
                nynr(0)=nynr(0)+1
                nynr(nynr(0))=ny
                nynp(nk,nv,nh,np) = ny
                npny(0,ny)=1 !mesh dof is node based
                npny(1,ny)=nk
                npny(2,ny)=nv
                npny(3,ny)=nh
                npny(4,ny)=np
                npny(5,ny)=1
             ENDDO !nk
          ENDDO !nv
       ENDDO !np
    ENDDO !njj
    num_depvar = ny

  END SUBROUTINE calculate_ny_maps

!!! ##########################################################################

  SUBROUTINE define_2d_elements(ELEMFILE)

    USE arrays,ONLY: elems_2d,elem_nodes_2d,elem_versn_2d,node_versn_2d,num_elems_2d
    USE geometry,ONLY: get_final_integer,get_four_nodes,element_connectivity_2d
    CHARACTER(len=*) :: ELEMFILE

    !     Local Variables
    INTEGER :: ierror,ne,nn,noelem,np,number_of_elements
    CHARACTER(len=132) :: ctemp1


    OPEN(10, file=ELEMFILE, status='old')

    read_number_of_elements : DO
       READ(unit=10, fmt="(a)", iostat=ierror) ctemp1
       IF(INDEX(ctemp1, "elements")> 0) THEN
          CALL get_final_integer(ctemp1,number_of_elements)
          EXIT read_number_of_elements
       ENDIF
    END DO read_number_of_elements

    num_elems_2d=number_of_elements
    IF(.NOT.ALLOCATED(elems_2d)) ALLOCATE(elems_2d(num_elems_2d))
    IF(.NOT.ALLOCATED(elem_nodes_2d)) ALLOCATE(elem_nodes_2d(4,num_elems_2d))
    IF(.NOT.ALLOCATED(elem_versn_2d)) ALLOCATE(elem_versn_2d(4,num_elems_2d))

    noelem=1

    read_an_element : DO
       !.......read element number
       READ(unit=10, fmt="(a)", iostat=ierror) ctemp1
       IF(INDEX(ctemp1, "Element")> 0) THEN
          CALL get_final_integer(ctemp1,ne) !get element number
          elems_2d(noelem)=ne
          noelem=noelem+1

          read_element_nodes : DO
             READ(unit=10, fmt="(a)", iostat=ierror) ctemp1
             IF(INDEX(ctemp1, "global")> 0) THEN !found the correct line
                CALL get_four_nodes(ne,ctemp1) !number of versions for node np
                ! note that only the ne'th data of elem_nodes_2d is passed to 'get_four_nodes'
                DO nn=1,4
                   np=elem_nodes_2d(nn,ne)
                   IF(node_versn_2d(np).GT.1)THEN
                      READ(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      READ(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      READ(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      CALL get_final_integer(ctemp1,elem_versn_2d(nn,ne)) !get version#
                   ELSE
                      elem_versn_2d(nn,ne)= 1
                   ENDIF !nversions
                ENDDO !nn
                EXIT read_element_nodes
             ENDIF !index
          END DO read_element_nodes

          IF(noelem.GT.number_of_elements) EXIT read_an_element
       ENDIF

    END DO read_an_element

    CLOSE(10)

    CALL element_connectivity_2d
    CALL line_segments_for_2d_mesh

  END SUBROUTINE define_2d_elements

!!! ##########################################################################

  SUBROUTINE define_xi_closest(data_elem,data_on_elem,ndata_on_elem,data_xi,first)
!!! find the closest xi location on a 2d mesh surface to each data point
    IMPLICIT NONE

!!! dummy arguments
    INTEGER :: data_elem(:),data_on_elem(:,:),ndata_on_elem(:)
    REAL(dp) :: data_xi(:,:)
    LOGICAL,INTENT(in) :: first
!!! local variables
    INTEGER :: i,n_check,ne_checklist(5),IT,ITMAX=20,nd,ne,neadj,nelast,neold,ni,nj
    REAL(dp) :: sqmax,sqnd,temp,xe(num_deriv_elem,num_coords),xi(3)
    REAL(dp),ALLOCATABLE :: sq(:)
    LOGICAL :: found

    INTEGER :: n_data
    CHARACTER(len=200) :: exfile
    CHARACTER(len=1) :: string_ne1
    CHARACTER(len=2) :: string_ne2

    ALLOCATE(sq(num_data))

    sqmax = 1.0e4_dp*1.0e4_dp

    !  initialise
    sq = 0.0_dp
    xi = 0.5_dp

!!! start by finding the closest centre of an element to each data point
    !    do ne = 1,num_elems_2d
    !       call xpxe(ne,xe)
    !       do nd = 1,num_data
    !          sqnd = 0.0_dp
    !          do nj=1,num_coords
    !             temp = pxi(1,xi,xe(1,nj))-data_xyz(nj,nd)
    !             sqnd = sqnd+temp**2
    !          enddo !nj
    !          if(data_elem(nd).eq.0.or.sqnd.lt.sq(nd)) then
    !             data_xi(1:2,nd) = xi(1:2)
    !             if(nd.eq.3976)then
    !                write(*,*) 'initial closest=',ne
    !             endif
    !             data_elem(nd) = ne
    !             sq(nd) = sqnd
    !          endif
    !       enddo ! nd
    !    enddo ! ne

    IF(first)THEN ! check every element for every data point
       DO nd = 1,num_data
          sqmax = 1.0e4_dp*1.0e4_dp
          DO ne = 1,num_elems_2d
             xi = 0.5_dp
             CALL xpxe(ne,xe)
             found = .FALSE.
             CALL project_orthogonal(nd,SQND,xe,xi,found)
             IF(ABS(xi(1)).GE.-zero_tol.AND.ABS(xi(1)).LT.1.0_dp+zero_tol.AND. &
                  ABS(xi(2)).GE.zero_tol.AND.ABS(xi(2)).LT.1.0_dp+zero_tol) THEN
                IF(sqnd.LT.sqmax)THEN
                   sqmax = sqnd
                   data_xi(1:2,nd) = xi(1:2)
                   data_elem(nd) = ne
                   SQ(nd) = SQND
                ENDIF
             ENDIF !FOUND
          ENDDO
       ENDDO ! nd
    ELSE

       DO nd = 1,num_data
          ne = data_elem(nd)
          ne_checklist(1) = ne
          n_check = 1
          IF(elem_cnct_2d(-1,0,ne).NE.0)THEN
             n_check = n_check + 1
             ne_checklist(n_check) = elem_cnct_2d(-1,1,ne)
          ENDIF
          IF(elem_cnct_2d(1,0,ne).NE.0)THEN
             n_check = n_check + 1
             ne_checklist(n_check) = elem_cnct_2d(1,1,ne)
          ENDIF
          IF(elem_cnct_2d(-2,0,ne).NE.0)THEN
             n_check = n_check + 1
             ne_checklist(n_check) = elem_cnct_2d(-2,1,ne)
          ENDIF
          IF(elem_cnct_2d(2,0,ne).NE.0)THEN
             n_check = n_check + 1
             ne_checklist(n_check) = elem_cnct_2d(2,1,ne)
          ENDIF
          sqmax = 1.0e4_dp*1.0e4_dp
          DO i = 1,n_check
             ne = ne_checklist(i)
             IF(i.EQ.1)THEN
                xi(1:2) = data_xi(1:2,nd)
             ELSE
                xi = 0.5_dp
             ENDIF
             CALL xpxe(ne,xe)
             found = .TRUE. !find nearest point in element
             CALL project_orthogonal(nd,sqnd,xe,xi,found)
             IF(ABS(xi(1)).GE.-zero_tol.AND.ABS(xi(1)).LT.1.0_dp+zero_tol.AND. &
                  ABS(xi(2)).GE.zero_tol.AND.ABS(xi(2)).LT.1.0_dp+zero_tol) THEN
                IF(sqnd.LT.sqmax)THEN
                   sqmax = sqnd
                   data_xi(1:2,nd) = xi(1:2)
                   data_elem(nd)=ne
                   sq(nd) = sqnd
                ENDIF
             ENDIF
          ENDDO !i
       ENDDO ! nd
    ENDIF

    ndata_on_elem=0
    DO nd = 1,num_data
       ne = data_elem(nd)
       IF(ne > 0) THEN
          ndata_on_elem(ne) = ndata_on_elem(ne)+1
          IF(ndata_on_elem(ne) > nmax_data_elem)THEN
             WRITE(*,'(''Number of data points on element'',i6,'' exceeds'',i6)') ne,nmax_data_elem
             WRITE(*,'(''Reduce the data point density: no advantage in using this many'')')
             STOP
          ENDIF
          data_on_elem(ne,ndata_on_elem(ne)) = nd
       ENDIF
    ENDDO

    DEALLOCATE(sq)

    exfile = 'temp.exdata'
    OPEN(10, file = exfile, status = 'replace')

    DO ne = 1,num_elems_2d
       !**   write the group name
       IF(ne.LT.10)THEN
          WRITE(string_ne1,'(i1)') ne
          WRITE(10,'( '' Group name: '',A)') 'datapoints_'//string_ne1
       ELSE
          WRITE(string_ne2,'(i2)') ne
          WRITE(10,'( '' Group name: '',A)') 'datapoints_'//string_ne2
       ENDIF
       WRITE(10,'(1X,''#Fields=1'')')
       WRITE(10,'(1X,''1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
       WRITE(10,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
       WRITE(10,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
       WRITE(10,'(1X,''  z.  Value index= 3, #Derivatives=0'')')

       DO n_data = 1,ndata_on_elem(ne)
          nd = data_on_elem(ne,n_data)
          WRITE(10,'(1X,''Node: '',I9)') nd
          WRITE(10,'(1X,3E13.5)')  (data_xyz(nj,nd),nj=1,num_coords)
       ENDDO !n_data
    ENDDO !ne
    CLOSE(10)

  END SUBROUTINE define_xi_closest


!!! ##########################################################################

  SUBROUTINE solve_geometry_fit(data_on_elem,ndata_on_elem,num_depvar,&
       elem_list,not_1,not_2,npny,nynp,nynr,nyny,data_xi,cyny,sobelov_wts,&
       fit_soln,fix_bcs)

    USE arrays,ONLY: node_xyz_2d
!!! dummy arguments
    INTEGER :: data_on_elem(:,:),ndata_on_elem(:),not_1,not_2,num_depvar,&
         elem_list(0:),npny(0:,:),nynp(:,:,:,:),nynr(0:),nyny(0:,:)
    REAL(dp) :: data_xi(:,:),cyny(0:,:),sobelov_wts(0:,:),fit_soln(:,:,:,:)
    LOGICAL :: fix_bcs(:)
!!! local variables
    INTEGER :: l,LGE2(3*16,2),ne,nh,nhs1,nhs2,NHST(2), &
         nk,no1,no2,no_nynr1,no_nynr2,noy1,noy2,np,nv,ny1,ny2,ny3,nyo1,nz,nzz
    INTEGER,ALLOCATABLE :: nony(:,:,:)
    INTEGER,ALLOCATABLE :: nyno(:,:,:)
    REAL(dp) :: co1,co2,ER(num_fit*num_deriv_elem),ES(3*16,3*16),&
         fit_soln_local(16,3),PG(16,6,9),WG(9)
    REAL(dp),ALLOCATABLE :: cony(:,:,:)
    REAL(dp),ALLOCATABLE :: cyno(:,:,:)
    REAL(dp),ALLOCATABLE :: GR(:)      ! right-hand-side vector
    REAL(dp),ALLOCATABLE :: GRR(:)     ! reduced right-hand-side vector
    REAL(dp),ALLOCATABLE :: incr_soln(:)         ! current solution returned from solver
    LOGICAL :: FIRST_A,UPDATE_MATRIX

    ! make all of these allocatable!
    REAL(dp),DIMENSION(nsize_gkk) :: GKK
    REAL(dp),DIMENSION(nsize_gkk) :: GK
    ! doesn't like allocating these! gives different answer for errors
    !    real(dp),allocatable :: GKK(:)
    !    real(dp),allocatable :: GK(:)
    REAL(dp),DIMENSION(3,nmax_data_elem) :: WDL
    REAL(dp),DIMENSION(2,nmax_data_elem) :: XIDL


    WG = [7.7160493827160628e-2_dp, 0.12345679012345677_dp, 7.7160493827160628e-2_dp,&
         0.12345679012345677_dp, 0.19753086419753044_dp, 0.12345679012345677_dp,&
         7.7160493827160628e-2_dp, 0.12345679012345677_dp, 7.7160493827160628e-2_dp]

    ALLOCATE(incr_soln(num_depvar))
    ALLOCATE(nony(0:1,num_depvar,2))
    ALLOCATE(nyno(0:5,num_depvar,2))
    ALLOCATE(cony(0:1,num_depvar,2))
    ALLOCATE(cyno(0:5,num_depvar,2))
    ALLOCATE(GR(num_depvar))
    ALLOCATE(GRR(num_depvar))
    !    allocate(GK(num_depvar*num_depvar))
    !    allocate(GKK(num_depvar*num_depvar))

    CALL gauss1(PG)

    UPDATE_MATRIX=.TRUE.
    FIRST_A=.TRUE.

    GR=0.0_dp

    DO l=1,elem_list(0) !loop over elements in the fit
       ne=elem_list(l)
       CALL melgef(LGE2,ne,NHST,nynp)
       ER=0.0_dp
       ES=0.0_dp
       CALL zpze_fit(ne,fit_soln_local,fit_soln) !gets fit_soln_local for element ne
       CALL zder(data_on_elem,ndata_on_elem,ne,data_xi,ER,PG,WDL,WG,&
            sobelov_wts,XIDL,fit_soln_local)
       CALL zdes(ndata_on_elem,ne,ES,PG,WDL,WG,sobelov_wts,XIDL)

       !*** Assemble element stiffness matrix into global system.
       DO nhs1=1,NHST(1) !3 dependent variables
          ny1=IABS(LGE2(nhs1,1))
          IF(ny1.EQ.0)THEN
             WRITE(*,'('' No dependent variable for node in element'',i6,'': are &
                  &you sure you have set up versions correctly?'')') ne
             STOP
          ENDIF
          GR(ny1)=GR(ny1)+ER(nhs1)
          DO nhs2=1,NHST(2) !3 dependent variables
             ny2=IABS(LGE2(nhs2,2))
             nz=ny1+(ny2-1)*num_depvar
             GK(nz)=GK(nz)+ES(nhs1,nhs2)
          ENDDO !nhs2
       ENDDO !nhs1
    ENDDO !l (ne)

    !*** Calculate solution mapping arrays for the current fit variable
    CALL globalf(nony,not_1,not_2,npny,nyno,nynp,nyny,cony,cyno,cyny,fix_bcs)

    IF(NOT_2.EQ.0) THEN
       WRITE(*,'('' >>The number of unknowns is zero'')')
       STOP
    ENDIF

    !----------------------- generate reduced system -----------------------

    GKK=0.0_dp
    GRR=0.0_dp

    !*** generate the reduced system of equations
    DO no_nynr1=1,nynr(0) !loop global rows of GK
       ny1=nynr(no_nynr1) !is row #
       DO noy1=1,nony(0,ny1,1) !loop over #no's attached to ny1
          no1=nony(noy1,ny1,1) !no# attached to row ny1
          co1=cony(noy1,ny1,1) !coupling coeff for row mapping
          !                     ie row_no1=a*row_ny1+b*row_ny2
          GRR(no1)=GRR(no1)+GR(ny1)*co1 !get reduced R.H.S.vector
          DO no_nynr2=1,nynr(0) !loop over #cols of GK
             ny2=nynr(no_nynr2) !is global variable #
             ny3=getnyr(npny,ny2,nynp)
             !local GK var #
             nz=ny1+(ny3-1)*num_depvar
             IF(nz.NE.0) THEN
                DO noy2=1,nony(0,ny2,2) !loop over #no's for ny2
                   no2=nony(noy2,ny2,2) !no# attached to ny2
                   co2=cony(noy2,ny2,2) !coup coeff col mapping
                   !                     i.e. var_no1=a*var_ny1+b*var_ny2
                   nzz=no1+(no2-1)*NOT_1
                   WRITE(*,*) 'nzz',nzz
                   IF(nzz.NE.0) GKK(nzz)=GKK(nzz)+GK(nz)*co1*co2
                ENDDO !noy2
             ENDIF
          ENDDO !no_nynr2
       ENDDO !noy1
    ENDDO !no_nynr1

    !-------------- solve reduced system of linear equations ---------------
    !Commented out since subroutines called further are temporarily unavailable
    WRITE(*,*) NOT_1,NOT_2, num_depvar, SIZE(GKK), GKK(1:10), SIZE(GRR), GRR(1:10)
    WRITE(*,*) SIZE(incr_soln), incr_soln(1:10)
    !pause
    !call direct_solver(NOT_1,NOT_1,NOT_2,num_depvar,GKK,GRR,incr_soln,FIRST_A)

    DO no1=1,NOT_2 ! for each unknown
       DO nyo1=1,nyno(0,no1,2)
          ny1=nyno(nyo1,no1,2) ! the dependent variable number
          co1=cyno(nyo1,no1,2) ! the weighting for mapped variables
          nk=npny(1,ny1)     ! derivative number
          nv=npny(2,ny1)     ! version number
          nh=npny(3,ny1)     ! dependent variable number
          np=npny(4,ny1)     ! node number
          fit_soln(nk,nv,nh,np) = fit_soln(nk,nv,nh,np) + incr_soln(no1)*co1
          ! current fit solution = previous + increment
       ENDDO !nyo1
    ENDDO !no1

    !*** Copy the fitted solution back into node_xyz_2d
    node_xyz_2d = fit_soln

    DEALLOCATE(incr_soln)
    DEALLOCATE(nony)
    DEALLOCATE(nyno)
    DEALLOCATE(cony)
    DEALLOCATE(cyno)
    DEALLOCATE(GR)
    DEALLOCATE(GRR)
    !    deallocate(GK)
    !    deallocate(GKK)

  END SUBROUTINE solve_geometry_fit

!!! ##########################################################################

  SUBROUTINE project_orthogonal(nd,SQ,xe,xi,inelem)

    USE arrays,ONLY: data_xyz
!!! dummy arguments
    INTEGER :: nd
    REAL(dp) :: sq,xe(num_deriv_elem,num_coords),xi(:)
    LOGICAL :: inelem
!!! local variables
    INTEGER :: IT,ITMAX,BOUND(2),it2,ni,nifix,nj
    REAL(dp) :: LOOSE_TOL=1.0e-6_dp
    REAL(dp) :: DELTA,DET,D2SQV2,D2SQVW2,D2SQXI(2,2),D2ZXI(3,2,2),DSQXI(2), &
         DSQXI1,DSQXI2,DSQV,DSQVW,DZ(3),DZXI(3,2),EVMIN,EVMAX,H(2), &
         MU,SQLIN,SQDIFF,SQDPRED,TEMP,TEMP1,TEMP2,TOL, &
         TOL2,V(2),V1,V2,VMAX=1.0_dp,W,XILIN(2),Z(3)
    LOGICAL :: CONVERGED,ENFORCE(2),FREE,NEWTON

    ITMAX = 10 ! max # iterations to use
    DELTA = VMAX/4.0_dp
    TOL = 5.0_dp*LOOSE_TOL !must be > sqrt(eps) or SQLIN<=SQ check may not work
    TOL2 = TOL**2
    SQ = 0.0_dp
    DO nj=1,num_coords
       Z(nj) = PXI(1,XI,XE(1,nj))
       DZ(nj) = Z(nj)-data_xyz(nj,nd)
       SQ = SQ+DZ(nj)**2
    ENDDO !nj
    IT=0
    CONVERGED=.FALSE.
    DO WHILE(.NOT.CONVERGED.AND.IT.LT.ITMAX)
       DSQXI = 0.0_dp
       DO nj=1,num_coords
          DZXI(nj,1) = PXI(2,XI,XE(1,nj))
          DZXI(nj,2) = PXI(4,XI,XE(1,nj))
          DSQXI(1) = DSQXI(1)+DZXI(nj,1)*DZ(nj)
          DSQXI(2) = DSQXI(2)+DZXI(nj,2)*DZ(nj)
       ENDDO !nj
       DO ni=1,2
          IF(dabs(xi(ni)) < zero_tol) THEN
             BOUND(ni)=1
             ENFORCE(ni)=DSQXI(ni).GE.0.0_dp
          ELSE IF(dabs(xi(ni)-1.0_dp) < zero_tol) THEN
             BOUND(ni)=-1
             ENFORCE(ni)=DSQXI(ni).LE.0.0_dp
          ELSE
             BOUND(ni)=0
             ENFORCE(ni)=.FALSE.
          ENDIF
       ENDDO !ni
       IF(ENFORCE(1).AND.ENFORCE(2)) EXIT

       D2SQXI = 0.0_dp
       DO nj=1,num_coords
          D2ZXI(nj,1,1) = PXI(3,XI,XE(1,nj))
          D2ZXI(nj,1,2) = PXI(5,XI,XE(1,nj))
          D2ZXI(nj,2,2) = PXI(5,XI,XE(1,nj))
          D2SQXI(1,1) = D2SQXI(1,1)+DZXI(nj,1)*DZXI(nj,1)+D2ZXI(nj,1,1)*DZ(nj)
          D2SQXI(1,2) = D2SQXI(1,2)+DZXI(nj,1)*DZXI(nj,2)+D2ZXI(nj,1,2)*DZ(nj)
          D2SQXI(2,2) = D2SQXI(2,2)+DZXI(nj,2)*DZXI(nj,2)+D2ZXI(nj,2,2)*DZ(nj)
       ENDDO !nj

       !       A Newton step is taken if the condition of the Hessian
       !       guarantees that the step will be within the trust region.
       !       Otherwise the Hessian is shifted towrds a diagonal matrix to
       !       shift the step towards steepest descent.  Usually it is a much
       !       better direction than steepest descent.  I think it is close to
       !       the best direction in the trust region.

       !***    Find the smallest eigen value of the Hessian.
       DSQXI2 = DSQXI(1)**2+DSQXI(2)**2
       DSQXI1 = DSQRT(DSQXI2)
       TEMP1 = (D2SQXI(1,1)+D2SQXI(2,2))/2.0_dp
       TEMP2 = DSQRT(((D2SQXI(1,1)-D2SQXI(2,2))/2.0_dp)**2+D2SQXI(1,2)**2)
       EVMIN = TEMP1-TEMP2
       EVMAX = TEMP1+TEMP2
       IF(DSQXI1.LT.TOL2) EXIT
       DO it2=1,ITMAX
          TEMP = DSQXI1/DELTA
          NEWTON = EVMIN.GE.TEMP
          IF(NEWTON) THEN !Newton is safe
             H(1) = D2SQXI(1,1)
             H(2) = D2SQXI(2,2)
             DET = EVMIN*EVMAX
          ELSE
             !***        Shift eigenvalues to restrict step
             MU = TEMP-EVMIN
             H(1) = D2SQXI(1,1)+MU
             H(2) = D2SQXI(2,2)+MU
             DET = TEMP*(EVMAX+MU)
          ENDIF
          V(1) = -(H(2)*DSQXI(1)-D2SQXI(1,2)*DSQXI(2))/DET
          V(2) = (D2SQXI(1,2)*DSQXI(1)-H(1)*DSQXI(2))/DET
          V2 = V(1)**2+V(2)**2
          DSQV = DSQXI(1)*V(1)+DSQXI(2)*V(2)
          !         This checks that numerical errors have not
          !         prevented the step direction being a descent direction.

          IF(DSQV**2.LT.DSQXI2*V2*TOL2) THEN !try a smaller trust region
             DELTA = DELTA/10.0_dp
          ELSE !step is good
             !***        Check feasible and limit step size
             FREE = .TRUE.
             DO ni=1,2
                IF(BOUND(ni)/=0 )THEN
                   IF(BOUND(ni)>0.EQV.V(ni)<0.0_dp) THEN
                      FREE=.FALSE.
                      nifix=ni
                   ENDIF
                ENDIF
             ENDDO
             W=1.0_dp
             IF(FREE) THEN
                V1=DSQRT(V2) !currently < DELTA
                D2SQV2=V(1)*(V(1)*D2SQXI(1,1)+2.0_dp*V(2)*D2SQXI(1,2)) &
                     +V(2)**2*D2SQXI(2,2)
                IF(.NOT.NEWTON) THEN
                   !               Try to step to estimate of minimum along line
                   IF(V1.GT.0.0_dp) THEN
                      W=DELTA/V1
                      IF(D2SQV2.GT.0.0_dp) THEN !minimum exists
                         W=DMIN1(W,-DSQV/D2SQV2) !minimum if within trust region
                      ENDIF
                   ENDIF
                ENDIF !newton
             ELSE
                IF(ENFORCE(2)) THEN !gradient suggests must use ni=1
                   nifix=2
                ELSE IF(ENFORCE(1)) THEN !gradient suggests must use ni=2
                   nifix=1
                   !else Gradient points into element.
                   !               Fix the direction that prevented the step.
                   !               P.d. Hessian guarantees there is only one of these
                   !               for this type of gradient.
                ENDIF
                ni=3-nifix
                IF(.NOT.INELEM) THEN
                   !***            If stepping predominantly out of element then exit
                   nifix=3-ni
                   IF(DABS(V(nifix)).GT.DABS(DSQXI(ni)/H(ni))) THEN
                      XI(nifix)=XI(nifix)+V(nifix)
                      EXIT
                   ENDIF
                ENDIF
                V(nifix)=0.0_dp
                IF(D2SQXI(ni,ni).GT.0.0_dp) THEN !minimum exists
                   V(ni)=-DSQXI(ni)/D2SQXI(ni,ni)
                   V1=DABS(V(ni))
                   NEWTON=V1.LE.DELTA
                ENDIF
                IF(.NOT.NEWTON) THEN
                   V(ni)=-DSIGN(DELTA,DSQXI(ni))
                   V1=DELTA
                ENDIF
                V2=V1*V1
                DSQV=DSQXI(ni)*V(ni)
                D2SQV2=V2*D2SQXI(ni,ni)
             ENDIF !free
             !***        First half of convergence test.
             !           Should be before boundary colllision check
             CONVERGED=V1*W.LT.TOL
             !***        Try the step.  (Name: XILIN is historical)
             XILIN(1)=XI(1)+V(1)*W
             XILIN(2)=XI(2)+V(2)*W
             !***        Test for boundary collision
             DO ni=1,2
                IF(XILIN(ni).LT.0.0_dp) THEN
                   XILIN(ni)=0.0_dp
                   W=XI(ni)/(-V(ni))
                   XILIN(3-ni)=XI(3-ni)+V(3-ni)*W
                ELSE IF(XILIN(ni).GT.1.0_dp) THEN
                   XILIN(ni)=1.0_dp
                   W=(1.0_dp-XI(ni))/V(ni)
                   XILIN(3-ni)=XI(3-ni)+V(3-ni)*W
                ENDIF
             ENDDO !ni
             !***        Calculate new distance
             SQLIN=0.0_dp
             DO nj=1,num_coords
                Z(nj)=PXI(1,XILIN,XE(1,nj))
                DZ(nj)=Z(nj)-data_xyz(nj,nd)
                SQLIN=SQLIN+DZ(nj)**2
             ENDDO
             !***        Second half of convergence test.
             CONVERGED=CONVERGED.AND.DABS(SQ-SQLIN)/(1.0_dp+SQ).LE.TOL
             IF(CONVERGED) GO TO 5
             DSQVW=DSQV*W !<0
             D2SQVW2=0.5_dp*D2SQV2*W*W !1/2 for computational efficiency
             SQDIFF=0.5_dp*(SQLIN-SQ) !1/2 because derivs are for SQ/2
             SQDPRED=SQDIFF-DSQVW-D2SQVW2
             !***        Exit loop if decrease is satisfactory
             IF(SQDIFF.LE.0.25_dp*DSQVW) THEN
                IF(NEWTON) THEN
                   DELTA=V1 !next step smaller unless this is increased
                ELSE IF(W.GE.1.0_dp) THEN
                   !               If the quadratic model is good increase trust region size
                   IF(SQDPRED.LT.-0.1_dp*SQDIFF) THEN
                      DELTA=DMIN1(VMAX,DELTA*2.0_dp)
                   ENDIF
                ENDIF
                GO TO 5
             ENDIF
             !***        Calculate new trust region size from an estimate of the
             !***        minimum along the step direction using a cubic approximation.
             TEMP=-3.0_dp*SQDPRED !<0
             DELTA=W*V1*(D2SQVW2-DSQRT(D2SQVW2**2+TEMP*DSQVW))/TEMP !>0
!!! note: was getting delta=0 because d2sqvw2=0. the following avoids, but is it correct?
             IF(ABS(delta).LE.zero_tol)THEN
                delta = 0.01_dp
             ENDIF
          ENDIF !DSQV**2.LT.DSQXI2*V2*TOL2
       ENDDO !it2

5      SQ=SQLIN
       XI(1)=XILIN(1)
       XI(2)=XILIN(2)
       IT=IT+1
    ENDDO

    IF(.NOT.inelem.AND.XI(1)>=0.0_dp.AND.XI(1)<=1.0_dp.AND.XI(2)>=0.0_dp.AND.XI(2)<=1.0_dp) THEN
       inelem=.TRUE.
    ENDIF

  END SUBROUTINE project_orthogonal

!!! ##########################################################################



END MODULE surface_fitting
