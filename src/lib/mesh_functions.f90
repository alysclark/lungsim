MODULE mesh_functions

!!! Subroutines and functions for general calculations. Not specific to any
!!! particular application, although it is expected these will generally be used
!!! for calculations to do with a mesh (not fields or solutions).
!!! Any function that is used by more than one module should appear in here.
!!! ALL subroutines and functions in this module are public.

  USE other_consts   !! pi

  IMPLICIT NONE

  PRIVATE

  PUBLIC  area_between_three_points,area_between_two_vectors,calc_branch_direction,&
       angle_btwn_vectors,calc_scale_factors_2d,check_colinear_points,cross_product,&
       distance_between_points,make_plane_from_3points,mesh_a_x_eq_b,ph3,pl1,&
       point_internal_to_surface,scalar_product_3,scalar_triple_product,scale_mesh,&
       unit_norm_to_plane_two_vectors,unit_norm_to_three_points,unit_vector,&
       vector_length,volume_internal_to_surface

CONTAINS

!!! list of subroutines

  ! calc_branch_direction
  ! ..... calculates direction of branch and stores in elem_direction

  ! calc_scale_factors_2d
  ! ..... calculate the scale factors for a 2d mesh

  ! make_plane_from_3points
  ! ..... finds the equation of a plane in 3D and a vector normal to the plane from three
  ! ..... non-colinear points.

  ! scale_mesh
  ! ..... multiply mesh (coordinates and derivatives) by a constant

!!! list of functions

  ! angle_btwn_vectors
  ! .... returns the angle between two vectors

  ! check_colinear_points
  ! .... returns true or false for whether 3 points are colinear

  ! cross_product ***
  ! .... returns the vector cross product of A*B

  ! distance_between_points ***
  ! .... calculates the distance between two arbitrary points

  ! mesh_a_x_eq_b ***
  ! .... solves a small matrix system

  ! scalar_product_3 ***
  ! .... dot product of two 3x1 vectors

  ! unit_vector ***
  ! .... Calculates the unit vector for an arbitrary 3x1 vector

  ! vector_length ***
  ! .... Calculates the length of a 3x1 vector

!!!#####################################################################

  SUBROUTINE calc_branch_direction(ne)

!!! calculates the direction of element ne and stores in elem_direction

    USE arrays!,only: elem_direction,elem_nodes,node_xyz

    INTEGER,INTENT(in) :: ne

    INTEGER :: np_end,np_start
    REAL(dp) :: length

    np_start = elem_nodes(1,ne)
    np_end = elem_nodes(2,ne)

    length = distance_between_points(node_xyz(1,np_end),node_xyz(1,np_start))
    elem_direction(1:3,ne)=(node_xyz(1:3,np_end)-node_xyz(1:3,np_start))/length

  END SUBROUTINE calc_branch_direction

!!! ##########################################################################

  SUBROUTINE calc_scale_factors_2d(sf_option)

!!! calculates the arclengths and scale factors for 2d surface elements,
!!! stores in scale_factors_2d

    USE arrays !only: arclength,elem_lines_2d,elem_nodes_2d,lines_2d,lines_in_elem,&
    !line_versn_2d,nodes_in_line,node_xyz_2d,num_elems_2d,num_lines_2d,scale_factors_2d,dp

    CHARACTER(len=4),INTENT(in) :: sf_option
!!! local variables
    INTEGER,PARAMETER :: num_deriv = 4
    INTEGER :: ido(num_deriv,2),IG(4),it,IT_count,ITMAX=20,k,N,NAE,ne,&
         ng,NGA=4,NI1(3),ni,ni2,nj,nk,nk2,nl,nn,nn2,NNK,no_nl,&
         np,ns,nv,NNL(2,4)
    REAL(dp) :: DA,SUM1,SUM2,SUM3,SUM4,W,WG_LOCAL(10),XA_LOCAL(4,3),XI,&
         XIGG(10),XN_LOCAL(2,3,4)
    LOGICAL :: FOUND

    XIGG = [0.6_dp,0.2113248654051_dp,0.7886751345948_dp,0.1127016653792_dp,&
         0.6_dp,0.8872983346207_dp,0.0694318442029_dp,0.3300094782075_dp,&
         0.6699905217924_dp,0.9305681557970_dp]
    WG_LOCAL = [1.0_dp,0.6_dp,0.6_dp,0.2777777777778_dp,0.4444444444444_dp,&
         0.2777777777778_dp,0.1739274225687_dp,0.3260725774313_dp,&
         0.3260725774313_dp,0.1739274225687_dp]
    IG = [0,1,3,6]
    ido = RESHAPE ([1,2,1,2,1,1,2,2],SHAPE(ido))
    NI1 = [1,2,1]
    NNL = RESHAPE([1,2,3,4,1,3,2,4],SHAPE(NNL))

    IF(.NOT.ALLOCATED(scale_factors_2d)) ALLOCATE(scale_factors_2d(16,num_elems_2d))

    SELECT CASE (sf_option)
    CASE ('unit')
       scale_factors_2d = 1.0_dp
    CASE('arcl')
       DO no_nl=1,num_lines_2d !loop over global lines
          nl=lines_2d(no_nl)
          ni = nodes_in_line(1,0,nl)
          DO n = 1,2                  !for each node on the line
             np = nodes_in_line(n+1,1,nl)       !np1,np2
             DO nj = 1,3
                nv = line_versn_2d(N,nj,nl)
                XN_LOCAL(1,nj,n) = node_xyz_2d(1,nv,nj,np)
                ne = lines_in_elem(1,nl)
                ni2 = 1+MOD(ni,2)
                nn = 1
                FOUND =.FALSE.
                DO WHILE((nn.LE.4).AND.(.NOT.FOUND))
                   IF(np.EQ.elem_nodes_2d(nn,ne))THEN
                      FOUND=.TRUE.
                   ELSE
                      nn=nn+1
                   ENDIF
                ENDDO
                DO nk=2,4           !dxi1, dxi2, d2xi1xi2
                   IF(IDO(nk,ni).EQ.2.AND.IDO(nk,ni2).EQ.1) THEN
                      XN_LOCAL(2,nj,n) = node_xyz_2d(nk,nv,nj,np)
                   ENDIF
                ENDDO !nk
             ENDDO !nj
          ENDDO !n

          SUM2=0.0_dp
          DO ng=1,NGA
             XI=XIGG(IG(NGA)+ng)
             W=WG_LOCAL(IG(NGA)+ng)
             DO nj=1,3
                DO k=1,2
                   XA_LOCAL(k,nj)=PL1(1,k,XI)*XN_LOCAL(1,nj,1) &
                        +PL1(2,k,XI)*XN_LOCAL(1,nj,2)
                ENDDO
             ENDDO

             SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
             SUM2=SUM2+W*DSQRT(SUM1)
          ENDDO !ng

          arclength(1:3,nl)=SUM2

          it=0
          iterative_loop : DO
             it=it+1
             IT_count=it
             SUM3=0.0_dp
             SUM4=0.0_dp
             DO ng=1,NGA
                XI=XIGG(IG(NGA)+ng)
                W=WG_LOCAL(IG(NGA)+ng)
                DO nj=1,3
                   DO k=1,2
                      XA_LOCAL(k,nj)=0.0_dp
                      DO n=1,2
                         XA_LOCAL(k,nj)=XA_LOCAL(k,nj)+ &
                              PH3(n,1,k,XI)*XN_LOCAL(1,nj,n) &
                              +PH3(n,2,k,XI)*XN_LOCAL(2,nj,n)*arclength(n,nl)
                      ENDDO
                   ENDDO
                   XA_LOCAL(3,nj)=0.0_dp
                   DO n=1,2
                      XA_LOCAL(3,nj)=XA_LOCAL(3,nj)+ &
                           PH3(n,2,2,XI)*XN_LOCAL(2,nj,n)
                   ENDDO
                   XA_LOCAL(4,nj)=0.0_dp
                   DO n=1,2
                      XA_LOCAL(4,nj)=XA_LOCAL(4,nj)+ &
                           PH3(n,2,1,XI)*XN_LOCAL(2,nj,n)
                   ENDDO
                ENDDO
                SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
                SUM2=0.0_dp
                DO nj=1,3
                   SUM2=SUM2+XA_LOCAL(2,nj)*XA_LOCAL(3,nj)
                ENDDO                  !nj
                SUM3=SUM3+W*DSQRT(SUM1)
                IF(SUM1.GT.1.0e-6_dp) SUM4=SUM4+W*SUM2/DSQRT(SUM1)
             ENDDO                     !ng
             DA=-(arclength(3,nl)-SUM3)/(1.0_dp-SUM4)
             IF(DABS(DA).GT.1.0e+6_dp) THEN
                arclength(3,nl)=1.0_dp
                EXIT iterative_loop
             ENDIF

             arclength(3,nl) = arclength(3,nl)+DA      !is new arclength
             arclength(1:2,nl) = arclength(3,nl)

             IF(it.EQ.ITMAX) EXIT iterative_loop

          ENDDO iterative_loop      !iteration
       ENDDO !loop over lines

       scale_factors_2d = 1.0_dp !initialise

       DO ne=1,num_elems_2d
          DO NAE=1,4
             nl = elem_lines_2d(NAE,ne)
             IF(nl /= 0)THEN
                ni = nodes_in_line(1,0,nl)
                ni2 = NI1(ni+1)
                DO N=1,2
                   nn=NNL(N,NAE)
                   ns=0
                   DO nn2=1,nn-1
                      DO nk2=1,num_deriv
                         ns=ns+1
                      ENDDO
                   ENDDO
                   DO nk=2,num_deriv
                      IF(IDO(nk,ni2).EQ.1) THEN
                         scale_factors_2d(nk+ns,ne) = arclength(N,nl)
                         IF(DABS(scale_factors_2d(nk+ns,ne)).LT.1.0e-6_dp) scale_factors_2d(nk+ns,ne) = 1.0_dp
                      ENDIF
                   ENDDO !nk
                ENDDO !N=1,2
             ENDIF
          ENDDO !NAE (nl)

          NNK=0
          DO ns=1,4
             scale_factors_2d(NNK+4,ne)=scale_factors_2d(NNK+2,ne)*scale_factors_2d(NNK+3,ne)
             NNK=NNK+4
          ENDDO !nn

       ENDDO !noelem (ne)
    END SELECT

  END SUBROUTINE calc_scale_factors_2d

!!!###############################################################

  SUBROUTINE make_plane_from_3points(NORML,NORMALTYPE,POINT1,POINT2,POINT3)

    !###    make_plane_from_3points finds the equation of a plane in three
    !###    dimensions and a vector normal to the plane from three
    !###    non-colinear points.
    !###    NORMALTYPE=1 for raw normal and plane equation
    !###    NORMALTYPE=2 for unit normal and plane equation
    !###    The coefficients represent aX + bY + cZ + d = 0
    !###    NORML(1)=a,NORML(2)=b,NORML(3)=c,NORML(4)=d
    USE arrays

    !     Parameter list
    INTEGER :: NORMALTYPE
    REAL(dp) :: POINT1(3),POINT2(3),POINT3(3),NORML(4)
    !     Local variables
    REAL(dp) :: DifF1(3),DifF2(3),NORMSIZE
    LOGICAL :: COLINEAR


    ! Check for colinearity
    COLINEAR=.FALSE.
    colinear = check_colinear_points(POINT1,POINT2,POINT3)
    IF(.NOT.COLINEAR) THEN
       DifF1(1:3)=POINT2(1:3)-POINT1(1:3)
       DifF2(1:3)=POINT2(1:3)-POINT3(1:3)

       NORML(1)=(DifF1(2)*DifF2(3))-(DifF1(3)*DifF2(2))
       NORML(2)=(DifF1(3)*DifF2(1))-(DifF1(1)*DifF2(3))
       NORML(3)=(DifF1(1)*DifF2(2))-(DifF1(2)*DifF2(1))

       IF(NORMALTYPE.EQ.2) THEN
          NORMSIZE = vector_length(NORML)
          NORML(1:3)=NORML(1:3)/NORMSIZE
       ENDIF

       NORML(4) = -scalar_product_3(NORML,POINT1)

    ELSE !Colinear

       WRITE(*,*) ' COLINEAR points in make_plane_from_3points '
       NORML = 0.0_dp
    ENDIF
  END SUBROUTINE make_plane_from_3points

!!!##################################################

  SUBROUTINE scale_mesh(scaling,TYPE)

    USE arrays!,only: node_xyz,node_xyz_2d,scale_factors_2d

    REAL(dp),INTENT(in) :: scaling
    CHARACTER(len=2),INTENT(in) :: TYPE

    SELECT CASE(TYPE)
    CASE('1d')
       node_xyz = node_xyz * scaling
    CASE('2d')
       node_xyz_2d = node_xyz_2d * scaling
       node_xyz_2d(4,:,:,:) = 0.0_dp
    END SELECT

    scale_factors_2d = 1.0_dp

  END SUBROUTINE scale_mesh

!!!##################################################

  FUNCTION area_between_two_vectors(vect_a,vect_b)

    !###
    USE arrays
    REAL(dp),INTENT(in) :: vect_a(3),vect_b(3)
    REAL(dp) :: cross(3)
    REAL(dp) :: area_between_two_vectors

    ! area = 1/2 x magnitude of the cross-product of vectors a and b

    cross = cross_product(vect_a,vect_b)
    area_between_two_vectors = 0.5_dp * SQRT(dot_PRODUCT(cross,cross))

  END FUNCTION area_between_two_vectors


!!!##################################################

  FUNCTION area_between_three_points(point_a,point_b,point_c)

    !###
    USE arrays
    REAL(dp),INTENT(in) :: point_a(3),point_b(3),point_c(3)
    REAL(dp) :: norm(3),vect_a(3),vect_b(3)
    REAL(dp) :: area_between_three_points

    ! area = 1/2 x magnitude of the cross-product of vectors a and b

    vect_a(1:3) = point_a(1:3) - point_b(1:3)
    vect_b(1:3) = point_a(1:3) - point_c(1:3)
    norm = cross_product(vect_a,vect_b)
    area_between_three_points = 0.5_dp * SQRT(dot_PRODUCT(norm,norm))

  END FUNCTION area_between_three_points


!!! ##########################################################################

  FUNCTION ph3(I,J,K,XI)
    USE arrays
!!! dummy arguments
    INTEGER :: I,I_J_K,J,K
    REAL(dp) :: XI
!!! local variables
    REAL(dp) :: ph3

    ! K is 1,2, or 3; J is 1 or 2; I is 1 or 2

    I_J_K = 100*I + 10*J + K

    SELECT CASE(I_J_K)
    CASE(111) !i=1,j=1,k=1
       PH3=(2.0_dp*XI-3.0_dp)*XI*XI+1.0_dp  ! 2xi^3-3xi^2+1
    CASE(121) !i=1,j=2,k=1
       PH3=((XI-2.0_dp)*XI+1.0_dp)*XI      ! xi^3-2xi^2+xi
    CASE(211) !i=2,j=1,k=1
       PH3=XI*XI*(3.0_dp-2.0_dp*XI)        ! -2xi^3+3xi^2
    CASE(221) !i=2,j=2,k=1
       PH3=XI*XI*(XI-1.0_dp)              ! xi^3-xi^2
    CASE(112) !i=1,j=1,k=2
       PH3=6.0_dp*XI*(XI-1.0_dp)           ! 6xi^2-6xi
    CASE(122) !i=1,j=2,k=2
       PH3=(3.0_dp*XI-4.0_dp)*XI+1.0_dp     ! 3xi^2-4xi+1
    CASE(212) !i=2,j=1,k=2
       PH3=6.0_dp*XI*(1.0_dp-XI)           ! -6xi^2+6xi
    CASE(222) !i=2,j=2,k=2
       PH3=XI*(3.0_dp*XI-2.0_dp)           ! 3xi^2-2xi
    CASE(113) !i=1,j=1,k=3
       PH3=12.0_dp*XI-6.0_dp               ! 12xi-6
    CASE(123) !i=1,j=2,k=3
       PH3=6.0_dp*XI-4.0_dp                ! 6xi-4
    CASE(213) !i=2,j=1,k=3
       PH3=6.0_dp-12.0_dp*XI               ! -12xi+6
    CASE(223) !i=2,j=2,k=3
       PH3=6.0_dp*XI-2.0_dp                ! 6xi-2
    END SELECT

  END FUNCTION ph3

!!! ##########################################################################

  FUNCTION pl1(I,K,XI)
    USE arrays
!!! dummy arguments
    INTEGER :: I,I_K,K
    REAL(dp) :: XI
!!! local variables
    REAL(dp) :: pl1

    I_K = 10*I + K

    SELECT CASE(I_K)
    CASE(11) !i=1,k=1
       PL1=1.0_dp-XI
    CASE(21) !i=2,k=1
       PL1=XI
    CASE(12) !i=1,k=2
       PL1=-1.0_dp
    CASE(22) !i=2,k=2
       PL1=1.0_dp
    CASE(30 :) !k=3
       PL1=0.0_dp
    END SELECT

    RETURN
  END FUNCTION pl1

!!!##################################################

  FUNCTION unit_norm_to_plane_two_vectors(vect_a,vect_b)
    USE arrays
    REAL(dp),INTENT(in) :: vect_a(3),vect_b(3)
    REAL(dp) :: magnitude,norm(3)
    REAL(dp) :: unit_norm_to_plane_two_vectors(3)

    norm = cross_product(vect_a,vect_b)
    magnitude = SQRT(dot_PRODUCT(norm,norm))
    unit_norm_to_plane_two_vectors = norm/magnitude

  END FUNCTION unit_norm_to_plane_two_vectors


!!!##################################################

  FUNCTION unit_norm_to_three_points(point_a,point_b,point_c)
    USE arrays
    REAL(dp),INTENT(in) :: point_a(3),point_b(3),point_c(3)
    REAL(dp) :: magnitude,norm(3),vect_a(3),vect_b(3)
    REAL(dp) :: unit_norm_to_three_points(3)

    vect_a(1:3) = point_a(1:3) - point_b(1:3)
    vect_b(1:3) = point_a(1:3) - point_c(1:3)
    norm = cross_product(vect_a,vect_b)
    magnitude = SQRT(dot_PRODUCT(norm,norm))
    unit_norm_to_three_points = norm/magnitude

  END FUNCTION unit_norm_to_three_points

!!!##################################################

  FUNCTION angle_btwn_vectors(U,V)
    USE arrays
    !###    ANGLE calculates the angle between two vectors

    REAL(dp),INTENT(in) :: U(3),V(3)

    REAL(dp) :: ANGLE,angle_btwn_vectors,N_U(3),N_V(3)

    N_U = unit_vector(U)
    N_V = unit_vector(V)
    ANGLE = scalar_product_3(N_U,N_V)
    ANGLE = MAX(-1.0_dp,ANGLE)
    ANGLE = MIN(1.0_dp,ANGLE)
    ANGLE = COS(ANGLE)

    angle_btwn_vectors=ANGLE

  END FUNCTION angle_btwn_vectors

!!!###############################################################

  FUNCTION check_colinear_points(POINT1,POINT2,POINT3)
    USE arrays ,ONLY : dp
    !###    check_colinear_points checks whether two vectors are colinear.

    !     Parameter list
    REAL(dp) :: POINT1(3),POINT2(3),POINT3(3)
    !     Local variables
    REAL(dp) :: ERR1(3),ERR2(3),LU,LV,U(3),V(3)
    REAL(dp),PARAMETER :: zero_tol = 1.0e-14_dp
    LOGICAL :: check_colinear_points


    check_colinear_points =.FALSE.
    U(1:3)=POINT2(1:3)-POINT1(1:3)
    V(1:3)=POINT3(1:3)-POINT1(1:3)
    LU = vector_length(U)
    LV = vector_length(V)
    ! If 2 of the points are the same then LU and LV
    ! can be zero causing div by zero below and resulting in
    ! the wrong answer (on Linux)
    IF((dabs(LU)>zero_tol).AND.(dabs(LV)>zero_tol)) THEN
       ERR1(1:3)=DABS(U(1:3)/LU-V(1:3)/LV)
       ERR2(1:3)=DABS(U(1:3)/LU+V(1:3)/LV)
       IF((ERR1(1).LE.ZERO_TOL.AND.ERR1(2).LE.ZERO_TOL.AND.ERR1(3).LE. &
            ZERO_TOL).OR.(ERR2(1).LE.ZERO_TOL.AND.ERR2(2).LE.ZERO_TOL.AND. &
            ERR2(3).LE.ZERO_TOL)) check_colinear_points=.TRUE.
    ELSE
       check_colinear_points=.TRUE.
    ENDIF
  END FUNCTION check_colinear_points


!!!###############################################################

  FUNCTION cross_product(A,B)
    USE arrays
    !###  cross_product returns the vector cross product of A*B in C.

    !     Parameter List
    REAL(dp),INTENT(in) :: A(3),B(3)

    REAL(dp) :: cross_product(3)

    cross_product(1) = A(2)*B(3)-A(3)*B(2)
    cross_product(2) = A(3)*B(1)-A(1)*B(3)
    cross_product(3) = A(1)*B(2)-A(2)*B(1)

  END FUNCTION cross_product

!!!###############################################################

  FUNCTION scalar_triple_product(A,B,C)
    USE arrays
    !###  scalar_triple_product returns A.(BxC)

    !     Parameter List
    REAL(dp),INTENT(in) :: A(3),B(3),C(3)

    REAL(dp) :: scalar_triple_product

    scalar_triple_product = A(1)*(B(2)*C(3)-B(3)*C(2)) + &
         A(2)*(B(3)*C(1)-B(1)*C(3)) + A(3)*(B(1)*C(2)-B(2)*C(1))

  END FUNCTION scalar_triple_product

!!!###############################################################

  FUNCTION distance_between_points(point1, point2)
    USE arrays
    !###    calculates the distance between two arbitrary points

    REAL(dp),INTENT(in) :: point1(3),point2(3)
    INTEGER :: i
    REAL(dp) :: distance_between_points

    distance_between_points = 0.0_dp
    DO i=1,3
       distance_between_points = distance_between_points + (point1(i)-point2(i))**2
    ENDDO
    distance_between_points = dsqrt(distance_between_points)

  END FUNCTION distance_between_points

!!!###############################################################

  FUNCTION mesh_a_x_eq_b(MATRIX,VECTOR)
    USE arrays
    REAL(dp) :: MATRIX(3,3),VECTOR(3)
    !Local variables
    INTEGER :: i,j,k,pivot_row
    REAL(dp) :: A(3,4),max,pivot_value,SOLUTION(3),TEMP(4)
    REAL(dp) :: mesh_a_x_eq_b(3)


    A(1:3,1:3) = MATRIX(1:3,1:3)
    A(1:3,4) = VECTOR(1:3)
    DO k=1,2
       max=0.0_dp
       DO i=k,3
          IF(DABS(A(i,k)).GT.max)THEN
             max=DABS(A(i,k))
             pivot_row=i
          ENDIF
       ENDDO !i
       IF(pivot_row.NE.k)THEN
          DO j=1,4
             TEMP(j)=A(k,j)
             A(k,j)=A(pivot_row,j)
             A(pivot_row,j)=TEMP(j)
          ENDDO !j
       ENDIF
       pivot_value = A(k,k)
       A(k,1:4) = A(k,1:4)/pivot_value
       DO i=k+1,3
          DO j=k+1,4
             A(i,j) = A(i,j)-A(i,k)*A(k,j)
          ENDDO
          A(i,k) = 0.0_dp
       ENDDO
    ENDDO !N
    A(3,4) = A(3,4)/A(3,3)
    A(2,4) = A(2,4)-A(3,4)*A(2,3)
    A(1,4) = A(1,4)-A(3,4)*A(1,3)-A(2,4)*A(1,2)

    SOLUTION(1:3) = A(1:3,4)

    mesh_a_x_eq_b = solution

  END FUNCTION mesh_a_x_eq_b

!!!##################################################

  FUNCTION scalar_product_3(A,B)
    USE arrays
    !### calculates scalar product of two vectors A,B of length 3.

    REAL(dp),INTENT(in) :: A(*),B(*)

    INTEGER :: i
    REAL(dp) :: scalar_product_3

    scalar_product_3 = 0.0_dp
    DO i=1,3
       scalar_product_3 = scalar_product_3 + A(i)*B(i)
    ENDDO

  END FUNCTION scalar_product_3

!!!###############################################################

  FUNCTION unit_vector(A)
    USE arrays
    !###  Calculates the unit vector for an arbitrary 3x1 vector

    REAL(dp),INTENT(in) :: A(*)
    REAL(dp) :: length_a,unit_vector(3)

    length_a = vector_length(A)
    IF(length_a.GT.1.0e-6_dp)THEN
       unit_vector(1:3) = A(1:3)/length_a
    ELSE
       WRITE(*,*) ' >>WARNING: Cannot normalise a zero length vector'
       WRITE(*,*) ' We recommend debugging, but hit enter to continue'
       READ(*,*)
    ENDIF

  END FUNCTION unit_vector

!!!##################################################

  FUNCTION vector_length(A)
    USE arrays
    !###  Calculates the length of a 3x1 vector

    REAL(dp),INTENT(in) :: A(*)
    REAL(dp) :: vector_length
    INTEGER :: i

    vector_length = 0.0_dp
    DO i=1,3
       vector_length = vector_length + A(i)*A(i)
    ENDDO
    vector_length = dsqrt(vector_length)

  END FUNCTION vector_length

!!!###############################################################

  FUNCTION volume_internal_to_surface(triangles,vertex_xyz)
    USE arrays
    ! calculates the volume enclosed by a list of surface elements

    INTEGER,INTENT(in) :: triangles(:,:)
    REAL(dp),INTENT(in) :: vertex_xyz(:,:)
    REAL(dp) :: volume_internal_to_surface

!!! Local Variables
    INTEGER :: ntri,num_triangles
    REAL(dp) :: volume,V1(3),V2(3),V3(3),P4(3)

    num_triangles = COUNT(triangles(:,:).NE.0)/3

    P4 = SUM(vertex_xyz,dim=2)/SIZE(vertex_xyz,dim=2)

    volume = 0.0_dp

    DO ntri = 1,num_triangles
       V1(1:3) = P4(1:3) - vertex_xyz(1:3,triangles(1,ntri))
       V2(1:3) = P4(1:3) - vertex_xyz(1:3,triangles(2,ntri))
       V3(1:3) = P4(1:3) - vertex_xyz(1:3,triangles(3,ntri))
       volume = volume + ABS(scalar_triple_product(V1,V2,V3))
    ENDDO

    volume_internal_to_surface = volume/6.0_dp

  END FUNCTION volume_internal_to_surface

!!!###############################################################

  FUNCTION point_internal_to_surface(triangles,point_xyz,vertex_xyz)
!!! Cast a line in positive x-direction from each data point and
!!! then work out how many triangular elements it crosses. If even it is in the
!!! shape and if odd it is outside the shape
    USE arrays
    INTEGER,INTENT(in) :: triangles(:,:)
    REAL(dp),INTENT(in) :: point_xyz(3),vertex_xyz(:,:)
    LOGICAL :: point_internal_to_surface

!!! Local Variables
    INTEGER :: ntri,num_triangles
    REAL(dp) :: area,area_triangle,cofm_surfaces(3),denominator,&
         norm_v(3),point(3),P1(3),P2(3),P3(3),u
    REAL(dp),PARAMETER :: dist_tol = 1.0e-4_dp, user_tol = 1.0e-14_dp
    LOGICAL :: cross_any

    num_triangles = COUNT(triangles(:,:).NE.0)/3

    cofm_surfaces = SUM(vertex_xyz,dim=2)/SIZE(vertex_xyz,dim=2)

    ! check whether the line that joins the centre of mass of the surface mesh and the point
    ! in question crosses ANY face. If it does, then point not inside.

    cross_any = .FALSE.

    DO ntri=1,num_triangles
       P1(1:3) = vertex_xyz(1:3,triangles(1,ntri))
       P2(1:3) = vertex_xyz(1:3,triangles(2,ntri))
       P3(1:3) = vertex_xyz(1:3,triangles(3,ntri))
       norm_v = unit_norm_to_three_points(P1,P2,P3) ! unit normal to triangle plane
       ! u = (a*x1+b*y1+c*z1+d)/(a*(x1-x2)+b*(y1-y2)+c*(z1-z2))
       denominator = norm_v(1)*(point_xyz(1)-cofm_surfaces(1)) + &
            norm_v(2)*(point_xyz(2)-cofm_surfaces(2)) + &
            norm_v(3)*(point_xyz(3)-cofm_surfaces(3))
       ! denominator is zero for line parallel to plane
       IF(ABS(denominator).GT.user_tol)THEN
          ! calculate the distance of the surface point from point_xyz
          u = (dot_PRODUCT(norm_v,point_xyz)-dot_PRODUCT(norm_v,P1))/denominator
          IF(u.GE.0.0_dp.AND.u.LE.1.0_dp)THEN ! POTENTIALLY crosses. Test further (angle)
             point = point_xyz + u*(cofm_surfaces-point_xyz) ! projection to surface
             area = area_between_two_vectors(P1-point,P2-point)+ &
                  area_between_two_vectors(P1-point,P3-point)+area_between_two_vectors(P2-point,P3-point)
             area_triangle = area_between_two_vectors(P1-P2,P1-P3)
             IF(ABS(area_triangle-area).LT.dist_tol) cross_any = .TRUE.
          ENDIF
       ENDIF
    ENDDO

    IF(.NOT.cross_any)THEN
       point_internal_to_surface = .TRUE.
    ELSE
       point_internal_to_surface = .FALSE.
    ENDIF

  END FUNCTION point_internal_to_surface



END MODULE mesh_functions
