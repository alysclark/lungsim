MODULE math_utilities
  !*Brief Description:* This module contains solvers required for lung problems
  !
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  !
  !
  USE arrays, ONLY: dp
  USE diagnostics, ONLY: enter_exit

  IMPLICIT NONE
  PRIVATE

  PUBLIC ax_cr,diagonal_pointer_cr,ilu_cr,lus_cr,mult_givens,rearrange_cr,bessel_complex
  PUBLIC sort_integer_list
  PUBLIC sort_real_list


CONTAINS

  SUBROUTINE bessel_complex(z,bessel0,bessel1)
    USE other_consts,ONLY: PI
    USE arrays, ONLY:dp

    COMPLEX(dp), INTENT(in) :: z
    COMPLEX(dp), INTENT(out) :: bessel0,bessel1

    REAL(dp) :: a(12),a1(10),b(12)
    REAL(dp) :: a0
    COMPLEX(dp) :: ci,cr,z1,ca,zr
    INTEGER :: k,k0
    !  complex ( kind = 8 ) ca
    !  complex ( kind = 8 ) cb
    !  complex ( kind = 8 ) ci
    !  complex ( kind = 8 ) cr
    !  complex ( kind = 8 ) cs
    !  complex ( kind = 8 ) ct
    !  complex ( kind = 8 ) cw
    !  integer ( kind = 4 ) k
    !  integer ( kind = 4 ) k0
    !  real ( kind = 8 ) pi
    !  real ( kind = 8 ) w0
    !  complex ( kind = 8 ) z
    !  complex ( kind = 8 ) z1
    !  complex ( kind = 8 ) z2
    !  complex ( kind = 8 ) zr
    !  complex ( kind = 8 ) zr2
    a = (/ &
         0.125e00_dp,           7.03125e-02_dp,&
         7.32421875e-02_dp,      1.1215209960938e-01_dp,&
         2.2710800170898e-01_dp, 5.7250142097473e-01_dp,&
         1.7277275025845e00_dp, 6.0740420012735e00_dp,&
         2.4380529699556e01_dp, 1.1001714026925e02_dp,&
         5.5133589612202e02_dp, 3.0380905109224e03_dp /)
    a1 = (/ &
         0.125e00_dp,            0.2109375e00_dp, &
         1.0986328125e00_dp,     1.1775970458984e01_dp, &
         2.1461706161499e002_dp, 5.9511522710323e03_dp, &
         2.3347645606175e05_dp,  1.2312234987631e07_dp, &
         8.401390346421e08_dp,   7.2031420482627e10_dp /)
    b = (/ &
         -0.375e00_dp,           -1.171875e-01_dp, &
         -1.025390625e-01_dp,     -1.4419555664063e-01_dp, &
         -2.7757644653320e-01_dp, -6.7659258842468e-01_dp, &
         -1.9935317337513e00_dp, -6.8839142681099e00_dp, &
         -2.7248827311269e01_dp, -1.2159789187654e02_dp, &
         -6.0384407670507e02_dp, -3.3022722944809e03_dp /)

    !
    ci = CMPLX (0.0_dp,1.0_dp,8)
    a0 = ABS (z)
    z1 = z
    !
    IF(a0.EQ.0.0_dp)THEN
       bessel0 = CMPLX(1.0_dp,0.0_dp,8)
       bessel1 = CMPLX(0.0_dp,0.0_dp,8)
    ENDIF

    IF(REAL(z).LT.0.0_dp) THEN
       z1 = -z
    ENDIF
    !
    IF( a0 <= 18.0D+00 ) THEN

       bessel0 =CMPLX(1.0_dp,0.0_dp,8)
       cr = CMPLX(1.0_dp,0.0_dp,8)
       DO k = 1,50
          cr = 0.25_dp*cr* z1**2/k**2
          bessel0 = bessel0+cr
          IF (ABS (cr/bessel0).LT.1.0e-15) THEN
             EXIT
          ENDIF
       ENDDO

       bessel1 =CMPLX(1.0_dp,0.0_dp,8)
       cr = CMPLX(1.0_dp,0.0_dp,8)
       DO k = 1,50
          cr = 0.25_dp*cr*z**2/(k*(k+1))
          bessel1 = bessel1+cr
          IF (ABS (cr/bessel1).LT.1.0e-15) THEN
             EXIT
          ENDIF
       ENDDO

       bessel1 = 0.5_dp*z1*bessel1

    ELSE

       IF ( a0 < 35.0D+00 ) THEN
          k0 = 12
       ELSE IF ( a0 < 50.0D+00 ) THEN
          k0 = 9
       ELSE
          k0 = 7
       END IF

       ca = EXP(z1)/SQRT(2.0_dp*pi*z1)
       bessel0 = CMPLX(1.0_dp,0.0_dp,8)
       zr = 1.0_dp/z1
       DO k = 1,k0
          bessel0 = bessel0 + a(k) * zr ** k
       ENDDO
       bessel0 = ca * bessel0
       bessel1 = CMPLX (1.0_dp,0.0_dp,8)
       DO k = 1,k0
          bessel1 =bessel1+b(k)*zr**k
       END DO
       bessel1 = ca * bessel1
    ENDIF


    IF ( REAL (z).LT.0.0_dp)THEN
       bessel1 = - bessel1
    ENDIF

  END SUBROUTINE bessel_complex


  !
  !###########################################################################
  !
  !*ax_cr:* Computes A*x for a matrix stored in sparse compressed row form
  SUBROUTINE ax_cr ( n, ia, ja, a, x, w )
    IMPLICIT NONE

    INTEGER ( kind = 4 ) n !the order of the system
    INTEGER ( kind = 4 ) ia(*) !ia(n+1) row indices
    INTEGER ( kind = 4 ) ja(*) !ja(nz_num) column indices
    REAL ( kind = 8 ) a(*) !a(nz_num) Matrix values
    REAL ( kind = 8 ) x(*) !x(n) Vector to be multiplied by A
    REAL ( kind = 8 ) w(*) !w(n) Value of A*x

    INTEGER ( kind = 4 ) i
    INTEGER ( kind = 4 ) k1
    INTEGER ( kind = 4 ) k2

    w(1:n) = 0.0D+00

    DO i = 1, n
       k1 = ia(i)
       k2 = ia(i+1) - 1
       w(i) = w(i) + dot_PRODUCT ( a(k1:k2), x(ja(k1:k2)) )
    END DO

    RETURN
  END SUBROUTINE ax_cr
  !
  !##############################################################################
  !
  ! *ILU_CR:* computes the incomplete LU factorization of a matrix. For a matrix
  ! stored in compressed row format.
  !    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
  !    of each row.
  !    Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.
  SUBROUTINE ilu_cr ( n, nz_num, ia, ja, a, ua, l )
    INTEGER ( kind = 4 ) n
    INTEGER ( kind = 4 ) nz_num
    INTEGER ( kind = 4 ) ia(*) !ia(n+1)
    INTEGER ( kind = 4 ) ja(*) !ja(nz_num)
    REAL ( kind = 8 ) a(*) !a(nz_num)
    INTEGER ( kind = 4 ) ua(*) !ua(n)
    REAL ( kind = 8 ) l(*) !l(nz_num)

    INTEGER ( kind = 4 ) i
    INTEGER ( kind = 4 ) iw(n)
    INTEGER ( kind = 4 ) j
    INTEGER ( kind = 4 ) jj
    INTEGER ( kind = 4 ) jrow
    INTEGER ( kind = 4 ) jw
    INTEGER ( kind = 4 ) k
    REAL ( kind = 8 ) tl


    !  Copy A.
    l(1:nz_num) = a(1:nz_num)

    DO i = 1, n ! for each row, up to max number of rows
       !  IW points to the nonzero entries in row I.
       iw(1:n) = -1
       DO k = ia(i), ia(i+1) - 1 !for each
          iw(ja(k)) = k
       END DO
       DO j = ia(i), ia(i+1) - 1
          jrow = ja(j)
          IF ( i <= jrow ) THEN
             EXIT
          END IF
          tl = l(j) * l(ua(jrow))
          l(j) = tl
          DO jj = ua(jrow) + 1, ia(jrow+1) - 1
             jw = iw(ja(jj))
             IF ( jw /= -1 ) THEN
                l(jw) = l(jw) - tl * l(jj)
             END IF
          END DO
       END DO
       ua(i) = j
       IF ( jrow /= i ) THEN
          WRITE ( *, '(a)' ) ' '
          WRITE ( *, '(a)' ) 'ILU_CR - Fatal error!'
          WRITE ( *, '(a)' ) '  JROW ~= I'
          WRITE ( *, '(a,i8)' ) '  JROW = ', jrow
          WRITE ( *, '(a,i8)' ) '  I    = ', i
          STOP
       END IF
       IF ( l(j) == 0.0D+00 ) THEN
          WRITE ( *, '(a)' ) ' '
          WRITE ( *, '(a)' ) 'ILU_CR - Fatal error!'
          WRITE ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
          WRITE ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
          STOP
       END IF
       l(j) = 1.0D+00 / l(j)
    END DO

    l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

    RETURN
  END SUBROUTINE ilu_cr
  !
  !##############################################################################
  !
  !*DIAGONAL_POINTER_CR:* finds diagonal entries in a sparse compressed row matrix.
  !    The array UA can be used to locate the diagonal elements of the matrix.
  !    It is assumed that every row of the matrix includes a diagonal element,
  !    and that the elements of each row have been ascending sorted.
  SUBROUTINE diagonal_pointer_cr ( n, ia, ja, ua )
    INTEGER ( kind = 4 ) n
    INTEGER ( kind = 4 ) ia(*) !ia(n+1)
    INTEGER ( kind = 4 ) ja(*) !ja(nz_num)
    INTEGER ( kind = 4 ) ua(*) !ua(n)

    INTEGER ( kind = 4 ) i
    INTEGER ( kind = 4 ) k

    ua(1:n) = -1

    DO i = 1, n
       DO k = ia(i), ia(i+1) - 1
          IF ( ja(k) == i ) THEN
             ua(i) = k
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE diagonal_pointer_cr

  !*****************************************************************************80

  SUBROUTINE lus_cr ( n, ia, ja, l, ua, r, z )
!!! LUS_CR applies the incomplete LU preconditioner.
    !    The linear system M * Z = R is solved for Z.  M is the incomplete
    !    LU preconditioner matrix, and R is a vector supplied by the user.
    !    So essentially, we're solving L * U * Z = R.
    !    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
    !    of each row.
    !    Input, real ( kind = 8 ) R(N), the right hand side.
    !    Output, real ( kind = 8 ) Z(N), the solution of the system M * Z = R.
    IMPLICIT NONE

    INTEGER ( kind = 4 ) n
    INTEGER ( kind = 4 ) ia(*) !ia(n+1)
    INTEGER ( kind = 4 ) ja(*) !ja(nz_num)
    REAL ( kind = 8 ) l(*) !l(nz_num)
    INTEGER ( kind = 4 ) ua(*) !ua(n)
    REAL ( kind = 8 ) r(*) !r(n)

    INTEGER ( kind = 4 ) i
    INTEGER ( kind = 4 ) j
    REAL ( kind = 8 ) w(n)
    REAL ( kind = 8 ) z(n)

    !  Copy R in.
    w(1:n) = r(1:n)

    !  Solve L * w = w where L is unit lower triangular.
    DO i = 2, n
       DO j = ia(i), ua(i) - 1
          w(i) = w(i) - l(j) * w(ja(j))
       END DO
    END DO

    !  Solve U * w = w, where U is upper triangular.
    DO i = n, 1, -1
       DO j = ua(i) + 1, ia(i+1) - 1
          w(i) = w(i) - l(j) * w(ja(j))
       END DO
       w(i) = w(i) / l(ua(i))
    END DO

    !  Copy Z out.
    z(1:n) = w(1:n)

    RETURN
  END SUBROUTINE lus_cr

  !*****************************************************************************80
  SUBROUTINE mult_givens ( c, s, k, g )
!!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
    !    In order to make it easier to compare this code with the Original C,
    !    the vector indexing is 0-based.
    !    Input, real ( kind = 8 ) C, S, the cosine and sine of a Givens
    !    rotation.
    !
    !    Input, integer ( kind = 4 ) K, indicates the location of the first
    !    vector entry.
    !
    !    Input/output, real ( kind = 8 ) G(1:K+1), the vector to be modified.
    !    On output, the Givens rotation has been applied to entries G(K) and G(K+1).

    IMPLICIT NONE

    REAL ( kind = 8 ) c
    REAL ( kind = 8 ) s
    INTEGER ( kind = 4 ) k
    REAL ( kind = 8 ) g(*) !g(1:k+1)

    REAL ( kind = 8 ) g1
    REAL ( kind = 8 ) g2

    g1 = c * g(k) - s * g(k+1)
    g2 = s * g(k) + c * g(k+1)

    g(k)   = g1
    g(k+1) = g2

    RETURN
  END SUBROUTINE mult_givens


  !*****************************************************************************80
  SUBROUTINE rearrange_cr ( n, ia, ja, a )
!!! REARRANGE_CR sorts a sparse compressed row matrix.
    !    This routine guarantees that the entries in the CR matrix
    !    are properly sorted.
    !
    !    After the sorting, the entries of the matrix are rearranged in such
    !    a way that the entries of each column are listed in ascending order
    !    of their column values.
    !    Input, integer ( kind = 4 ) N, the order of the system.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
    !
    !    Input, integer ( kind = 4 ) IA(N+1), the compressed row indices.
    !
    !    Input/output, integer ( kind = 4 ) JA(NZ_NUM), the column indices.
    !    On output, these may have been rearranged by the sorting.
    !
    !    Input/output, real ( kind = 8 ) A(NZ_NUM), the matrix values.  On output,
    !    the matrix values may have been moved somewhat because of the sorting.
    !
    IMPLICIT NONE

    INTEGER ( kind = 4 ) n
    INTEGER ( kind = 4 ) ia(*) !ia(n+1)
    INTEGER ( kind = 4 ) ja(*) !ja(nz_num)
    REAL ( kind = 8 ) a(*) !a(nz_num)

    INTEGER ( kind = 4 ) i
    INTEGER ( kind = 4 ) i4temp
    INTEGER ( kind = 4 ) k
    INTEGER ( kind = 4 ) l
    REAL ( kind = 8 ) r8temp

    DO i = 1, n

       DO k = ia(i), ia(i+1) - 2
          DO l = k + 1, ia(i+1) - 1

             IF ( ja(l) < ja(k) ) THEN
                i4temp = ja(l)
                ja(l)  = ja(k)
                ja(k)  = i4temp

                r8temp = a(l)
                a(l)   = a(k)
                a(k)   = r8temp
             END IF

          END DO
       END DO

    END DO

    RETURN
  END SUBROUTINE rearrange_cr


  !######################################################################

  !
  !*sort_integer_list:* sorts a list of integer values into a non-decreasing order.
  ! sorts N integer IDATA values into a non-decreasing sequence using IHEAPSORT
  ! (N>50) or ISHELLSORT (N>50) and then  removes all duplicates from the list. On
  ! exit N contains the number of unique elements in the list.
  !
  SUBROUTINE sort_integer_list(N,IDATA)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SORT_INTEGER_LIST" :: SORT_INTEGER_LIST
    INTEGER :: IDATA(:),N

    !Local Variables
    INTEGER :: count1,index,itemp,nolist
    LOGICAL :: CONTINUE

    CHARACTER(len=60) :: sub_name

    sub_name = 'sort_integer_list'
    CALL enter_exit(sub_name,1)

    !order the array non-decreasing
    DO nolist=2,N
       count1=0
       CONTINUE=.TRUE.
       DO WHILE(CONTINUE)
          IF(IDATA(nolist-count1).LT.IDATA(nolist-count1-1))THEN
             itemp=IDATA(nolist-1)
             IDATA(nolist-1)=IDATA(nolist)
             IDATA(nolist)=itemp
             count1=count1+1
             IF(nolist-count1-1.EQ.0) CONTINUE=.FALSE.
          ELSE
             CONTINUE=.FALSE.
          ENDIF
       ENDDO !while
    ENDDO !N

    !eliminate duplicate entries
    index=0
    DO nolist=2,N
       IF(IDATA(nolist).EQ.IDATA(nolist-1)) THEN
          index=index+1
       ELSE
          IDATA(nolist-index)=IDATA(nolist)
       ENDIF
    ENDDO !nolist

    N=N-index

    CALL enter_exit(sub_name,2)

  END SUBROUTINE sort_integer_list

!!!#########################################################################
  !*sort_real_list:* sorts a list of real values into a non-decreasing order
  ! using a bubble sort algorithm.

  SUBROUTINE sort_real_list(n,RDATA,INDEX)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SORT_REAL_LIST" :: SORT_REAL_LIST

    INTEGER :: INDEX(*),n
    REAL(dp) :: RDATA(*)

    !Local Variables
    INTEGER :: FLAG,i,ITEMP,j,k
    REAL(dp) :: TEMP

    CHARACTER(len=60) :: sub_name

    sub_name = 'sort_real_list'
    CALL enter_exit(sub_name,1)

    IF(N.LE.1) THEN
    ELSE
       FLAG=n
       DO i=1,n
          k=FLAG-1
          FLAG=0
          DO j=1,k
             IF(RDATA(j).GT.RDATA(j+1)) THEN
                TEMP=RDATA(j)
                RDATA(j)=RDATA(j+1)
                RDATA(j+1)=TEMP
                ITEMP=INDEX(j)
                INDEX(j)=INDEX(j+1)
                INDEX(j+1)=ITEMP
                FLAG=j
             ENDIF
          ENDDO
          IF(FLAG.EQ.0) THEN
             WRITE(*,*) 'warning in rsort'
          ENDIF
       ENDDO
    ENDIF

    CALL enter_exit(sub_name,2)

  END SUBROUTINE sort_real_list


END MODULE math_utilities
