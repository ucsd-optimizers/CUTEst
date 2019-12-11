! THIS VERSION: CUTEST 1.4 - 26/02/2016 AT 08:00 GMT.

!-*-*-*-*-  C U T E S T  C I N T _ C G R D H    S U B R O U T I N E  -*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 21st August 2013

      SUBROUTINE CUTEST_Cint_cgrdh( status, n, m, X, Y, grlagf, G, jtrans,     &
                                    lj1, lj2, J_val, lh1, H_val )
      USE CUTEST
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_Bool
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( IN ) :: n, m, lj1, lj2, lh1
      INTEGER, INTENT( OUT ) :: status
      LOGICAL ( KIND = C_Bool ), INTENT( IN ) :: grlagf, jtrans
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: G
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lj1, lj2 ) :: J_val
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lh1, n ) :: H_val

!  ---------------------------------------------------------------
!  compute both the gradients of the objective, or Lagrangian, and
!  general constraint functions and the Hessian matrix of the
!  Lagrangian function of a problem initially written in Standard
!  Input Format (SIF)

!  G	 is an array which gives the value of the gradient of
!	 the objective function evaluated at X (grlagf = .FALSE.)
!        of of the Lagrangian function evaluated at X and Y
!        (GRLAGF = .TRUE.)

!  J_val is a two-dimensional array of dimension ( lj1, lj2 )
!	 which gives the value of the Jacobian matrix of the
!	 constraint functions, or its transpose, evaluated at X.
!	 If jtrans is .TRUE., the i,j-th component of the array
!        will contain the i-th derivative of the j-th constraint
!        function. Otherwise, if jtrans is .FALSE., the i,j-th
!        component of the array will contain the j-th derivative
!        of the i-th constraint function

!  H_val is a two-dimensional array which gives the value of the
!        Hessian matrix of the Lagrangian function evaluated at
!        X and Y. The i,j-th component of the array will contain
!        the derivative with respect to variables X(i) and X(j)
!  ---------------------------------------------------------------

      LOGICAL :: grlagf_fortran, jtrans_fortran

      grlagf_fortran = grlagf
      jtrans_fortran = jtrans
      CALL CUTEST_cgrdh( status, n, m, X, Y, grlagf_fortran, G,                &
                         jtrans_fortran, lj1, lj2, J_val, lh1, H_val )

      RETURN

!  end of subroutine CUTEST_Cint_cgrdh

      END SUBROUTINE CUTEST_Cint_cgrdh

!-*-*-*-*-*-*-  C U T E S T    C G R D H    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 29th December 2012

      SUBROUTINE CUTEST_cgrdh( status, n, m, X, Y, grlagf, G, jtrans,          &
                               lj1, lj2, J_val, lh1, H_val )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( IN ) :: n, m, lj1, lj2, lh1
      INTEGER, INTENT( OUT ) :: status
      LOGICAL, INTENT( IN ) :: grlagf, jtrans
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: G
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lj1, lj2 ) :: J_val
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lh1, n ) :: H_val

!  ---------------------------------------------------------------
!  compute both the gradients of the objective, or Lagrangian, and
!  general constraint functions and the Hessian matrix of the
!  Lagrangian function of a problem initially written in Standard
!  Input Format (SIF)

!  G	 is an array which gives the value of the gradient of
!	 the objective function evaluated at X (grlagf = .FALSE.)
!        of of the Lagrangian function evaluated at X and Y
!        (GRLAGF = .TRUE.)

!  J_val is a two-dimensional array of dimension ( lj1, lj2 )
!	 which gives the value of the Jacobian matrix of the
!	 constraint functions, or its transpose, evaluated at X.
!	 If jtrans is .TRUE., the i,j-th component of the array
!        will contain the i-th derivative of the j-th constraint
!        function. Otherwise, if jtrans is .FALSE., the i,j-th
!        component of the array will contain the j-th derivative
!        of the i-th constraint function

!  H_val is a two-dimensional array which gives the value of the
!        Hessian matrix of the Lagrangian function evaluated at
!        X and Y. The i,j-th component of the array will contain
!        the derivative with respect to variables X(i) and X(j)
!  ---------------------------------------------------------------

      CALL CUTEST_cgrdh_threadsafe( CUTEST_data_global,                        &
                                    CUTEST_work_global( 1 ),                   &
                                    status, n, m, X, Y, grlagf, G,             &
                                    jtrans, lj1, lj2, J_val, lh1, H_val )
      RETURN

!  end of subroutine CUTEST_cgrdh

      END SUBROUTINE CUTEST_cgrdh

!-*-*-  C U T E S T    C G R D H _ t h r e a d e d   S U B R O U T I N E  -*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 29th December 2012

      SUBROUTINE CUTEST_cgrdh_threaded( status, n, m, X, Y, grlagf, G, jtrans, &
                                        lj1, lj2, J_val, lh1, H_val, thread )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( IN ) :: n, m, lj1, lj2, lh1, thread
      INTEGER, INTENT( OUT ) :: status
      LOGICAL, INTENT( IN ) :: grlagf, jtrans
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: G
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lj1, lj2 ) :: J_val
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lh1, n ) :: H_val

!  ---------------------------------------------------------------
!  compute both the gradients of the objective, or Lagrangian, and
!  general constraint functions and the Hessian matrix of the
!  Lagrangian function of a problem initially written in Standard
!  Input Format (SIF)

!  G	 is an array which gives the value of the gradient of
!	 the objective function evaluated at X (grlagf = .FALSE.)
!        of of the Lagrangian function evaluated at X and Y
!        (GRLAGF = .TRUE.)

!  J_val is a two-dimensional array of dimension ( lj1, lj2 )
!	 which gives the value of the Jacobian matrix of the
!	 constraint functions, or its transpose, evaluated at X.
!	 If jtrans is .TRUE., the i,j-th component of the array
!        will contain the i-th derivative of the j-th constraint
!        function. Otherwise, if jtrans is .FALSE., the i,j-th
!        component of the array will contain the j-th derivative
!        of the i-th constraint function

!  H_val is a two-dimensional array which gives the value of the
!        Hessian matrix of the Lagrangian function evaluated at
!        X and Y. The i,j-th component of the array will contain
!        the derivative with respect to variables X(i) and X(j)
!  ---------------------------------------------------------------

!  check that the specified thread is within range

      IF ( thread < 1 .OR. thread > CUTEST_data_global%threads ) THEN
        IF ( CUTEST_data_global%out > 0 )                                      &
          WRITE( CUTEST_data_global%out, "( ' ** CUTEST error: thread ', I0,   &
         &  ' out of range [1,', I0, ']' )" ) thread, CUTEST_data_global%threads
        status = 4 ; RETURN
      END IF

!  evaluate using specified thread

      CALL CUTEST_cgrdh_threadsafe( CUTEST_data_global,                        &
                                    CUTEST_work_global( thread ),              &
                                    status, n, m, X, Y, grlagf, G,             &
                                    jtrans, lj1, lj2, J_val, lh1, H_val )
      RETURN

!  end of subroutine CUTEST_cgrdh_threaded

      END SUBROUTINE CUTEST_cgrdh_threaded

!-*-  C U T E S T    C G R D H _ t h r e a d s a f e   S U B R O U T I N E  -*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, November 1991
!   fortran 2003 version released in CUTEst, 24th November 2012

      SUBROUTINE CUTEST_cgrdh_threadsafe( data, work, status, n, m, X, Y,      &
                                          grlagf, G, jtrans, lj1, lj2, J_val,  &
                                          lh1, H_val )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      TYPE ( CUTEST_work_type ), INTENT( INOUT ) :: work
      INTEGER, INTENT( IN ) :: n, m, lj1, lj2, lh1
      INTEGER, INTENT( OUT ) :: status
      LOGICAL, INTENT( IN ) :: grlagf, jtrans
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: G
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lj1, lj2 ) :: J_val
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lh1, n ) :: H_val

!  ---------------------------------------------------------------
!  compute both the gradients of the objective, or Lagrangian, and
!  general constraint functions and the Hessian matrix of the
!  Lagrangian function of a problem initially written in Standard
!  Input Format (SIF)

!  G	 is an array which gives the value of the gradient of
!	 the objective function evaluated at X (grlagf = .FALSE.)
!        of of the Lagrangian function evaluated at X and Y
!        (GRLAGF = .TRUE.)

!  J_val is a two-dimensional array of dimension ( lj1, lj2 )
!	 which gives the value of the Jacobian matrix of the
!	 constraint functions, or its transpose, evaluated at X.
!	 If jtrans is .TRUE., the i,j-th component of the array
!        will contain the i-th derivative of the j-th constraint
!        function. Otherwise, if jtrans is .FALSE., the i,j-th
!        component of the array will contain the j-th derivative
!        of the i-th constraint function

!  H_val is a two-dimensional array which gives the value of the
!        Hessian matrix of the Lagrangian function evaluated at
!        X and Y. The i,j-th component of the array will contain
!        the derivative with respect to variables X(i) and X(j)
!  ---------------------------------------------------------------

!  local variables

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, jj, ll, nnzh, iendgv, icon
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, ifstat, igstat, alloc_status
      REAL ( KIND = wp ) :: ftt, gi, gii, scalee
      REAL :: time_in, time_out
      LOGICAL :: nontrv
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )
      EXTERNAL :: RANGE

      IF ( work%record_times ) CALL CPU_TIME( time_in )

!  check input parameters

      IF ( data%numcon > 0 ) THEN
        IF ( jtrans ) THEN
          IF ( lj1 < n .OR. lj2 < m ) THEN
            IF ( lj1 < n .AND. data%out > 0 ) WRITE( data%out, 2020 ) n
            IF ( lj2 < m .AND. data%out > 0 ) WRITE( data%out, 2030 ) m
            status = 2 ; GO TO 990
          END IF
        ELSE
          IF ( lj1 < m .OR. lj2 < n ) THEN
            IF ( lj1 < m .AND. data%out > 0 ) WRITE( data%out, 2020 ) m
            IF ( lj2 < n .AND. data%out > 0 ) WRITE( data%out, 2030 ) n
            status = 2 ; GO TO 990
          END IF
        END IF
      END IF
      IF ( lh1 < n ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2040 ) n
        status = 2 ; GO TO 990
      END IF

!  there are non-trivial group functions

      DO i = 1, MAX( data%nel, data%ng )
        work%ICALCF( i ) = i
      END DO

!  evaluate the element function values

      CALL ELFUN( work%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, work%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function values

      CALL ELFUN( work%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, work%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  3, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the group argument values ft

      DO ig = 1, data%ng
        ftt = - data%B( ig )

!  include the contribution from the linear element

        DO j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
          ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
        END DO

!  include the contributions from the nonlinear elements

        DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
          ftt = ftt + data%ESCALE( j ) * work%FUVALS( data%IELING( j ) )
        END DO
        work%FT( ig ) = ftt

!  record the derivatives of trivial groups

        IF ( data%GXEQX( ig ) ) THEN
          work%GVALS( ig, 2 ) = 1.0_wp
          work%GVALS( ig, 3 ) = 0.0_wp
        END IF
      END DO

!  evaluate the group derivative values

      IF ( .NOT. data%altriv ) THEN
        CALL GROUP( work%GVALS, data%ng, work%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, work%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                   .TRUE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  for unconstrained problems, skip construction of work%W_ws( ig ) and skip
!  specialized construction of gradient and Jacobian.  Call ELGRD instead

!  change the group weightings to include the contributions from the
!  Lagrange multipliers

      IF ( data%numcon > 0 ) THEN
        DO ig = 1, data%ng
          i = data%KNDOFC( ig )
          IF ( i == 0 ) THEN
            work%GSCALE_used( ig ) = data%GSCALE( ig )
          ELSE
            work%GSCALE_used( ig ) = data%GSCALE( ig ) * Y( i )
          END IF
        END DO

!  compute the gradient values. Initialize the gradient and Jacobian (or its
!  transpose) as zero

        G( : n ) = 0.0_wp
        IF ( jtrans ) THEN
          J_val( : n, : m ) = 0.0_wp
        ELSE
          J_val( : m, : n ) = 0.0_wp
        END IF

!  consider the ig-th group

        DO ig = 1, data%ng
          ig1 = ig + 1
          icon = data%KNDOFC( ig )
          istrgv = data%ISTAGV( ig )
          iendgv = data%ISTAGV( ig1 ) - 1
          nelow = data%ISTADG( ig )
          nelup = data%ISTADG( ig1 ) - 1
          nontrv = .NOT. data%GXEQX( ig )

!  compute the first derivative of the group

          gi = data%GSCALE( ig )
          gii = work%GSCALE_used( ig )
          IF  ( nontrv ) THEN
            gi = gi * work%GVALS( ig, 2 )
            gii = gii * work%GVALS( ig, 2 )
          END IF

!  this is the first gradient evaluation or the group has nonlinear elements

          IF ( work%firstg .OR. nelow <= nelup ) THEN
            work%W_ws( data%ISVGRP( istrgv : iendgv ) ) = 0.0_wp

!  loop over the group's nonlinear elements

             DO ii = nelow, nelup
               iel = data%IELING( ii )
               k = data%INTVAR( iel )
               l = data%ISTAEV( iel )
               nvarel = data%ISTAEV( iel + 1 ) - l
               scalee = data%ESCALE( ii )
               IF ( data%INTREP( iel ) ) THEN

!  the iel-th element has an internal representation

                 nin = data%INTVAR( iel + 1 ) - k
                 CALL RANGE( iel, .TRUE., work%FUVALS( k ),                    &
                             work%W_el, nvarel, nin, data%ITYPEE( iel ),       &
                             nin, nvarel )
!DIR$ IVDEP
                 DO i = 1, nvarel
                   j = data%IELVAR( l )
                   work%W_ws( j ) = work%W_ws( j ) + scalee * work%W_el( i )
                   l = l + 1
                 END DO
               ELSE

!  the iel-th element has no internal representation

!DIR$ IVDEP
                 DO i = 1, nvarel
                   j = data%IELVAR( l )
                   work%W_ws( j ) = work%W_ws( j ) + scalee * work%FUVALS( k )
                   k = k + 1 ; l = l + 1
                 END DO
               END IF
             END DO

!  include the contribution from the linear element

!DIR$ IVDEP
            DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
              j = data%ICNA( k )
              work%W_ws( j ) = work%W_ws( j ) + data%A( k )
            END DO

!  Allocate a gradient

!DIR$ IVDEP
             DO i = istrgv, iendgv
               ll = data%ISVGRP( i )

!  the group belongs to the objective function

               IF ( icon == 0 ) THEN
                 G( ll ) = G( ll ) + gi * work%W_ws( ll )

!  the group defines a constraint

               ELSE
                 IF ( JTRANS ) THEN
                   J_val( ll, icon ) = gi * work%W_ws( ll )
                 ELSE
                   J_val( icon, ll ) = gi * work%W_ws( ll )
                 END IF
                 IF ( grlagf ) G( ll ) = G( ll ) + gii * work%W_ws( ll )
               END IF

!  if the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC

               IF ( nontrv ) THEN
                 jj = work%ISTAJC( ll )
                 work%FUVALS( data%lgrjac + jj ) = work%W_ws( ll )

!  increment the address for the next nonzero in the column of the Jacobian
!  for variable ll

                 work%ISTAJC( ll ) = jj + 1
               END IF
             END DO

!  This is not the first gradient evaluation and there is only a linear element

          ELSE

!  allocate a gradient

!DIR$ IVDEP
            DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
              ll = data%ICNA( k )

!  the group belongs to the objective function

              IF ( icon == 0 ) THEN
                G( ll ) = G( ll ) + gi * data%A( k )

!  the group defines a constraint

              ELSE
                IF ( JTRANS ) THEN
                  J_val( ll, icon ) = gi * data%A( k )
                ELSE
                  J_val( icon, ll ) = gi * data%A( k )
                END IF
                IF ( grlagf ) G( ll ) = G( ll ) + gii * data%A( k )
              END IF
            END DO

!  the group is non-trivial; increment the starting addresses for the groups
!  used by each variable in the (unchanged) linear element to avoid resetting
!  the nonzeros in the Jacobian

            IF ( nontrv ) THEN
!DIR$ IVDEP
              DO i = istrgv, iendgv
                ll = data%ISVGRP( i )
                work%ISTAJC( ll ) = work%ISTAJC( ll ) + 1
              END DO
            END IF
          END IF
        END DO

!  reset the starting addresses for the lists of groups using each variable to
!  their values on entry

        DO i = n, 2, - 1
           work%ISTAJC( i ) = work%ISTAJC( i - 1 )
        END DO
        work%ISTAJC( 1 ) = 1

!  compute the gradient value

      ELSE
        CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,         &
               data%nvrels, data%nnza, data%nvargp, work%firstg, data%ICNA,    &
               data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,             &
               data%IELVAR, data%INTVAR, data%A, work%GVALS( : , 2 ),          &
               work%FUVALS, data%lnguvl, work%FUVALS( data%lggfx + 1 ),        &
               data%GSCALE, data%ESCALE, work%FUVALS( data%lgrjac + 1 ),       &
               data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE, &
               work%ISTAJC, work%W_ws, work%W_el, RANGE )

!  store the gradient value

        DO i = 1, n
          G( i ) = work%FUVALS( data%lggfx + i )
        END DO
      END IF
      work%firstg = .FALSE.

!  assemble the Hessian

      IF ( data%numcon > 0 ) THEN
        CALL CUTEST_assemble_hessian(                                          &
               n, data%ng, data%nel, data%ntotel, data%nvrels, data%nnza,      &
               data%maxsel, data%nvargp, data%ISTADH,                          &
               data%ICNA, data%ISTADA, data%INTVAR, data%IELVAR, data%IELING,  &
               data%ISTADG, data%ISTAEV, data%ISTAGV, data%ISVGRP, data%A,     &
               work%FUVALS, data%lnguvl, work%FUVALS, data%lnhuvl,             &
               work%GVALS( : , 2 ), work%GVALS( :  , 3 ), work%GSCALE_used,    &
               data%ESCALE, data%GXEQX, data%ITYPEE, data%INTREP, RANGE,       &
               0, data%out, data%out, .TRUE., .FALSE.,                         &
               n, status, alloc_status, bad_alloc,                             &
               work%array_status, work%lh_row, work%lh_col, work%lh_val,       &
               work%H_row, work%H_col, work%H_val, work%ROW_start,             &
               work%POS_in_H, work%USED, work%FILLED,                          &
               work%lrowst, work%lpos, work%lused, work%lfilled,               &
               work%W_ws, work%W_el, work%W_in, work%H_el, work%H_in,          &
               nnzh = nnzh )
      ELSE
        CALL CUTEST_assemble_hessian(                                          &
               n, data%ng, data%nel, data%ntotel, data%nvrels, data%nnza,      &
               data%maxsel, data%nvargp, data%ISTADH,                          &
               data%ICNA, data%ISTADA, data%INTVAR, data%IELVAR, data%IELING,  &
               data%ISTADG, data%ISTAEV, data%ISTAGV, data%ISVGRP, data%A,     &
               work%FUVALS, data%lnguvl, work%FUVALS, data%lnhuvl,             &
               work%GVALS( : , 2 ), work%GVALS( :  , 3 ), data%GSCALE,         &
               data%ESCALE, data%GXEQX, data%ITYPEE, data%INTREP, RANGE,       &
               0, data%out, data%out, .TRUE., .FALSE.,                         &
               n, status, alloc_status, bad_alloc,                             &
               work%array_status, work%lh_row, work%lh_col, work%lh_val,       &
               work%H_row, work%H_col, work%H_val, work%ROW_start,             &
               work%POS_in_H, work%USED, work%FILLED,                          &
               work%lrowst, work%lpos, work%lused, work%lfilled,               &
               work%W_ws, work%W_el, work%W_in, work%H_el, work%H_in,          &
               nnzh = nnzh )
      END IF

!  check for errors in the assembly

      IF ( status > 0 ) GO TO 990

!  initialize the dense Hessian matrix

      H_val( : n, : n ) = 0.0_wp

!  transfer the matrix from co-ordinate to dense storage and symmetrize the
!  martix

      DO k = 1, nnzh
        i = work%H_row( k ) ; j = work%H_col( k )
        H_val( i, j ) = work%H_val( k ) ; H_val( j, i ) = work%H_val( k )
      END DO

!  update the counters for the report tool

      work%nc2og = work%nc2og + 1
      work%nc2cg = work%nc2cg + work%pnc
      work%nc2oh = work%nc2oh + 1
      work%nc2ch = work%nc2ch + work%pnc
      status = 0
      GO TO 990

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CGRDH: error flag raised during SIF evaluation' )" )
      status = 3

!  update elapsed CPU time if required

  990 CONTINUE
      IF ( work%record_times ) THEN
        CALL CPU_TIME( time_out )
        work%time_cgrdh = work%time_cgrdh + time_out - time_in
      END IF
      RETURN

!  non-executable statements

 2020 FORMAT( ' ** SUBROUTINE CGRDH: Increase the leading dimension of J_val', &
              ' to ', I0 )
 2030 FORMAT( ' ** SUBROUTINE CGRDH: Increase the second dimension of J_val',  &
              ' to ', I0 )
 2040 FORMAT( ' ** SUBROUTINE CGRDH: Increase the leading dimension of H_val', &
              ' to ', I0 )

!  end of subroutine CUTEST_cgrdh_threadsafe

      END SUBROUTINE CUTEST_cgrdh_threadsafe
