module DNOPT_mod
  logical           :: unconstrained
end module DNOPT_mod

program           DNOPT_main

  use DNOPT_mod
  implicit none

  character*10      :: pName

  integer           :: lenrw, leniw, lencw, minrw, miniw, mincw
  integer           :: start, m, n, nb, mLCon, mNCon, nnJac, nnObj, &
                       iObjA, ldA, ldcj, ldH, ldJCon,               &
                       inform, nInf, nNames, nFirst
  integer           :: Errors, i, ic, iP, iS, j, jslack, neq, status

  double precision  :: objAdd, Obj, sInf, fObj, ObjNP, normG
  double precision  :: cpu(2), calls(7)
  logical           :: secondDeriv


  integer,          allocatable :: iu(:), iw(:), iw0(:)
  double precision, allocatable :: ru(:), rw(:), rw0(:)
  character*8,      allocatable :: cu(:), cw(:), cw0(:)

  integer,          allocatable :: state(:)
  double precision, allocatable :: x0(:), x(:), bl(:), bu(:), y(:), &
                                   c(:), H(:,:), fCon(:), gObj(:),  &
                                   A(:,:), Jcon(:,:), cJac(:,:)

  character*8,      allocatable :: Names(:)
  character*10,     allocatable :: vname(:), cname(:)
  logical,          allocatable :: equatn(:), linear(:)


  integer,          parameter   :: iPrint = 15, iSpecs = 4,    &
                                   iSumm  = 0,  nout   = 6,    &
                                   input  = 55, io_buffer = 11

  integer,          parameter   :: lvlHes = 86

  double precision, parameter   :: zero = 0.0d+0, big = 1.0d+20

  external                      :: DNOPT_evalcj, DNOPT_evalfg, DNOPT_evalh, DNOPT_hv


  status = 0
  inform = 0

  ! -----------------------
  ! Open problem data file.
  ! -----------------------
  open ( input, file = 'OUTSDIF.d', form = 'formatted', status = 'old' )
  rewind ( input )


  ! ---------------------------
  ! Compute problem dimensions.
  ! ---------------------------
  call CUTEST_cdimen( status, input, n, m )
  if (status /= 0) go to 920

  unconstrained = m .eq. 0

  ! -------------------------------
  ! Allocate space for the problem.
  ! -------------------------------
  nb     = n+m+1
  ldA    = m+1 !mLCon
  ldH    = nb
  ldJcon = m   !mNCon
  ldcJ   = m

  allocate ( state(nb), names(nb), x0(nb), x(nb), y(nb),    &
             bl(nb), bu(nb), gObj(n), fCon(ldjCon),         &
             A(ldA,n), Jcon(ldJcon,n), cJac(m,n), H(ldH,n), &
             STAT=status )

  allocate( c(m+1), equatn(m+1), linear(m+1), vname(n), cname(m+1), STAT=status )
  if (status /= 0) go to 990


  ! -----------------------------
  ! Start setting up the problem.
  ! -----------------------------
  if (unconstrained) then
     call CUTEST_usetup &
          ( status, input, nout, io_buffer, &
          n, x, bl(1:n), bu(1:n) )
  else
     call CUTEST_csetup &
          ( status, input, nout, io_buffer, n, m, &
          x, bl, bu, y, bl(n+1), bu(n+1), equatn, &
          linear, 0, 2, 1 )
  end if
  close (input)
  if (status /= 0) go to 910

  ! Compute number of nonlinear variables, and linear/equality constraints.
  if (unconstrained) then
     nnObj = n
     nnJac = 0
     mLCon = 0
     mNCon = 0
  else
     call CUTEST_cstats ( status, nnObj, nnJac, neq, mLCon )
     if (status /= 0) go to 910

     ! Compute the number of nonlinear constraints.
     mNCon = m - mLCon
  end if

  ! Compute the objective and constraints at x = 0.
  y(1:n) = zero
  if (unconstrained) then
     call CUTEST_ufn ( status, n, y, Obj )
  else
     call CUTEST_cfn ( status, n, m, y, Obj, c )
     if (status /= 0) go to 910
  end if

  ! Fix the bounds on the equality constraints
  ! Set the bounds on the linear   constraints.
  do j = 1, m
     if (equatn(j)) then
        bl(n+j) = zero
        bu(n+j) = zero
     end if
  end do

  do j = 1, m
     jslack = n + j
     if (j > mNCon) then
        bu(jslack) = bu(jslack) - c(j)
        bl(jslack) = bl(jslack) - c(j)
     end if

     ! ! If possible, set slack variables to be nonbasic at zero.
     ! x(jslack) = max( zero,       bl(jslack) )
     ! x(jslack) = min( x( jslack), bu(jslack) )
  end do


  ! Construct the linear constraint matrix.
  if ( .not. unconstrained ) then
     call cutest_ccfg ( status, n, m, x, c, .false., ldcj, n, cJac, .true. )
     if (status /= 0) go to 910

     do j = 1, n
        do i = 1, mNCon
           Jcon(i,j) = cJac(i,j)
        end do

        do i = 1, mLCon
           ic     = mNCon + i
           A(i,j) = cJac(ic,j)
        end do
     end do
     if (status /= 0) go to 910
  end if

  ! Set linear objective row (if any).
  iObjA = 0
  if (nnObj < n) then
     call cutest_cofg ( status, n, x, fObj, gObj, .true. )
     if (status /= 0) go to 910

     normG = zero
     do j = nnObj+1, n
        normG = normG + abs(gObj(j))
     end do

     if (normG > zero) then
        mLCon   = mLCon + 1
        m       = m     + 1
        bl(n+m) = -big
        bu(n+m) =  big
        iObjA   =  mLCon

        do j   = 1, n
           A(mLCon,j) = zero
        end do
        do j   = nnObj+1, n
           A(mLCon,j) = gObj(j)
        end do
     end if
  end if

  ! Set names, initialize some vectors.
  if (unconstrained) then
     call CUTEST_unames( status, n, pName, vname )
  else
     call CUTEST_cnames( status, n, m, pName, vname, cname )
  end if
  if (status /= 0) go to 910

  Names(1:n) = vname(1:n)
  if (unconstrained) then
     if (iObjA > 0) then
        Names(n+m)     = 'objectiv'
     end if
  else
     if (iObjA == 0) then
        Names(n+1:n+m) = cname(1:m)
     else
        Names(n+1:n+m-1) = cname(1:m-1)
        Names(n+m)       = 'objectiv'
     end if
  end if
  nNames = n + m


  state(1:nb) = 0
  x0(1:n)     = x(1:n)
  y(1:nb)     = zero
  objAdd      = zero
  if (nnObj == 0) objAdd = objAdd - Obj

  ! -------------------------------------
  ! Allocate initial workspace for DNOPT.
  ! -------------------------------------
  lencw  = 500
  leniw  = 500
  lenrw  = 500

  allocate ( iw(leniw), rw(lenrw), cw(lencw), STAT=status )


  ! ------------------------------
  ! Open DNOPT input/output files.
  ! Initialize DNOPT
  ! Read a Specs file.
  ! Solve the problem.
  ! ------------------------------
  open ( iSpecs, file = 'DNOPT.SPC', form = 'formatted', status = 'old' )
  open ( iPrint, file = trim(pName)//'.out', status = 'unknown' )


  call dnBegin( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

  call dnSpec( iSpecs, inform, cw, lencw, iw, leniw, rw, lenrw )
  if (inform .ne. 101 .and. inform .ne. 107) then
     if (inform .eq. 131 .and. nout .gt. 0) write( nout, 2010 )
     stop
  end if

  ! ------------------------------------------------------------------
  ! Compute an estimate of the storage needed by DNOPT.
  ! Copy the first 500 elements of cw, iw and rw into
  ! cw, iw and rw.  The required values are mincw, miniw and minrw.
  ! The default upper limits on the DNOPT workspace must be updated
  ! with these values.
  ! ------------------------------------------------------------------
  call dnMem                                   &
       ( inform,                               &
         mLCon, mNCon, n, nnJac, nnObj, iObjA, &
         mincw, miniw, minrw,                  &
         cw, lencw, iw, leniw, rw, lenrw )

  ! Allocate the full arrays  cw(lencw), iw(leniw), rw(lenrw).
  ! Copy the first 500 elements of cwtmp, iwtmp and rwtmp into
  ! cw, iw and rw.
  lencw = mincw
  leniw = miniw
  lenrw = minrw

  allocate( cw0(leniw) )
  allocate( iw0(leniw) )
  allocate( rw0(lenrw) )

  cw0(1:500) = cw(1:500)
  iw0(1:500) = iw(1:500)
  rw0(1:500) = rw(1:500)

  call move_alloc( from=cw0, to=cw )
  call move_alloc( from=iw0, to=iw )
  call move_alloc( from=rw0, to=rw )

  Errors =   0
  iP     =   0
  iS     =   0

  secondDeriv = iw(lvlHes) == 2

  call dnSetInt &
       ( 'Total character workspace', lencw, iP, iS, Errors, &
       cw, lencw, iw, leniw, rw, lenrw )
  call dnSetInt &
       ( 'Total integer   workspace', leniw, iP, iS, Errors, &
       cw, lencw, iw, leniw, rw, lenrw )
  call dnSetInt &
       ( 'Total real      workspace', lenrw, iP, iS, Errors, &
       cw, lencw, iw, leniw, rw, lenrw )

  ! ------------------------------------------------------------------
  ! Solve the problem, using a Cold start.
  ! ------------------------------------------------------------------
  start = 0              ! cold start

  if (secondDeriv) then
     call dnoptH                                     &
          ( start, n, mLCon, mNCon, nnJac, nnObj,    &
            pName, Names, nNames, iObjA, objAdd,     &
            DNOPT_evalcj, DNOPT_evalfg, DNOPT_evalh, &
            state, A, ldA, bl, bu,                   &
            fObj, gObj, fCon, JCon, ldJCon, H, ldH,  &
            objNP, nInf, sInf, x, y,                 &
            inform, mincw, miniw, minrw,             &
            cw, lencw, iw, leniw, rw, lenrw,         &
            cw, lencw, iw, leniw, rw, lenrw )
  else
     call dnopt                                      &
          ( start, n, mLCon, mNCon, nnJac, nnObj,    &
            pName, Names, nNames, iObjA, objAdd,     &
            DNOPT_evalcj, DNOPT_evalfg,              &
            state, A, ldA, bl, bu,                   &
            fObj, gObj, fCon, JCon, ldJCon, H, ldH,  &
            objNP, nInf, sInf, x, y,                 &
            inform, mincw, miniw, minrw,             &
            cw, lencw, iw, leniw, rw, lenrw,         &
            cw, lencw, iw, leniw, rw, lenrw )
  end if

!!$      if (inform .eq. 13) then
!!$         do j = 1, n
!!$!           x(j) = x(j) + abs(x(j))*1.0d-1
!!$            x(j) = x(j) + 0.5d+0*(x0(j) - x(j))
!!$         end do
!!$
!!$         start = 0              ! cold start
!!$         call dnopt                                  &
!!$           ( start, n, mLCon, mNCon, nnJac, nnObj,   &
!!$             pName, Names, nNames, iObjA, objAdd,    &
!!$             DNOPT_evalcj, DNOPT_evalfg,             &
!!$             state, A, ldA, bl, bu,                  &
!!$             fObj, gObj, fCon, JCon, ldJCon, H, ldH, &
!!$             objNP, nInf, sInf, x, y,                &
!!$             inform, mincw, miniw, minrw,            &
!!$             cw, lencw, iw, leniw, rw, lenrw,        &
!!$             cw, lencw, iw, leniw, rw, lenrw )
!!$      end if

  call CUTEST_creport( status, calls, cpu )
  if (status /= 0) go to 910

  write(nout,2000) pName, n, m, calls(1), calls(2), calls(5), &
                   calls(6), inform, fObj, cpu(1), cpu(2)

  call statsNP                                 &
       ( m, mNCon, n,                          &
         mNCon, nnObj, nnJac,                  &
         pName,                                &
         inform,                               &
         nInf, sInf, iObjA, objAdd, objNP, x,  &
         cw, lencw, iw, leniw, rw, lenrw )

  ! Try to handle abnormal DNOPT inform codes gracefully.

  if (inform .ge. 50 .and. &
       ( iPrint .gt. 0 .or. iSumm .gt. 0)) then

     if (iPrint .gt. 0) write ( iPrint, 3000) inform
     if (iSumm  .gt. 0) write ( iSumm , 3000) inform
  end if

  call dnEnd ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )


910 continue
  if (allocated( iu     )) deallocate ( iu     )
  if (allocated( ru     )) deallocate ( ru     )
  if (allocated( cu     )) deallocate ( cu     )

  if (allocated( cw     )) deallocate ( cw     )
  if (allocated( iw     )) deallocate ( iw     )
  if (allocated( rw     )) deallocate ( rw     )

  if (allocated( state  )) deallocate ( state  )
  if (allocated( names  )) deallocate ( names  )
  if (allocated( x      )) deallocate ( x      )
  if (allocated( y      )) deallocate ( y      )
  if (allocated( bl     )) deallocate ( bl     )
  if (allocated( bu     )) deallocate ( bu     )
  if (allocated( gObj   )) deallocate ( gObj   )
  if (allocated( fCon   )) deallocate ( fCon   )

  if (allocated( A      )) deallocate ( A      )
  if (allocated( Jcon   )) deallocate ( Jcon   )
  if (allocated( cJac   )) deallocate ( cJac   )
  if (allocated( H      )) deallocate ( H      )

  if (allocated( c      )) deallocate ( c      )
  if (allocated( equatn )) deallocate ( equatn )
  if (allocated( linear )) deallocate ( linear )
  if (allocated( vname  )) deallocate ( vname  )
  if (allocated( cname  )) deallocate ( cname  )

  close( iSpecs )
  close( iPrint )

  if (status /= 0) go to 920

  call CUTEST_cterminate( status )

  stop

920 continue
  write(nout, "( ' CUTEst error, status = ', i0, ', stopping' )") status
  stop

990 continue
  write(nout, "( ' Allocation error, status = ', I0 )" ) status
  stop

  ! Non-executable statements.
2010 format( /, ' ** PROGRAM DNOPT_main: No Specs file found.' )
3000 format( /, ' WARNING!  Abnormal DNOPT termination code:', &
                ' INFORM = ', I2 )
2000 format( /, 24('*'), ' CUTEst statistics ', 24('*') //             &
        ,' Code used               :  DNOPT',    /                     &
        ,' Problem                 :  ', A10,    /                     &
        ,' # variables             =      ', I10 /                     &
        ,' # constraints           =      ', I10 /                     &
        ,' # objective functions   =        ', F8.2 /                  &
        ,' # objective gradients   =        ', F8.2 /                  &
        ,' # constraints functions =        ', F8.2 /                  &
        ,' # constraints gradients =        ', F8.2 /                  &
        ' Exit code               =      ', I10 /                     &
        ,' Final f                 = ', E15.7 /                        &
        ,' Set up time             =      ', 0P, F10.2, ' seconds' /   &
        ' Solve time              =      ', 0P, F10.2, ' seconds' //  &
        66('*') / )

end program DNOPT_main

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine DNOPT_evalfg ( mode, nnObj, x, fObj, gObj, nState, &
                          cu, lencu, iu, leniu, ru, lenru )

  use DNOPT_mod, only : unconstrained

  implicit none
  integer          :: mode, nnObj, nState, lencu, leniu, lenru
  integer          :: iu(leniu)
  double precision :: fObj
  double precision :: x(nnObj), gObj(nnObj), ru(lenru)
  character*8      :: cu(lencu)

  logical          :: needG
  integer          :: status

  if (mode .eq. 0) then
     needG = .false.
  else
     needG = .true.
  end if

  if (unconstrained) then
     call cutest_uofg ( status, nnObj, x, fObj, gObj, needG )
  else
     call cutest_cofg ( status, nnObj, x, fObj, gObj, needG )
  end if

  if (status .ne. 0) THEN
     write( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") status
     stop
  end if

end subroutine DNOPT_evalfg

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine DNOPT_evalcj ( mode, mNCon, nnJac, x, fCon, &
                          Jcon, ldJ, nState, &
                          cu, lencu, iu, leniu, ru, lenru )

  implicit none
  integer          :: mode, mNCon, nnJac, ldJ, nState
  integer          :: lencu, leniu, lenru
  integer          :: iu(leniu)
  double precision :: x(nnJac), fCon(mNCon), Jcon(ldJ,*), ru(lenru)
  character*8      :: cu(lencu)

  logical          :: grad
  integer          :: status

  if ( mNCon > 0 ) then
     if (mode == 0) then
        grad = .false.
     else
        grad = .true.
     end if

     call cutest_ccfg ( status, nnJac, mNCon, x, fCon, .false., &
                        ldJ, nnJac, Jcon, grad )
  end if

end subroutine DNOPT_evalcj

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine DNOPT_evalH ( mode, nnH, mNCon,     &
                         x, y, H, ldH, nState, &
                         cu, lencu, iu, leniu, ru, lenru )
  use DNOPT_mod, only : unconstrained

  implicit none
  integer          :: mode, nnH, mNCon, ldH, nState
  integer          :: lencu, leniu, lenru
  integer          :: iu(leniu)
  double precision :: x(nnH), y(mNCon), H(ldH,*), ru(lenru)
  double precision, parameter   :: one = 1.0d+0
  character*8      :: cu(lencu)

  integer          :: status

  if (unconstrained) then
     call cutest_udh( status, nnH, x, ldH, H )
  else
     if (mNCon == 0) then
        call cutest_cidh( status, nnH, x, 0, ldH, H )
     else
        call dscal     ( mNCon, (-one), y, 1 )
        call cutest_cdh( status, nnH, mNCon, x, y, ldH, H )
        call dscal     ( mNCon, (-one), y, 1 )
     end if

  end if

end subroutine DNOPT_evalH

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine DNOPT_hv ( mode, nnH, mNCon, H, ldH, x, y, v, Hv, nState, &
                      cu, lencu, iu, leniu, ru, lenru )
  use DNOPT_mod, only : unconstrained

  implicit none
  integer          :: mode, nnH, mNCon, ldH, nState
  integer          :: lencu, leniu, lenru
  integer          :: iu(leniu)
  double precision :: x(nnH), y(mNCon), v(nnH), Hv(nnH), H(ldH,*), &
                      ru(lenru)
  double precision, parameter   :: one = 1.0d+0
  character*8      :: cu(lencu)

  integer          :: status
  logical          :: gotH

  if ( mode == 0 ) then
     ! Compute Hessian at x,y
     gotH = .false.
  else
     gotH = .true.
  end if

  ! Compute Hv = H*v
  if ( unconstrained ) then
     call cutest_uhprod( status, nnH, gotH, x, v, Hv )
  else
     call dscal        ( mNCon, (-one), y, 1 )
     call cutest_chprod( status, nnH, mNCon, gotH, x, y, v, Hv )
     call dscal        ( mNCon, (-one), y, 1 )
  end if

end subroutine DNOPT_hv

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
