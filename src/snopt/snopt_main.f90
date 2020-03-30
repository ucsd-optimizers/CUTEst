
program snopt_main

  use            CUTEst_interface_double
  implicit       none
  integer,       parameter   :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  character                  :: Start*8, pName*10, filename*20

  integer(ip)                :: leniu, lenru, leniw, lenrw, lencw,           &
                                mincw, miniw, minrw
  integer(ip)                :: iObj, n, m, m0, m1, nm, neH, neJ,            &
                                neG, nlocG, nNames, nnCon, nnH, nnJac, nnObj,&
                                nLC, nEq, nS, nInf
  integer(ip)                :: imax1, imax2, imax3, rmax1, rmax2, rmax3
  integer(ip)                :: iuStart, ruStart, j, jslack
  integer(ip)                :: liGfun, lenH, lenJ, ljGvar, lx, ly
  integer(ip)                :: lcolH, lrowH, llocG, llocH, lfCon,           &
                                lG, lgObj, lgCon, lvalH, nlocH
  integer(ip)                :: Errors, INFO, status, tries, alloc_stat

  real(rp)                   :: fobj, objAdd, sInf
  real(rp)                   :: CPU( 2 ), CALLS( 7 )

  logical                    :: Constrained, UseHessian

  integer(ip)                :: iPrint, iSumm

  character,     allocatable :: cNames(:)*12
  character,     allocatable :: Names(:)*8
  integer(ip),   allocatable :: hs(:), indH(:), indJ(:), locH(:), locJ(:)
  real(rp),      allocatable :: bl(:), bu(:), x(:), x0(:), pi(:), rc(:)
  real(rp),      allocatable :: valJ(:), c(:)
  real(rp),      allocatable :: valH(:)
  logical,       allocatable :: equation(:), linear(:)

  integer(ip),   allocatable :: iw(:), iu(:), iw0(:)
  real(rp),      allocatable :: rw(:), ru(:), rw0(:)
  character(8),  allocatable :: cw(:), cw0(:)

  integer(ip),   parameter   :: iCutest  = 55,  iOut   = 6, io_buffer = 11
  integer(ip),   parameter   :: iPP      = 14
  integer(ip),   parameter   :: maxtries = 10
  real(rp),      parameter   :: zero     = 0.0d+0, infBnd = 1.0d+20

  external                   :: SNOPT_evalCon, SNOPT_evalObj, SNOPT_userfun
  external                   :: SNOPTCH_userfun

  alloc_stat = 0

  !-----------------------
  ! Open problem data file
  !-----------------------
  open  ( iCutest, file = 'OUTSDIF.d', form = 'FORMATTED', status = 'old' )
  rewind( iCutest )

  !----------------------------------
  ! Get dimensions and allocate space
  !----------------------------------
  call CUTEST_cdimen( status, iCutest, n, m )
  if (status /= 0) go to 910

  UseHessian  = .false.
  Constrained = m > 0
  m0 = m
  m1 = m + 1

  !-------------------------
  ! Start setting up problem
  !-------------------------
  nm = n + m1

  allocate( hs(nm), bl(nm), bu(nm), x(nm), x0(nm), pi(m1), rc(nm), &
            Names(nm), cNames(nm),                                 &
            c(m1), equation(m), linear(m), stat=alloc_stat )
  if (alloc_stat /= 0) GO TO 900

  ! Set bounds and initial x

  if (Constrained) then ! The problem has general constraints.
     call CUTEST_csetup( status, iCutest, iOut, io_buffer, &
                         n, m, x(1:n), bl(1:n), bu(1:n),   &
                         pi, bl(n+1:nm), bu(n+1:nm),       &
                         equation, linear, 0, 2, 1 )
  else            ! The problem is unconstrained.
     call CUTEST_usetup( status, iCutest, iOut, io_buffer, &
                         n, x, bl(1:n), bu(1:n) )
  end if
  if (status /= 0) go to 910

  ! Fix the bounds on the equality constraints
  do j = 1, m
     if (equation(j)) then
        bl(n+j) = zero
        bu(n+j) = zero
     end if
  end do

  ! Compute number of nonlinear variables, and linear/equality constraints.
  if (Constrained) then
     call CUTEST_cstats( status, nnObj, nnJac, nEq, nLC )
  else
     nnObj = n
     nnJac = 0
     nEq   = 0
     nLC   = 0
  end if
  if (status /= 0) go to 910

  ! Compute the objective and constraint functions at x = 0.
  rc(:) = zero
  if (Constrained) then
     call CUTEST_cfn( status, n, m, rc, fobj, c )
  else
     call CUTEST_ufn( status, n, rc, fobj )
  end if
  if (status /= 0) go to 910

  ! Compute the number of nonlinear constraints.

  nnCon = m - nLC
  nnH   = max( nnObj, nnCon )
  nlocG = nnJac + 1
  nlocH = nnH   + 1

  ! Set names
  if (Constrained) then
     call CUTEST_cnames( status, n, m, pName, cNames(1:n), cNames(n+1:n+m) )
  else
     call CUTEST_unames( status, n, pName, cNames(1:n) )
  end if
  if (status /= 0) go to 910

  iuStart = 7 ! address of available integer workspace.
  ruStart = 1 ! address of available real    workspace.

  !------------------------------------------------------------------------
  ! Compute all dimensions.
  !------------------------------------------------------------------------

  if (Constrained) then
     ! =======================================
     ! The problem has general constraints.
     ! neJ  includes n obj gradient components.
     ! ========================================
     call cutest_cdimsj( status, neJ ) ! includes n obj gradient components
     if (status /= 0) go to 910

     if (nnH > 0) then
        call cutest_cdimsh( status, neH )
        if (alloc_stat /= 0) go to 900
     else
        neH = 0
     end if
  else
     ! ===================================
     ! The problem is unconstrained.
     ! ===================================
     neJ  = 1 ! dummy free row added here.

     if (nnH > 0) then
        call cutest_udimsh( status, neH )
        if (alloc_stat /= 0) go to 900
     else
        neH = 0
     end if
  end if

  ! Integer and real workspace for funobj

  lx     = ruStart           !    x(1:n)
  lgObj  = lx      + n       ! gObj(1:n)
  rmax1  = lgObj   + n - 1

  imax1  = iuStart - 1

  ! Integer and real workspace for funCon

  lx     = ruStart           !     x(1:n)
  lfcon  = lx      + n       !  fCon(1:m)
  lgCon  = lfCon   + m1      !  gCon(1:neJ)
  lG     = lgCon   + neJ     !     G(1:neJ)
  rmax2  = lG      + neJ - 1

  llocG  = iuStart           !  locG(1:nlocG)
  liGfun = llocG   + nlocG   ! iGfun(1:neJ)
  ljGvar = liGfun  + neJ     ! iGvar(1:neJ)
  imax2  = ljGvar  + neJ - 1

  ! Integer and real workspace for funHes

  lx     = ruStart           !     x(1:n)
  ly     = lx      + n       !     y(1:m)
  lvalH  = ly      + m1      !  valH(1:neH)
  rmax3  = lvalH   + neH - 1

  llocH  = iustart           !  locH(1:nlocH)
  lrowH  = llocH   + nlocH   !  rowH(1:neH)
  lcolH  = lrowH   + neH     !  colH(1:neH)
  imax3  = lcolH   + neH - 1

  if (UseHessian) then
     lenru  = max(rmax1,rmax2,rmax3)
     leniu  = max(imax1,imax2,imax3)
  else
     lenru  = max(rmax1,rmax2)
     leniu  = max(imax1,imax2)
  end if

  !------------------------------------------------------------------------
  ! All dimensions are known (or over-estimated)
  ! Allocate workspace.
  !------------------------------------------------------------------------

  allocate( valJ(neJ), locJ(n+1), indJ(neJ), stat=alloc_stat )
  if (alloc_stat /= 0) go to 900

  if (UseHessian .and. nnH > 0) then
     allocate( valH(neH), locH(nnH+1), indH(neH), stat=alloc_stat )
     if (alloc_stat /= 0) go to 900
  end if

  allocate( iu(leniu), ru(lenru), stat=alloc_stat )
  if (alloc_stat /= 0) GO TO 900

  !------------------------------------------------------------------------
  ! Build the Jacobian and Hessian.
  !------------------------------------------------------------------------

  if (Constrained) then
     ! ===================================
     ! The problem has general constraints.
     ! ===================================
     Names(1:n+m) = cNames(1:n+m)(1:8)

     lenJ  = neJ
     call buildJac( status, n, m1, nnObj, nnCon, nnJac, lenJ,     &
                    iObj, locJ, indJ, valJ, nlocG, iu(llocG),     &
                    iu(liGfun), iu(ljGvar), ru, x, pi, neJ, neG )
     if (status /= 0) go to 910

     ! Set the bounds on the linear objective row (if there is one).
     if (iObj > 0) then
        m       =  m1
        bl(n+m) = -infBnd
        bu(n+m) =  InfBnd
        c(m)    =  zero

        Names(n+m) = 'OBJ     '
     end if

     nNames = n + m
     nNames = 1

     ! Set bounds on linear constraints.
     do j      = 1, m
        jslack = n + j
        if (j > nnCon) then
           bu(jslack) = bu(jslack) - c(j)
           bl(jslack) = bl(jslack) - c(j)
        end if

        ! If possible, set slack variables to be nonbasic at zero.
        x(jslack) = max( zero, bl(jslack) )
        x(jslack) = min( x( jslack), bu(jslack) )
     end do
  else
     ! ====================================================================
     ! The problem is unconstrained.
     ! Add a dummy free row.
     ! ====================================================================
     ! This part needs changing so that any linear objective
     ! gradients are held in Jcol. (Unlikely scenario though)

     iObj =  0
     neG  =  0

     ! Add a dummy linear objective row

     m           = m1
     valJ(neJ)   = zero
     indJ(neJ)   = 1
     locJ(1)     = 1
     locJ(2:n+1) = neJ + 1
     locJ(n+1)   = neJ + 1

     bl(n+m)     = -infBnd
     bu(n+m)     =  InfBnd
     c(m)        =  zero

     nNames      = n + m
     Names(n+m)  = pName(1:8)
     nNames      = 1
  end if

  ! Build the Hessian

  if (UseHessian .and. nnH > 0) then
     lenH  = neH
     rc(1:n+nnCon) = zero
     call buildHess( status, n, m, nnH, lenH,                  &
                     nlocH, locH, indH, valH,                  &
                     iu(lrowH), iu(lcolH), ru, x, rc(n+1), neH  )
  end if

  ! All the workspace has been allocated and assigned.

  iu(1)  = iuStart ! address of available integer workspace.
  iu(2)  = ruStart ! address of available real    workspace.
  iu(3)  = m0
  iu(4)  = n
  iu(5)  = neJ
  iu(6)  = neH

  hs(:)  = 0
  pi(:)  = zero
  objAdd = zero

  if (nnObj == 0) then
     objAdd = objAdd - fobj
  end if

  Start = 'cold'
  x0(:) = x(1:nm)

  !-----------------------------------------------------------------------------
  ! Ok, we're done with CUTEst stuff.
  !-----------------------------------------------------------------------------

  ! Allocate workspace
  leniw = 20000000
  lenrw = 40000000
  lencw = 500
  allocate( cw(lencw), iw(leniw), rw(lenrw), stat=alloc_stat )
  if (alloc_stat /= 0) go to 900

  filename = trim(pName)//'.out'
  call snInitF &
       (filename, 'screen', iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw)

  ! Read specs file
  call snSpecF('SNOPT.SPC', INFO, cw, lencw, iw, leniw, rw, lenrw)
  if (INFO /= 101  .and.  INFO /= 107) go to 920


  ! Check workspace
  call snMemB &
       (INFO, m, n, neJ, neG, nnCon, nnJac, nnObj, &
        mincw, miniw, minrw, cw, lencw, iw, leniw, rw, lenrw)

  do tries = 1, maxtries

     if (leniw .le. miniw) then
        ! Not enough integer space
        if (tries == 1) then
           leniw = miniw
        else
           leniw = 10*miniw
        end if
        allocate( iw0(leniw), stat=alloc_stat )
        if (alloc_stat /= 0) GO TO 900

        iw0(1:500) = iw(1:500)

        call move_alloc( from=iw0, to=iw )

        call snSeti &
           ( 'Total integer   workspace', leniw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw )
     end if

     if (lenrw .le. minrw) then
        ! Not enough real space
        if (tries == 1) then
           lenrw = minrw
        else
           lenrw = 10*minrw
        end if
        allocate( rw0(lenrw), stat=alloc_stat )
        if (alloc_stat /= 0) GO TO 900

        rw0(1:500) = rw(1:500)

        call move_alloc( from=rw0, to=rw )

        call snSeti &
           ( 'Total real      workspace', lenrw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw )
     end if

     if (lencw .le. mincw) then
        ! Not enough character space
        if (tries == 1) then
           lencw = mincw
        else
           lencw = 10*mincw
        end if
        allocate( cw0(lencw), stat=alloc_stat )
        if (alloc_stat /= 0) GO TO 900

        cw0(1:500) = cw(1:500)

        call move_alloc( from=cw0, to=cw )

        call snSeti &
           ( 'Total character workspace', lencw, 0, 0, Errors, &
             cw, lencw, iw, leniw, rw, lenrw )
     end if

     ! Solve the problem

     if (UseHessian) then
        call snOptCH( start, m, n  , neJ  , neH   , nNames,         &
                      nnCon, nnObj , nnJac, nnH   ,                 &
                      iObj , objAdd, pName(1:8),                    &
                      SNOPTCH_userfun,                              &
                      valJ , indJ  , locJ ,                         &
                      valH , indH  , locH ,                         &
                      bl   , bu    , Names,                         &
                      hs   , x     , pi   , rc    ,                 &
                      INFO , mincw , miniw, minrw ,                 &
                      nS   , nInf  , sInf , fObj  ,                 &
                      cw   , lencw , iu   , leniu , ru    , lenru,  &
                      cw   , lencw , iw   , leniw , rw    , lenrw )
     else
!!$        call snoptC ( start, m , n , neJ  , nNames,                 &
!!$                      nnCon, nnObj , nnJac,                         &
!!$                      iObj , objAdd, pName(1:8),                    &
!!$                      SNOPT_userfun,                                &
!!$                      valJ , indJ  , locJ , bl    , bu    , names , &
!!$                      hs   , x     , pi   , rc    ,                 &
!!$                      INFO , mincw , miniw, minrw ,                 &
!!$                      nS   , nInf  , sInf , fobj  ,                 &
!!$                      cw   , lencw , iu   , leniu , ru    , lenru , &
!!$                      cw   , lencw , iw   , leniw , rw    , lenrw )

        call snoptB ( start, m , n , neJ  , nNames,                 &
                      nnCon, nnObj , nnJac,                         &
                      iObj , objAdd, pName(1:8),                    &
                      SNOPT_evalCon, SNOPT_evalObj,                 &
                      valJ , indJ  , locJ , bl    , bu    , names , &
                      hs   , x     , pi   , rc    ,                 &
                      INFO , mincw , miniw, minrw ,                 &
                      nS   , nInf  , sInf , fobj  ,                 &
                      cw   , lencw , iu   , leniu , ru    , lenru , &
                      cw   , lencw , iw   , leniw , rw    , lenrw )
     endif

     if (INFO/10 /= 8) then
        exit
     end if

     ! Ran out of memory. Try to allocate some more.
     ! Start from scratch.
     x(:)  = x0(1:nm)
     hs(:) = 0
     pi(:) = zero
  end do

  call s4npGetStats( m, n,                                  &
                     nnCon, nnObj, nnJac,                   &
                     pName,                                 &
                     INFO,                                  &
                     nS, nInf, sInf, iObj, objAdd, fobj, x, &
                     cw, lencw, iw, leniw, rw, lenrw )

  CALLS = 0.0
  if (Constrained) then
     call CUTEST_creport( status, CALLS, CPU )
  else
     call CUTEST_ureport( status, CALLS, CPU )
  end if
  write ( iOut, 2000 ) pName, n, m, CALLS( 1 ), CALLS( 2 ),                     &
                       CALLS( 5 ), CALLS( 6 ), INFO, fobj, CPU( 1 ), CPU( 2 )

900 continue
  if (allocated( Names  )  ) deallocate( Names  )
  if (allocated( cNames )  ) deallocate( cNames )

  if (allocated( bl )      ) deallocate( bl  )
  if (allocated( bu )      ) deallocate( bu  )
  if (allocated( hs )      ) deallocate( hs  )
  if (allocated( x  )      ) deallocate( x   )
  if (allocated( x0 )      ) deallocate( x0  )
  if (allocated( pi )      ) deallocate( pi  )
  if (allocated( rc )      ) deallocate( rc  )

  if (allocated( c  )      ) deallocate( c   )

  if (allocated( indJ )    ) deallocate( indJ )
  if (allocated( locJ )    ) deallocate( locJ )
  if (allocated( valJ )    ) deallocate( valJ )

  if (allocated( indH )    ) deallocate( indH )
  if (allocated( locH )    ) deallocate( locH )
  if (allocated( valH )    ) deallocate( valH )

  if (allocated( iw )      ) deallocate( iw  )
  if (allocated( rw )      ) deallocate( rw  )
  if (allocated( cw )      ) deallocate( cw  )

  if (allocated( iw0)      ) deallocate( iw0 )
  if (allocated( rw0)      ) deallocate( rw0 )
  if (allocated( cw0)      ) deallocate( cw0 )

  if (allocated( iu)       ) deallocate( iu  )
  if (allocated( ru)       ) deallocate( ru  )

  if (allocated( equation )) deallocate( equation )
  if (allocated( linear )  ) deallocate( linear   )

  close (iCutest)

  if (Constrained) then
     call CUTEST_cterminate( status )
  else
     call CUTEST_uterminate( status )
  end if

  if (alloc_stat /= 0) GO TO 990

  stop

910 continue
  write( iOut, "( ' CUTEst error, status = ', i0, ', stopping' )") status
  stop

920 continue
  write( iOut, "( ' Specs file error, INFO= ', i0 )" ) INFO
  stop

990 continue
  write( iOut, "( ' Allocation error, status = ', i0 )" ) status
  stop

2000 format( /, 24('*'), ' CUTEst statistics ', 24('*') //                    &
         ,' Package used            : SNOPT',    /                            &
         ,' Problem                 :  ', A10,    /                           &
         ,' # variables             =      ', I10 /                           &
         ,' # constraints           =      ', I10 /                           &
         ,' # objective functions   =        ', F8.2 /                        &
         ,' # objective gradients   =        ', F8.2 /                        &
         ,' # constraints functions =        ', F8.2 /                        &
         ,' # constraints gradients =        ', F8.2 /                        &
         ,' Exit code               =      ', I10 /                           &
         ,' Final f                 = ', E15.7 /                              &
         ,' Set up time             =      ', 0P, F10.2, ' seconds' /         &
         ,' Solve time              =      ', 0P, F10.2, ' seconds' //        &
         66('*') / )

end program snopt_main

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine buildJac( status, n, m1, nnObj, nnCon, nnJac, lenJ, &
                     iObj, locJ, indJ, Jcol, nlocG,            &
                     locG, iGfun, jGvar, G, x, y, neJ, neG )

  use          CUTEst_interface_double
  implicit     none
  integer,     parameter   :: ip = kind( 1 ), rp = kind( 1.0d+0 )

  integer(ip), intent(in)  :: n, m1, nnObj, nnCon, nnJac, lenJ, nlocG
  real(rp),    intent(in)  :: x(n), y(m1)
  integer(ip), intent(out) :: status, iObj, neJ, neG, locJ(n+1), indJ(lenJ), &
                              locG(nlocG), iGfun(lenJ), jGvar(lenJ)
  real(rp),    intent(out) :: Jcol(lenJ), G(lenJ)

  !===========================================================================
  ! On entry, nlocG = nnJac + 1
  !===========================================================================
  real(rp)                 :: normObj
  integer(ip)              :: i, j, k, l, neJac
  real(rp),    parameter   :: zero     = 0.0d+0

  ! Compute  neJ, G, jGvar, iGfun and status

  call cutest_csgr( status, n, m1, x, y, .false., neJ, lenJ, &
                    G, jGvar, iGfun )
  if (status /= 0) return

  ! Initialize
  locJ(1:n)     = 0
  locG(1:nnJac) = 0

  ! Find the norm of the constant objective gradients.
  ! If the norm is nonzero, define a linear objective row.
  normObj = zero
  do k = 1, neJ
     j = jGvar(k)
     i = iGfun(k)
     if (i .eq. 0  .and.  j .gt. nnObj) then ! linear obj
        normObj = normObj + abs(G(k))
     end if
  end do

  ! normObj = 1.0d+0

  if (normObj > zero) then
     iObj = m1
  else
     iObj = 0
  end if

  ! Count the Jacobian entries in column j.
  ! Initialize locJ.  Store the column counts in locJ(1:n).
  ! Don't include nonlinear objective function entries.
  neJac = neJ
  do k  = 1, neJ
     j  = jGvar(k)
     i  = iGfun(k)
     if (i .eq. 0 .and. (j .le. nnObj .or. iObj .eq. 0)) then
        neJac   = neJac   - 1
     else
        locJ(j) = locJ(j) + 1
     end if
  end do
  locJ(n+1) = neJac + 1

  ! Set locJ(j) to point to start of column j.
  do j = n, 1, -1
     locJ(j) = locJ(j+1) - locJ(j)
  end do

  ! Load nonlinear Jacobian entries into Jcol and indJ.
  ! Use locJ to keep track of position for each variable b.
  ! Also count nonlinear Jacobian entries in each column and
  ! store count in locG(j).
  neG = 0
  do k = 1, neJ
     j = jGvar(k)
     i = iGfun(k)
     if (i .gt. 0 .and. i .le. nnCon .and. j .le. nnJac) then
        l       = locJ(j)
        Jcol(l) = G(k)
        indJ(l) = i
        locJ(j) = l   + 1
        locG(j) = locG(j) + 1
        neG     = neG + 1
     end if
  end do

  ! Load the constant Jacobian entries, including linear objective
  ! function entries,  in Jcol and indJ.
  ! locJ(j) points to the next available position in column j.
  do k = 1, neJ
     j = jGvar(k)
     i = iGfun(k)
     if (i .eq. 0 .and. j .gt. nnObj .and. iObj > 0) then ! linear obj
        l       = locJ(j)
        Jcol(l) = G(k)
        indJ(l) = iObj
        locJ(j) = l + 1
     else if (i .gt. nnCon .or. (i .gt. 0 .and. j .gt. nnJac)) then
        l       = locJ(j)
        Jcol(l) = G(k)
        indJ(l) = i
        locJ(j) = l + 1
     end if
  end do

  ! Reset locJ  and set locG.
  ! locG(j)  points  to the start of column j in Gcon,
  ! the nonlinear Jacobian.  These pointers are used by usrfgh.
  locG(nlocG) = neG + 1

  do  j = n, 2, -1
     locJ(j) = locJ(j-1)
  end do
  locJ(1) = 1

  do  j = nnJac, 2, -1
     locG(j) = locG(j+1) - locG(j)
  end do
  locG(1) = 1

  neJ = neJac

  if (neJ == 0) then
     neJ       = 1
     Jcol(neJ) = zero
     indJ(neJ) = 1
     locJ(n+1) = neJ + 1
  end if

end subroutine buildJac

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine buildHess( status, n, m, nnH, lenH,            &
                      nlocH, locH, indH, Hcol,            &
                      rowH, colH, valH, x, y, neH  )

  use           CUTEst_interface_double
  implicit      none
  integer,      parameter     :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  integer(ip),  intent(in)    :: lenH, m, n, nlocH, nnH
  real(rp),     intent(in)    :: x(nnH), y(m)

  integer(ip),  intent(out)   :: rowH(lenH), colH(lenH)
  real(rp),     intent(out)   :: valH(lenH)

  integer(ip),  intent(out)   :: status, neH, locH(nnH+1), indH(lenH)
  real(rp),     intent(out)   :: Hcol(lenH)

  !=============================================================================
  ! Build the Hessian of the Lagrangian.
  !
  !=============================================================================
  integer(ip)                 :: i, ii, j, k, nval
  real(rp),    parameter      :: zero = 0.0d+0

  ! Compute  neH, valH, rowH, colH and status

  if (m > 0) then
     call cutest_csh( status, n, m, x, y, neH, lenH, valH, rowH, colH )
  else
     call cutest_ush( status, n, x, neH, lenH, valH, rowH, colH )
  end if
  if (status /= 0) go to 910

  ! Convert coordinate form to sparse-by-col form
  nval = 0
  locH = 0
  do k = 1, neH
     i = rowH(k)
     j = colH(k)

     locH(j) = locH(j) + 1
     nval    = nval    + 1
  end do

  locH(nlocH) = nval + 1

  do k = nnH, 1, -1
     locH(k) = locH(k+1) - locH(k)
  end do

  indH = 0
  Hcol = zero
  do k = 1, neH
     j  = colH(k)
     i  = rowH(k)
     ii = locH(j)

     Hcol(ii) = valH(k)
     indH(ii) = i
     locH(j)  = ii + 1
  end do

  locH(2:nlocH) = locH(1:nlocH-1)
  locH(1)       = 1

  return

910 continue
  write( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") status
  stop

end subroutine buildHess

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine SNOPT_userfun( mode, nnObj, nnCon, nnJac, nnH,          &
                          neG, x,                                  &
                          fObj, gObj, fCon, gCon, statusUser,      &
                          cu, lencu, iu, leniu, ru, lenru )

  use           CUTEst_interface_double
  implicit      none
  integer,      parameter     :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  integer(ip),  intent(in)    :: mode, nnObj, nnCon, nnJac, nnH, neG, &
                                 statusUser, lencu, leniu, lenru
  real(rp),     intent(in)    :: x(nnH)

  character(8), intent(inout) :: cu(lencu)
  integer(ip),  intent(inout) :: iu(leniu)
  real(rp),     intent(inout) :: ru(lenru)

  real(rp),     intent(out)   :: fObj, gObj(nnObj), fCon(nnCon), gCon(neG)

  !=============================================================================
  ! For SNOPTCH
  !=============================================================================
  if (nnObj .gt. 0) then
     call SNOPT_evalObj( mode, nnObj,                        &
                         x, fObj, gObj, statusUser,          &
                         cu, lencu, iu, leniu, ru, lenru )
  end if

  if (nnCon .gt. 0) then
     call SNOPT_evalCon( mode, nnCon, nnJac, neG,            &
                         x, fCon, gCon, statusUser,          &
                         cu, lencu, iu, leniu, ru, lenru )
  end if

end subroutine SNOPT_userfun

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine SNOPTCH_userfun( mode, nnObj, nnCon, nnJac, nnH,          &
                            neG, neH, x, yCon,                       &
                            fObj, gObj, fCon, gCon, Hcol, statusUser,&
                            cu, lencu, iu, leniu, ru, lenru )

  use  CUTEst_interface_double
  implicit      none

  integer,      parameter     :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  integer(ip),  intent(in)    :: mode, nnObj, nnCon, nnJac, nnH, neG, neH, &
                                 statusUser, lencu, leniu, lenru
  real(rp),     intent(in)    :: x(nnH), yCon(nnCon)

  character(8), intent(inout) :: cu(lencu)
  integer(ip),  intent(inout) :: iu(leniu)
  real(rp),     intent(inout) :: ru(lenru)

  real(rp),     intent(out)   :: fObj, gObj(nnObj), fCon(nnCon), gCon(neG),&
                                 Hcol(neH)

  !=============================================================================
  ! For SNOPTCH
  !=============================================================================
  if (nnObj .gt. 0) then
     call SNOPT_evalObj( mode, nnObj,                        &
                         x, fObj, gObj, statusUser,          &
                         cu, lencu, iu, leniu, ru, lenru )
  end if

  if (nnCon .gt. 0) then
     call SNOPT_evalCon( mode, nnCon, nnJac, neG,            &
                         x, fCon, gCon, statusUser,          &
                         cu, lencu, iu, leniu, ru, lenru )
  end if

  if (nnH .gt. 0) then
     call SNOPT_evalHes( nnCon, nnH, neH,                    &
                         x, yCon, Hcol, statusUser,          &
                         cu, lencu, iu, leniu, ru, lenru )
  end if

end subroutine SNOPTCH_userfun

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine SNOPT_evalObj( mode, nnObj, x, fObj, gObj, statusUser, &
                          cu, lencu, iu, leniu, ru, lenru )

  use           CUTEst_interface_double
  implicit      none
  integer,      parameter     :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  integer(ip),  intent(in)    :: mode, nnObj, statusUser, lencu, leniu, lenru
  real(rp),     intent(in)    :: x(nnObj)

  character(8), intent(inout) :: cu(lencu)
  integer(ip),  intent(inout) :: iu(leniu)
  real(rp),     intent(inout) :: ru(lenru)

  real(rp),     intent(out)   :: fObj, gObj(nnObj)

  !=============================================================================
  ! On entry, nnObj > 0.
  !=============================================================================
  logical                     :: needG
  integer(ip)                 :: ruStart, lx, lgObj, n, m, status
  real(rp),      parameter    :: zero = 0.0d+0

  needG   = mode > 0

  ruStart = iu(2) ! address of available real    workspace.
  m       = iu(3)
  n       = iu(4)

  lx      = ruStart           !    x(1:n)
  lgObj   = lx      + n       ! gObj(1:n)

  ru(lx:lx+2*n-1)   = zero
  ru(lx:lx+nnObj-1) = x

  if (m > 0) then
     call cutest_cofg( status, n, ru(lx), fObj, ru(lgObj), needG )
  else
     call cutest_uofg( status, n, ru(lx), fObj, ru(lgObj), needG )
  end if

  if (needG) then
     gObj = ru(lgObj:lgObj+nnObj-1)
  end if

  if (status /= 0) then
     write( 6, "( ' CUTEst error, status = ', i0, ', stopping' )")  status
     stop
  end if

end subroutine SNOPT_evalObj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine SNOPT_evalCon( mode, nnCon, nnJac, neG, x, &
                          fCon, gCon, statusUser, &
                          cu, lencu, iu, leniu, ru, lenru )

  use  CUTEst_interface_double
  implicit none

  integer,      parameter     :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  integer(ip),  intent(in)    :: mode, nnCon, nnJac, neG, statusUser, &
                                 lencu, leniu, lenru
  real(rp),     intent(in)    :: x(nnJac)

  character(8), intent(inout) :: cu(lencu)
  integer(ip),  intent(inout) :: iu(leniu)
  real(rp),     intent(inout) :: ru(lenru)

  real(rp),     intent(out)   :: fCon(nnCon), gCon(neG)

  !=============================================================================
  ! On entry, nnCon > 0.
  !=============================================================================
  logical                     :: needG
  integer(ip)                 :: iuStart, ruStart, liGfun, ljGvar, llocG, lG, &
                                 m, nlocG, n, neJ, lfcon, lgcon, lx
  real(rp),     parameter     :: zero = 0.0d+0

  needG   = mode > 0

  iuStart = iu(1) ! address of available integer workspace.
  ruStart = iu(2) ! address of available real    workspace.
  m       = iu(3)
  n       = iu(4)
  neJ     = iu(5)

  nlocG   = nnJac + 1

  lx      = ruStart           !     x(1:n)
  lfcon   = lx      + n       !  fCon(1:m)
  lgCon   = lfCon   + m       !  gCon(1:neJ)
  lG      = lgCon   + neJ     !     G(1:neJ)

  llocG   = iuStart           !  locG(1:nlocG)
  liGfun  = llocG   + nlocG   ! iGfun(1:neJ)
  ljGvar  = liGfun  + neJ     ! iGvar(1:neJ)

  ru(lx:lx+n-1)     = zero
  ru(lx:lx+nnJac-1) = x

  call SNOPT_evalCon0( needG, m, n, nnCon, nnJac, neJ,                 &
                       ru(lx), ru(lfcon), ru(lgcon), statusUser,       &
                       nlocG, iu(llocG), iu(liGfun), iu(ljGvar), ru(lG) )

  fCon(1:nnCon) = ru(lfcon:lfcon+nnCon-1)
  if (needG) then
     gCon(1:neG) = ru(lgcon:lgcon+neG-1)
  end if

end subroutine SNOPT_evalCon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine SNOPT_evalCon0( needG, m, n, nnCon, nnJac, lenG, &
                           x, fCon, gCon, statusUser,       &
                           nlocG, locG, iGfun, jGvar, G )

  use           CUTEst_interface_double
  implicit      none
  integer,      parameter     :: ip = kind( 1 ), rp = kind( 1.0D+0 )
  logical,      intent(in)    :: needG
  integer(ip),  intent(in)    :: m, n, nnCon, nnJac, lenG, statusUser, nlocG
  real(rp),     intent(in)    :: x(n)

  integer(ip),  intent(inout) :: locG(nlocG)

  integer(ip),  intent(out)   :: iGfun(lenG), jGvar(lenG)
  real(rp),     intent(out)   :: fCon(m), gCon(lenG), G(lenG)

  !===========================================================================
  ! Does the work for SNOPT_evalCon.
  !
  !===========================================================================
  integer(ip)                 :: j, k, l, nnzJ, status

  ! Evaluate the problem constraints. On input, nnCon > 0.
  ! The Jacobian is stored in sparse format.
  call cutest_ccfsg( status, n, m, x, fCon, nnzJ, lenG, &
                     G, jGvar, iGfun, needG )
  if (status /= 0) go to 910

  if (needG) then
     ! Copy the Jacobian from CSCFG, stored in (G,iGfun,jGvar)
     ! into the SNOPT Jacobian stored by columns in gCon.
     ! locG(j) points to the next available position in column j.
     do l      = 1, nnzJ
        j      = jGvar(l)
        if (j <= nnJac .and. iGfun(l) <= nnCon) then
           k       = locG(j)
           gCon(k) = G(l)
           locG(j) = k + 1
        end if
     end do

     ! Reset locG.
     do j = nnJac, 2, -1
        locG(J) = locG(J-1)
     end do
     locG(1) = 1
  end if

  return

910 continue
  write( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") status
  stop

end subroutine SNOPT_evalCon0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine SNOPT_evalHes( nnCon, nnH, neH,                           &
                          x, yCon, Hcol, statusUser,                 &
                          cu, lencu, iu, leniu, ru, lenru )

  use           CUTEst_interface_double
  implicit      none
  integer,      parameter     :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  integer(ip),  intent(in)    :: nnCon, nnH, neH, statusUser,        &
                                 lencu, leniu, lenru
  real(rp),     intent(in)    :: x(nnH), yCon(nnCon)

  real(rp),     intent(out)   :: Hcol(neH)

  character(8), intent(inout) :: cu(lencu)
  integer(ip),  intent(inout) :: iu(leniu)
  real(rp),     intent(inout) :: ru(lenru)

  !=============================================================================
  ! Compute the Hessian at the given x.
  ! On entry, nnH > 0.
  !=============================================================================
  integer(ip)                 :: iuStart, ruStart, m, n, nlocH
  integer(ip)                 :: llocH, lcolH, lrowH, lvalH, lx, ly

  iuStart = iu(1) ! address of available integer workspace.
  ruStart = iu(2) ! address of available real    workspace.
  m       = iu(4)
  n       = iu(5)

  nlocH   = nnH     + 1

  lx      = ruStart          !     x(1:n)
  ly      = lx      + n      !     y(1:m)
  lvalH   = ly      + m      !  valH(1:neH)

  llocH   = iustart          !  locH(1:nlocH)
  lrowH   = llocH   + nlocH  !  rowH(1:neH)
  lcolH   = lrowH   + neH    !  colH(1:neH)

  ru(lx:lx+nnH-1) =  x

  if (m > 0) then
     ru(ly:ly+nnCon-1) = -yCon
  end if

  call SNOPT_evalHes0( m, n, nnH, neH, nlocH, iu(llocH), Hcol,        &
                       ru(lvalH), iu(lrowH), iu(lcolH), ru(lx), ru(ly) )

end subroutine SNOPT_evalHes

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine SNOPT_evalHes0( m, n, nnH, neH, nlocH, locH, Hcol,       &
                           valH, rowH, colH, x, y )
  use           CUTEst_interface_double
  implicit      none
  integer,      parameter     :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  integer(ip),  intent(in)  :: m, n, nnH, neH, nlocH
  real(rp),     intent(in)  :: y(m), x(n)

  integer(ip),  intent(out) :: locH(nlocH), rowH(neH), colH(neH)
  real(rp),     intent(out) :: Hcol(neH), valH(neH)

  !=============================================================================
  ! Does the work for SNOPT_evalHes.
  !
  !=============================================================================
  integer(ip)               :: j, k, l, lenH, nnzH, status

  ! Compute nnzH, valH, rowH, colH

  lenH = neH
  if (m > 0) then
     call cutest_csh( status, n, m, x, y, nnzH, lenH, valH, rowH, colH )
  else
     call cutest_ush( status, n,    x,    nnzH, lenH, valH, rowH, colH )
  end if
  if (status /= 0) go to 910

  ! Copy the Hessian from csh or ush, stored in (valH,rowH,colH)
  ! into the SNOPT Hessian stored by columns in Hcol.
  ! locH(j) points to the next available position in column j.

  do l = 1, nnzH
     j = colH(l)
     k = locH(j)
     Hcol(k) = valH(l)
     locH(j) = k + 1
  end do

  ! Reset locH.
  do j = nnH, 2, -1
     locH(j) = locH(j-1)
  end do
  locH(1) = 1

  return

910 continue
  write( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") status
  stop

end subroutine SNOPT_evalHes0
