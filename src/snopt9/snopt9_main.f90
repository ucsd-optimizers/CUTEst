module cutefun_module
  use snopt_module
  implicit none

  type, extends(snoptFunB) :: cuteFun
     integer(ip) :: n, m, ne
     integer(ip) :: neH, lenH

     integer(ip), allocatable :: iGfun(:), jGvar(:), locG(:)
     real(rp),    allocatable :: G(:)

     integer(ip), allocatable :: rowH(:), colH(:)
     real(rp),    allocatable :: valH(:), cx(:), cy(:)
   contains
     procedure, pass(this) :: funobj => SNOPT_evalfg
     procedure, pass(this) :: funcon => SNOPT_evalcj
  end type cuteFun

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lrhb_usrfun(mode, n, x, f, g, status)
    use snopt_module, only : ip, rp
    implicit none
    integer(ip),  intent(in)    :: n, status
    integer(ip),  intent(inout) :: mode
    real(rp),     intent(in)    :: x(n)
    real(rp),     intent(inout) :: f, g(n)
    !===========================================================================
    logical :: needG

    needG = .true.
    call cutest_uofg(status, n, x, f, g, needG)

  end subroutine lrhb_usrfun

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SNOPT_evalfg &
       (mode, this, nnObj, x, fObj, gObj, nState)
    class(cuteFun) :: this
    integer(ip),  intent(in)    :: nnObj, nState
    real(rp),     intent(in)    :: x(nnObj)
    integer(ip),  intent(inout) :: mode
    real(rp),     intent(inout) :: fObj, gObj(nnObj)

    !===========================================================================
    logical     :: needG
    integer(ip) :: status
    real(rp)    :: ggObj(this%n)

    if (mode == 0) then
       needG = .false.
    else
       needG = .true.
    end if

    this%cx(:) = 0.0
    this%cx(1:nnObj) = x(1:nnObj)
    call cutest_cofg (status, this%n, this%cx, fObj, ggObj, needG)

    if (needG) &
         gObj(:) = ggObj(1:nnObj)

    if (status /= 0) then
       write(6, "(' CUTEst error, status = ', i0, ', stopping')")  status
       stop
    end if

  end subroutine SNOPT_evalfg

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SNOPT_evalcj &
       (mode, this, nnCon, nnJac, negCon, x, fCon, gCon, nState)
    class(cuteFun) :: this
    integer(ip), intent(in)    :: nnJac, nnCon, negCon, nState
    real(rp),    intent(in)    :: x(nnJac)
    integer(ip), intent(inout) :: mode
    real(rp),    intent(inout) :: fcon(nnCon), gcon(negCon)

    !===========================================================================
    integer(ip) :: j, k, l, nnzJ, status
    logical     :: needG

    if (nnCon > 0) then
       needG  = mode > 0

       this%cx(:) = 0.0d+0
       this%cx(1:nnJac) = x(1:nnJac)

       ! Evaluate the problem constraints. On input, nnCon > 0.
       ! The Jacobian is stored in sparse format.
       call cutest_ccfsg &
            (status, this%n, this%m, this%cx, this%cy, nnzj, this%ne, &
             this%G, this%jGvar, this%iGfun, needG)
       if (status /= 0) go to 910

       fCon(1:nnCon) = this%cy(1:nnCon)

       if (needG) then
          ! Copy the Jacobian from CSCFG, stored in (G,iGfun,jGvar)
          ! into the SNOPT Jacobian stored by columns in gCon.
          ! locG(j) points to the next available position in column j.
          do l       = 1, nnzJ
             j       = this%jGvar(l)
             if (j <= nnJac .and. this%iGfun(l) <= nnCon) then
                k       = this%locG(j)
                gCon(k) = this%G(l)
                this%locG(j) = k + 1
             end if
          end do

          ! Reset locG.
          do j = nnJac, 2, -1
             this%locG(j) = this%locG(j-1)
          end do
          this%locG(1) = 1

       end if
    end if

    return

910 write(6, "(' CUTEst error, status = ', i0, ', stopping')") status
    stop

  end subroutine SNOPT_evalcj

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!$  subroutine SNOPT_getLagrangian (H, nnL, nnCon, x, y, nstate)
!!$    !use  SNOPT_mod, only : m0, n, lenH, valH, rowH, colH, cx, cy
!!$    class(snLj), target        :: H
!!$    integer(ip), intent(in)    :: nnL, nnCon, nstate
!!$    real(rp),    intent(in)    :: x(nnL), y(nnCon)
!!$
!!$    !===========================================================================
!!$    integer(ip) :: j, status, nnz, nnH
!!$
!!$    if (m0 > 0) then
!!$       cx(:) = 0.0
!!$       cy(:) = 0.0
!!$       cx(1:nnL) = x
!!$       cy(1:nnCon) = y
!!$
!!$       if (nnCon == 0) then
!!$          call cutest_cish(status, n, cx, 0, nnz, lenH, valH, rowH, colH)
!!$       else
!!$          call cutest_csh(status, n, m0, cx, cy, nnz, lenH, valH, rowH, colH)
!!$       end if
!!$
!!$    else
!!$       cx(:) = 0.0
!!$       cx(1:nnL) = x
!!$       call cutest_ush(status, n, cx, nnz, lenH, valH, rowH, colH)
!!$    end if
!!$
!!$    nnH = maxval(colH(1:nnz))
!!$    call populateH(H, nnH, nnz, valH, rowH, colH, .true.)
!!$    H%gotL = .true.
!!$
!!$  end subroutine SNOPT_getLagrangian

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module cutefun_module


program snopt_main
  use cutefun_module
  use snopt_module
  implicit none

  character     :: &
       Start*8, pname*10, filename*20
  logical       :: &
       noConstraints
  integer(ip)   :: &
       Errors, INFO, status, alloc_stat


  integer(ip),  pointer :: &
       m, n, neH, lenH, ne
  integer(ip) :: &
       nm, neJ, neG, nnH, iObj, nNames, nnCon, nnJac, nnObj, &
       nlc, neq, j, jslack, lj, nlocG
  real(rp)      :: &
       fObj, ObjAdd

  integer(ip)   :: &
       nFuns, nObjs, nCons, nfCon1, nfObj1
  real(rp)      :: &
       runTime, cpu(2), calls(7)

  type(cuteFun), target :: &
       funs
  type(snAj) :: &
       Jac     ! constraint matrix
  type(solver_data) :: &
       opts

  character(10), allocatable :: &
       cNames(:)
  character(8),  allocatable :: &
       Names(:)
  integer(ip), allocatable :: &
       hs(:), eType(:)
  real(rp), allocatable :: &
       bl(:), bu(:), x(:), pi(:), rc(:), g(:)
  real(rp), allocatable :: &
       c(:)
  logical, allocatable :: &
       equation(:), linear(:)

  integer(ip), parameter :: iCutest = 55,  iOut   = 6, io_buffer = 11
  real(rp),    parameter :: zero    = 0.0, infBnd = 1.0d+20

  !------------------------------------------------------------------------

  m    => funs%m
  n    => funs%n
  ne   => funs%ne
  neH  => funs%neH
  lenH => funs%lenH

  !-----------------------
  ! Open problem data file
  !-----------------------
  open  (iCutest, file = 'OUTSDIF.d', form = 'FORMATTED', status = 'old')
  rewind(iCutest)


  !----------------------------------
  ! Get dimensions and allocate space
  !----------------------------------
  call CUTEST_cdimen(status, iCutest, n, m)
  if (status /= 0) go to 910

  noConstraints = m == 0


  !-------------------------
  ! Start setting up problem
  !-------------------------
  if (noConstraints) then
     allocate &
          (bl(n), bu(n), x(n), g(n), Names(n), cNames(n), stat=alloc_stat)
     if (alloc_stat /= 0) GO TO 990

     call CUTEST_usetup &
          (status, iCutest, iOut, io_buffer, n, x, bl(1:n), bu(1:n))
     if (status /= 0) go to 910

     ! Get names
     call CUTEST_unames(status, n, pname, cNames(1:n))
     if (status /= 0) go to 910

     fObj = zero

  else
     nm = n + m
     allocate &
          (hs(nm+1), eType(nm+1), bl(nm+1), bu(nm+1), x(nm+1), pi(m+1), rc(nm+1), &
           Names(nm+1), cNames(nm+1), c(m+1), stat=alloc_stat)
     if (alloc_stat /= 0) GO TO 990

     allocate(equation(m), linear(m), stat=alloc_stat)
     call CUTEST_csetup &
          (status, iCutest, iOut, io_buffer, &
           n, m, x(1:n), bl(1:n), bu(1:n), pi, bl(n+1:nm), bu(n+1:nm), &
           equation, linear, 0, 2, 1)

     ! Fix bounds for equality constraints
     do j = 1, m
        if (equation(j)) then
           bl(n+j) = zero
           bu(n+j) = zero
        end if
     end do

     if (allocated(equation)) deallocate(equation)
     if (allocated(linear) ) deallocate(linear)
     if (status /= 0) go to 910

     ! Compute number of nonlinear variables,
     ! number of nonlinear constraints, and
     ! linear/equality constraints.
     call CUTEST_cstats (status, nnObj, nnJac, neq, nlc)
     nnCon = m - nlc
     if (status /= 0) go to 910


     ! Compute the objective and constraints at x = 0.
     rc(1:n) = zero
     call CUTEST_cfn (status, n, m, rc, fObj, c)
     if (status /= 0) go to 910


     ! Get names
     call CUTEST_cnames(status, n, m, pname, cNames(1:n), cNames(n+1:nm))
     if (status /= 0) go to 910


     ! Construct the structures for Jacobian.
     call cutest_cdimsj (status, ne)
     if (status /= 0) go to 910

     nlocG  = nnJac + 1

     allocate(funs%locG(nlocG), funs%iGfun(ne), funs%jGvar(ne), funs%G(ne))
     allocate(Jac%val(ne), Jac%loc(n+1), Jac%ind(ne))

     funs%locG(:)  = 0
     funs%iGfun(:) = 0
     funs%jGvar(:) = 0
     funs%G(:)     = zero

     lj    = ne
     call buildJac &
          (status, n, m, nnObj, nnCon, nnJac, lj, &
           iObj, Jac%loc, Jac%ind, Jac%val, &
           nlocG, funs%locG, funs%iGfun, funs%jGvar, funs%G, &
           x, pi, ne, neG)
     if (status /= 0) go to 910

     if (iObj > 0) then
        m       = m+1
        bl(n+m) = -infBnd
        bu(n+m) =  infBnd
        c(m)    =  zero

        Names(n+m) = 'OBJ     '
     end if

     Names(1:n+m) = cNames(1:n+m)
     nNames = n + m

     ! Set bounds on linear constraints.
     do j = 1, m
        jslack = n + j
        if (j > nnCon) then
           bu(jslack) = bu(jslack) - c(j)
           bl(jslack) = bl(jslack) - c(j)
        end if

        ! If possible, set slack variables to be nonbasic at zero.
        x(jslack) = max(zero, bl(jslack))
        x(jslack) = min(x(jslack), bu(jslack))
     end do


!!$     if (secondDeriv) then
!!$        ! Get structure of Hessian of Lagrangian
!!$        call cutest_cdimsh (status, lenH)
!!$        allocate (rowH(lenH), colH(lenH), valH(lenH))
!!$
!!$        call cutest_cshp (status, n, neH, lenH, rowH, colH)
!!$        nnH = maxval (colH(1:neH))
!!$     end if

     Jac%m   = m
     Jac%n   = n
     Jac%nnz = ne

     eType(:) = 0
     hs(:)    = 0
     pi(:)    = zero
     rc(:)    = zero

     if (nnObj == 0) ObjAdd = Objadd - fObj
  end if


  Start = 'cold'
  ObjAdd = zero

  !-----------------------------------------------------------------------------
  ! Ok, we're done with CUTEst stuff.
  !-----------------------------------------------------------------------------
  filename = trim(pname) // '.out'

  call snopt_initialize(filename, 'screen', opts)
  call snopt_specs(opts, 'SNOPT9.SPC', info)

  if (noConstraints) then
     ! LRHB
!!$     call snopt &
!!$          (pname(1:10), n, objAdd, x, bl, bu, lrhb_usrfun, fobj, g, info, opts)

     call CUTEST_ureport(status, calls, cpu)
     write(iOut, 2010) &
          pname, n, m, calls(1), calls(2), info, fObj, cpu(1), cpu(2)
     call CUTEST_uterminate(status)

  else
     ! SNOPT9
     allocate (funs%cx(n), funs%cy(m))

     ! Solve the problem
     Start = 'cold'

     call snopt &
          (start , pname(1:8), m, n, nnCon, nnObj, nnJac, ObjAdd, iObj, &
           funs, Jac, bl, bu, hs, eType, x, pi, rc, nNames, Names, &
           opts, info)

     call CUTEST_creport(status, calls, cpu)
     write(iOut, 2000) &
          pname, n, m, calls(1), calls(2), calls(5), calls(6), info, &
          fobj, cpu(1), cpu(2)
     call CUTEST_cterminate(status)
  end if
  call snopt_end(opts)


  ! Deallocate
  call Jac%trash
  if (allocated(funs%cx)) deallocate(funs%cx)
  if (allocated(funs%cy)) deallocate(funs%cy)

  if (allocated(funs%iGfun)) deallocate(funs%iGfun)
  if (allocated(funs%jGvar)) deallocate(funs%jGvar)
  if (allocated(funs%locG))  deallocate(funs%locG)
  if (allocated(funs%G))     deallocate(funs%G)

  if (allocated(Names)) deallocate(Names)
  if (allocated(cNames)) deallocate(cNames)
  if (allocated(Etype)) deallocate(eType)
  if (allocated(hs)   ) deallocate(hs)
  if (allocated(bl)   ) deallocate(bl)
  if (allocated(bu)   ) deallocate(bu)
  if (allocated(x)    ) deallocate(x)
  if (allocated(pi)   ) deallocate(pi)
  if (allocated(rc)   ) deallocate(rc)
  if (allocated(g)    ) deallocate(g)
  if (allocated(c)    ) deallocate(c)


  stop

910 continue
  write(iOut, "(' CUTEst error, status = ', i0, ', stopping')") status
  stop

990 continue
  write(iOut, "(' Allocation error, status = ', i0)") status
  stop

2000 format(/, 24('*'), ' CUTEst statistics ', 24('*') //                    &
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
         66('*') /)

2010 format(/, 24('*'), ' CUTEst statistics ', 24('*') //                    &
         ,' Package used            : SNOPT',    /                            &
         ,' Problem                 :  ', A10,    /                           &
         ,' # variables             =      ', I10 /                           &
         ,' # constraints           =      ', I10 /                           &
         ,' # objective functions   =        ', F8.2 /                        &
         ,' # objective gradients   =        ', F8.2 /                        &
         ,' Exit code               =      ', I10 /                           &
         ,' Final f                 = ', E15.7 /                              &
         ,' Set up time             =      ', 0P, F10.2, ' seconds' /         &
         ,' Solve time              =      ', 0P, F10.2, ' seconds' //        &
         66('*') /)
3000 format(a8, 3x, i2, 2x, f12.3, 3x, i12)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine buildJac &
       (status, n, m, nnObj, nnCon, nnJac, lj, &
        iObj, locJ, indJ, A, nlocG, locG, iGfun, jGvar, G, x, y, ne, neG)

    integer(ip), intent(in)  :: n, m, nnObj, nnCon, nnJac, lj, nlocG
    real(rp),    intent(in)  :: x(n), y(m)

    integer(ip), intent(out) :: iObj, status, ne, neG, locJ(n+1), indJ(lj), &
                                locG(nlocG), iGfun(lj), jGvar(lj)
    real(rp),    intent(out) :: A(lj), G(lj)

    !===========================================================================
    ! On entry, nlocG = nnJac + 1
    !===========================================================================
    integer(ip) :: i, j, k, l, neJ
    real(rp)    :: normObj


    call cutest_csgr &
         (status, n, m, x, y, .false., ne, lj, G, jGvar, iGfun)
    if (status /= 0) return

    ! Initialize
    locJ(1:n)     = 0
    locG(1:nnJac) = 0

    ! Find the norm of the constant objective gradients.
    ! If the norm is nonzero, define a linear objective row.
    normObj = 0.0
    do k = 1, ne
       j = jGvar(k)
       i = iGfun(k)
       if (i == 0  .and.  j > nnObj) then ! linear obj
          normObj = normObj + abs(G(k))
       end if
    end do

    ! normObj = 1.0d+0

    if (normObj > 0.0) then
       iObj = m+1
    else
       iObj = 0
    end if

    ! Count the Jacobian entries in column j.
    ! Initialize locJ.  Store the column counts in locJ(1:n).
    ! Don't include nonlinear objective function entries.
    neJ = ne
    do l = 1, ne
       j = jGvar(l)
       i = iGfun(l)
       if (i == 0 .and. (j <= nnObj .or. iObj == 0)) then
          neJ = neJ - 1
       else
          locJ(j) = locJ(j) + 1
       end if
    end do
    locJ(n+1) = neJ + 1

    ! Set locJ(j) to point to start of column j.
    do j = n, 1, -1
       locJ(j) = locJ(j+1) - locJ(j)
    end do

    ! Load nonlinear Jacobian entries into A and indJ.
    ! Use locJ to keep track of position for each variable b.
    ! Also count nonlinear Jacobian entries in each column and
    ! store count in locG(j).
    neG = 0
    do k = 1, ne
       j = jGvar(k)
       i = iGfun(k)
       if (i > 0 .and. i <= nnCon .and. j <= nnJac) then
          l       = locJ(j)
          A(l)    = G(k)
          indJ(l) = i
          locJ(j) = l   + 1
          locG(j) = locG(j) + 1
          neG     = neG + 1
       end if
    end do

    ! Load the constant Jacobian entries, including linear objective
    ! function entries,  in A and indJ.
    ! locJ(j) points to the next available position in column j.
    do k = 1, ne
       j = jGvar(k)
       i = iGfun(k)
       if (i == 0  .and.  j > nnObj .and. iObj > 0) then ! linear obj
          l       = locJ(j)
          A(l)    = G(k)
          indJ(l) = iObj
          locJ(j) = l + 1
       else if (i > nnCon .or. (i > 0 .and. j > nnJac)) then
          l       = locJ(j)
          A(l)    = G(k)
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

    ne = neJ

    if (ne == 0) then
       ne        = 1
       A(ne)     = 0.0
       indJ(ne)  = 1
       locJ(n+1) = ne + 1
    end if

  end subroutine buildJac

end program snopt_main
