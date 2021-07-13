module cutefun_module
  use sdsqp_module
  implicit none

  type, extends(sdsqp_NLP) :: cutest_NLP
     ! Arrays for cutest calls
     real(rp),    allocatable :: bl(:), bu(:), x0(:), y0(:)

     integer(ip), allocatable :: iGfun(:), jGvar(:), locG(:)
     real(rp),    allocatable :: G(:)
  end type cutest_NLP

  type, extends(sdsqp_QP) :: cutest_QP
     real(rp),    allocatable :: bl(:), bu(:), x0(:), y0(:)
  end type cutest_QP

  public

  logical, public, parameter :: isQP = .false.
  !logical, public, parameter :: isQP = .true.
contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_bounds(prob, m, n, bl, bu)
    class(sdsqp_problem) :: prob
    integer(ip), intent(in)  :: m, n
    real(rp),    intent(out) :: bl(n+m), bu(n+m)
    !===========================================================================
    !===========================================================================

    select type(prob)
    type is (cutest_QP)
       bl(1:n+m) = prob%bl(1:n+m)
       bu(1:n+m) = prob%bu(1:n+m)
    type is (cutest_NLP)
       bl(1:n+m) = prob%bl(1:n+m)
       bu(1:n+m) = prob%bu(1:n+m)
    end select

  end subroutine SDSQP_bounds

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_start(prob, m, n, x, y)
    class(sdsqp_problem) :: prob
    integer(ip), intent(in)  :: m, n
    real(rp),    intent(out) :: x(n+m), y(m)
    !===========================================================================
    !===========================================================================

    select type(prob)
    type is (cutest_QP)
       x(1:n+m) = prob%x0(1:n+m)
       y(1:m)   = prob%y0(1:m)
    type is (cutest_NLP)
       x(1:n+m) = prob%x0(1:n+m)
       y(1:m)   = prob%y0(1:m)
    end select

  end subroutine SDSQP_start

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_jacobian &
       (prob, m, n, nncon, nnjac, nnobj, neJ, Jac, neG)
    class(sdsqp_NLP) :: prob
    integer(ip), intent(in)  :: m, n, nncon, nnjac, nnobj
    integer(ip), intent(out) :: neJ
    integer(ip), intent(out), optional :: neG
    class(snA),  intent(out), allocatable :: Jac
    !===========================================================================
    ! Return the structure for the Jacobian matrix.
    !===========================================================================
    integer(ip) :: i, j, k, l, lenj, nejac, status, iobj, nlocg
    integer(ip), allocatable :: &
         rowJ(:), colJ(:)
    real(rp), allocatable :: &
         valJ(:), fcon(:)

    iobj = 0
    if (m > 0) then

       ! Get Jacobian structure from CUTEst
       call cutest_cdimsj (status, lenj)
       if (status /= 0) return

       ! Allocate Jacobian structure
       allocate(snAj :: Jac)
       call Jac%init(m, n, lenj)

       select type(prob)
       type is (cutest_NLP)

!!$       select type(Jac)
!!$       class is (snAij)
!!$          call cutest_ccfsg &
!!$               (status, n, m, prob%x0, fcon, nej, lenj, &
!!$               Jac%val, Jac%col, Jac%row, .true.)

          select type(Jac)
          class is (snAj)

             ! Extra workspace for Jacobian
             nlocG  = nnJac + 1
             allocate(prob%locG(nlocG), prob%iGfun(lenj), prob%jGvar(lenj), prob%G(lenj))
             prob%locG(:)  = 0
             prob%iGfun(:) = 0
             prob%jGvar(:) = 0
             prob%G(:)     = 0.0d+0

             allocate(valJ(lenj), rowJ(lenj), colJ(lenj), fcon(m))

             call cutest_ccfsg &
                  (status, n, m, prob%x0, fcon, nej, lenj, &
                  valJ, colJ, rowJ, .true.)

!!$     call cutest_csgr &
!!$          (status, n, m, prob%x0, prob%y0, .false., neJ, lenJ, valJ, colJ, rowJ)
             if (status /= 0) return

             Jac%loc(1:n+1) = 0
             prob%locG(1:nlocG) = 0

             ! Count the Jacobian entries in column j.
             ! Initialize locJ.  Store the column counts in locJ(1:n).
             ! Don't include nonlinear objective function entries.
             nejac = nej
             do k  = 1, neJ
                j  = colJ(k)
                i  = rowJ(k)
                if (i == 0 .and. (j <= nnObj .or. iObj == 0)) then
                   neJac   = neJac   - 1
                else
                   Jac%loc(j) = Jac%loc(j) + 1
                end if
             end do
             Jac%loc(n+1) = neJac + 1

             ! Set locJ(j) to point to start of column j.
             do j = n, 1, -1
                Jac%loc(j) = Jac%loc(j+1) - Jac%loc(j)
             end do

             ! Load nonlinear Jacobian entries into Jcol and indJ.
             ! Use locJ to keep track of position for each variable b.
             ! Also count nonlinear Jacobian entries in each column and
             ! store count in locG(j).
             neG = 0
             do k = 1, neJ
                j = colJ(k)
                i = rowJ(k)
                if (i > 0 .and. i <= nnCon .and. j <= nnJac) then
                   l       = Jac%loc(j)
                   Jac%val(l) = valJ(k)
                   Jac%ind(l) = i
                   Jac%loc(j) = l   + 1

                   prob%locG(j) = prob%locG(j) + 1
                   neG     = neG + 1
                end if
             end do

             ! Load the constant Jacobian entries, including linear objective
             ! function entries,  in Jcol and indJ.
             ! locJ(j) points to the next available position in column j.
             do k = 1, neJ
                j = colJ(k)
                i = rowJ(k)
                if (i == 0 .and. j > nnObj .and. iObj > 0) then ! linear obj
                   l       = Jac%loc(j)
                   Jac%val(l) = valJ(k)
                   Jac%ind(l) = iObj
                   Jac%loc(j) = l + 1
                else if (i > nnCon .or. (i > 0 .and. j > nnJac)) then
                   l          = Jac%loc(j)
                   Jac%val(l) = valJ(k)
                   Jac%ind(l) = i
                   Jac%loc(j) = l + 1
                end if
             end do

             ! Reset locJ  and set locG.
             ! locG(j)  points  to the start of column j in Gcon,
             ! the nonlinear Jacobian.  These pointers are used by usrfgh.
             prob%locG(nlocG) = neG + 1

             do  j = n, 2, -1
                Jac%loc(j) = Jac%loc(j-1)
             end do
             Jac%loc(1) = 1

             do  j = nnJac, 2, -1
                prob%locG(j) = prob%locG(j+1) - prob%locG(j)
             end do
             prob%locG(1) = 1
             neJ = neJac
          end select
       end select
       neG = 0

    end if

    if (allocated(valJ)) deallocate(valJ)
    if (allocated(rowJ)) deallocate(rowJ)
    if (allocated(colJ)) deallocate(colJ)
    if (allocated(fcon)) deallocate(fcon)

  end subroutine SDSQP_jacobian

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_hessian &
       (prob, m, n, nncon, nnjac, nnobj, nnH, neH, H)
    class(sdsqp_NLP) :: prob
    integer(ip), intent(in)  :: m, n, nncon, nnjac, nnobj
    integer(ip), intent(out) :: neH, nnH
    class(snH),  intent(out), allocatable :: H

    !===========================================================================
    ! Get structure of Hessian of Lagrangian
    !===========================================================================
    integer(ip) :: lenH, status

    if (m > 0) then
       call cutest_cdimsh(status, lenH)
    else
       call cutest_udimsh(status, lenH)
    end if

    allocate(snHij :: H)
    call H%init(n, lenH)

    select type(H)
    type is (snHij)
       if (m > 0) then
          call cutest_cshp(status, n, neH, lenH, H%col, H%row)
       else
          call cutest_ushp(status, n, neH, lenH, H%col, H%row)
       end if
       nnH = maxval(H%col(1:neH))
    end select

  end subroutine SDSQP_hessian

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_evalfg &
       (prob, mode, status, n, nnobj, x, fObj, gObj)
    class(sdsqp_NLP) :: prob
    integer(ip),  intent(in)    :: n, nnobj, status
    real(rp),     intent(in)    :: x(n)
    integer(ip),  intent(inout) :: mode
    real(rp),     intent(inout) :: fObj, gObj(n)

    !===========================================================================
    logical     :: needG
    integer(ip) :: stat

    if (mode == 0) then
       needG = .false.
    else
       needG = .true.
    end if

    if (prob%m > 0) then
       call cutest_cofg (stat, n, x, fObj, gObj, needG)
    else
       call cutest_uofg(stat, n, x, fobj, gobj, needG)
    end if

    if (stat /= 0) then
       write(6, "(' CUTEst error, status = ', i0, ', stopping')")  stat
       stop
    end if

  end subroutine SDSQP_evalfg

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_evalcj &
       (prob, mode, status, m, n, nncon, nnjac, ne, x, fcon, Jac)
    class(sdsqp_NLP) :: prob
    integer(ip), intent(in)    :: m, n, nncon, nnjac, ne, status
    real(rp),    intent(in)    :: x(n)
    integer(ip), intent(inout) :: mode
    real(rp),    intent(inout) :: fcon(m)
    class(snA),  intent(inout) :: Jac
    !===========================================================================
    integer(ip) :: j, k, l, neg, leng, stat
    logical     :: needG

    if (m == 0) return

    needG  = mode > 0
    lenG   = ne

    select type(prob)
    type is (cutest_NLP)
       ! Evaluate the problem constraints.
       ! The Jacobian is stored in sparse format.
       call cutest_ccfsg &
            (stat, n, m, x, fCon, neg, leng, &
             prob%G, prob%jGvar, prob%iGfun, needG)

!!$       call cutest_ccfsg &
!!$            (stat, n, m, x, fCon, neg, leng, &
!!$             Jac%val, Jac%col, Jac%row, needG)
       if (stat /= 0) go to 910

    select type(Jac)
    type is (snAj)
       if (needG) then
          ! Copy the Jacobian from CSCFG, stored in (G,iGfun,jGvar)
          ! into the SNOPT Jacobian stored by columns in gCon.
          ! locG(j) points to the next available position in column j.
          do l = 1, neg
             j = prob%jGvar(l)
             k = prob%locG(j)
             Jac%val(k) = prob%G(l)
             prob%locG(j) = k + 1
          end do

          ! Reset locG.
          do j = n, 2, -1
             prob%locG(j) = prob%locG(j-1)
          end do
          prob%locG(1) = 1
       end if

    end select
    end select

    return

910 write(6, "(' CUTEst error, status = ', i0, ', stopping')") stat
    stop

  end subroutine SDSQP_evalcj

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_evalh &
       (prob, mode, status, m, n, nncon, nnjac, nnobj, ne, x, y, H)
    class(sdsqp_NLP) :: prob
    integer(ip), intent(in)    :: m, n, nncon, nnjac, nnobj, ne, status
    real(rp),    intent(in)    :: x(n), y(m)
    integer(ip), intent(inout) :: mode
    class(snH),  intent(inout) :: H

    !===========================================================================
    integer(ip) :: stat, nnz, lenh

    lenh = ne

    select type(H)
    type is (snHij)
       if (m > 0) then
!          if (nnCon > 0) then
          call cutest_csh(stat, n, m, x, -y, nnz, lenH, H%val, H%col, H%row)
!          else
!             call cutest_cish(stat, n, x, 0, nnz, lenH, H%val, H%col, H%row)
!          end if
       else
          call cutest_ush(stat, n, x, nnz, lenH, H%val, H%col, H%row)
       end if
    end select

  end subroutine SDSQP_evalh

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_constraintsQP &
       (prob, m, n, neJ, Jac, lenb, b)
    class(sdsqp_QP)     :: prob
    integer(ip), intent(in)  :: m, n
    integer(ip), intent(out) :: neJ, lenb
    real(rp),    intent(out) :: b(m)
    class(snA),  intent(out), allocatable :: Jac
    !===========================================================================
    integer(ip) :: i, j, k, ii, lenj, status
    real(rp)    :: x(n)
    integer(ip), allocatable :: row(:), col(:)
    real(rp),    allocatable :: val(:)

    ! Get Jacobian structure from CUTEst
    call cutest_cdimsj (status, lenj)

    allocate(row(lenj), col(lenj), val(lenj))
    x(1:n) = 0.0d+0

    lenb = 0
    call CUTEST_ccfsg &
         (status, n, m, x, b, nej, lenj, val, col, row, .true.)
    !b(1:m) = -b(1:m)

    ! Allocate Jacobian structure
    allocate(snAj :: Jac)
    call Jac%init(m, n, nej)

    select type(Jac)
    class is (snAj)
       ! Convert to sparse-by-col
       ! First count up elements in each column and store in loc.
       Jac%loc(1:n+1) = 0
       do k = 1, nej
          j = col(k)
          i = row(k)
          Jac%loc(j) = Jac%loc(j) + 1
       end do
       Jac%loc(n+1) = nej + 1

       ! Set the column pointers.
       do k = n, 1, -1
          Jac%loc(k) = Jac%loc(k+1) - Jac%loc(k)
       end do

       ! Put the row indices and values in and let loc(j) track the position of the
       ! elements for the jth column.
       Jac%ind(1:nej) = 0
       Jac%val(1:nej) = 0.0d+0
       do k = 1, nej
          j  = col(k)
          i  = row(k)
          ii = Jac%loc(j)

          Jac%val(ii) = val(k)
          Jac%ind(ii) = i
          Jac%loc(j)  = ii + 1
       end do

       ! Reset the column pointers.
       Jac%loc(2:n) = Jac%loc(1:n-1)
       Jac%loc(1)   = 1
    end select

    deallocate(val, row, col)

  end subroutine SDSQP_constraintsQP

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_objectiveQP &
       (prob, m, n, nnH, neH, H, fobj, cobj)

    class(sdsqp_QP) :: prob
    integer(ip), intent(in)  :: m, n
    integer(ip), intent(out) :: nnH, neH
    real(rp),    intent(out) :: fobj, cobj(n)
    class(snH),  intent(out), allocatable :: H
    !===========================================================================
    integer(ip) :: lenh, status
    real(rp) :: x(n)

    ! Hessian
    call cutest_cdimsh(status, lenh)

    allocate(snHij :: H)
    call H%init(nnH, lenH)

    select type(H)
    class is (snHij)
       ! Swap row/col for lower-tri
       x(1:n) = 0.0d+0
       call CUTEST_cish(status, n, x, 0, neH, lenh, H%val, H%col, H%row)

    end select
    !nnH = maxval(H%col(1:neH))
    nnH = n

    ! Objective constant and gradient
    x(1:n) = 0.0d+0
    call CUTEST_cofg(status, n, x, fobj, cObj, .true.)

  end subroutine SDSQP_objectiveQP

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module cutefun_module

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program sdsqp_main
  use cutefun_module
  use sdsqp_module
  implicit none

  character     :: &
       Start*8, pname*10, filename*20
  integer(ip)   :: &
       Errors, INFO, status, alloc_stat

  integer(ip) :: &
       m, n, nnCon, nnjac, nnobj
  integer(ip) :: &
       nm, nNames, nlc, neG, nejac, neq
  real(rp)      :: &
       fObj, ObjAdd

  integer(ip)   :: &
       i, j, jslack, k, l, lenJ, lenH, &
       nFuns, nObjs, nCons, nfCon1, nfObj1
  real(rp)      :: &
       runTime, cpu(2), calls(7)

  type(cutest_NLP), target :: &
       prob
!!$  type(cutest_QP), target :: &
!!$       prob

  type(sdsqp_options) :: &
       opts

  character(len=10), allocatable :: &
       Names(:)
  real(rp), allocatable :: &
       x(:), c(:), z0(:)
  logical, allocatable :: &
       equation(:), linear(:)

  integer(ip), parameter :: iCutest = 55,  iOut   = 6, io_buffer = 11
  real(rp),    parameter :: zero    = 0.0, infBnd = 1.0d+20

  !------------------------------------------------------------------------

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

!  if (m == 0) go to 900
  if (m > 1000 .or. n > 1000) go to 900

  !-------------------------
  ! Start setting up problem
  !-------------------------
  if (m > 0) then
     nm = n + m
     allocate &
          (prob%bl(nm), prob%bu(nm), prob%x0(nm), prob%y0(m), &
           z0(nm), Names(nm), c(m), stat=alloc_stat)
     if (alloc_stat /= 0) GO TO 990
     allocate(equation(m), linear(m), stat=alloc_stat)

     call CUTEST_csetup &
          (status, iCutest, iOut, io_buffer, &
           n, m, prob%x0(1:n), prob%bl(1:n), prob%bu(1:n), prob%y0(1:m), &
           prob%bl(n+1:nm), prob%bu(n+1:nm), &
           equation, linear, 0, 2, 1)

     ! Fix bounds for equality constraints
     do j = 1, m
        if (equation(j)) then
           prob%bl(n+j) = zero
           prob%bu(n+j) = zero
        end if
     end do

     if (allocated(equation)) deallocate(equation)
     if (allocated(linear) ) deallocate(linear)
     if (status /= 0) go to 910


     ! Compute number of nonlinear variables, number of nonlinear constraints, and
     ! linear/equality constraints.
     call CUTEST_cstats (status, nnobj, nnjac, neq, nlc)
     if (status /= 0) go to 910

     !nnCon = m - nlc
     nncon = m
     nnjac = n
     nnobj = n

     ! Compute the objective and constraints at x = 0.
     z0(1:nm) = zero
     call CUTEST_cfn (status, n, m, z0, fObj, c)
     if (status /= 0) go to 910

     ! Set bounds on linear constraints.
     do j = 1, m
        jslack = n + j
        !     if (j > nnCon) then! .or. isQP) then
        if (j > nnCon .or. isQP) then
           prob%bu(jslack) = prob%bu(jslack) - c(j)
           prob%bl(jslack) = prob%bl(jslack) - c(j)
        end if

        ! If possible, set slack variables to be nonbasic at zero.
        prob%x0(jslack) = max(zero, prob%bl(jslack))
        prob%x0(jslack) = min(prob%x0(jslack), prob%bu(jslack))
     end do

     !if (nnObj == 0) ObjAdd = Objadd - fObj

     ! Get names
     call CUTEST_cnames(status, n, m, pname, Names(1:n), Names(n+1:nm))
     if (status /= 0) go to 910
     nNames = n + m

  else
     ! unconstrained
     nm = n
     allocate &
          (prob%bl(n), prob%bu(n), prob%x0(n), z0(n), Names(n), c(n), stat=alloc_stat)
     if (alloc_stat /= 0) GO TO 990

     call CUTEST_usetup &
          (status, iCutest, iOut, io_buffer, &
           n, prob%x0(1:n), prob%bl(1:n), prob%bu(1:n))

     call CUTEST_unames(status, n, pname, Names(1:n))
     if (status /= 0) go to 910
     nNames = n

     nncon = 0
     nnjac = n
     nnobj = n

     z0(1:n) = zero
     call cutest_uofg(status, n, z0, fobj, c, .false.)
  end if

  prob%name  = trim(pname)
  prob%m     = m
  prob%n     = n

  prob%get_bounds   => SDSQP_bounds
  prob%get_start    => SDSQP_start

  if (isQP) then
!!$     prob%get_jacobian => SDSQP_constraintsQP
!!$     prob%get_hessian  => SDSQP_objectiveQP
  else
     !prob%nncon = nncon
     prob%nncon = m
     prob%nnjac = nnjac
     prob%nnobj = nnobj

     prob%get_jacobian => SDSQP_jacobian
     prob%get_hessian  => SDSQP_hessian
     prob%eval_obj     => SDSQP_evalfg
     prob%eval_con     => SDSQP_evalcj
     prob%eval_hes     => SDSQP_evalh
  end if

  !-----------------------------------------------------------------------------
  ! Ok, we're done with CUTEst stuff.
  !-----------------------------------------------------------------------------
  Start = 'cold'
  ObjAdd = zero

  filename = trim(pname) // '.out'
  call sdsqp_initialize(filename, 'screen', opts)
  call sdsqp_specs(opts, 'SDSQP.SPC', info)

  ! Solve the problem
  call sdsqp(prob, opts, info, fobj, x)

  ! Problem stats
  call CUTEST_creport(status, calls, cpu)
  write(iOut, 2000) &
       pname, n, m, calls(1), calls(2), calls(5), calls(6), info, &
       fobj, cpu(1), cpu(2)
  call CUTEST_cterminate(status)

  call sdsqp_end(opts)

  ! Deallocate
  if (allocated(Names)) deallocate(Names)
  if (allocated(c))     deallocate(c)
  if (allocated(x))     deallocate(x)
  if (allocated(z0))    deallocate(z0)

  if (allocated(prob%bl)) deallocate(prob%bl)
  if (allocated(prob%bu)) deallocate(prob%bu)
  if (allocated(prob%x0)) deallocate(prob%x0)
  if (allocated(prob%y0)) deallocate(prob%y0)

!!$  if (allocated(prob%iGfun)) deallocate(prob%iGfun)
!!$  if (allocated(prob%jGvar)) deallocate(prob%jGvar)
!!$  if (allocated(prob%locG))  deallocate(prob%locG)
!!$  if (allocated(prob%G))     deallocate(prob%G)

900 stop

910 continue
  write(iOut, "(' CUTEst error, status = ', i0, ', stopping')") status
  stop

990 continue
  write(iOut, "(' Allocation error, status = ', i0)") status
  stop

2000 format(/, 24('*'), ' CUTEst statistics ', 24('*') //                    &
         ,' Package used            : SDSQP',    /                            &
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

end program sdsqp_main
