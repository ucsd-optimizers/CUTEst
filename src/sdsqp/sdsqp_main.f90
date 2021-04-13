module cutefun_module
  use sdsqp_module
  implicit none

  type, extends(sdsqp_problem) :: cutestprob
     ! Arrays for cutest calls
     integer(ip), allocatable :: iGfun(:), jGvar(:), locG(:)
     real(rp),    allocatable :: G(:)
     integer(ip), allocatable :: rowH(:), colH(:)
     real(rp),    allocatable :: valH(:)
  end type cutestprob

  public

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_evalfg &
       (prob, mode, n, nnobj, x, fObj, gObj, status)
    class(sdsqp_problem) :: prob
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

    call cutest_cofg (stat, n, x, fObj, gObj, needG)

    if (stat /= 0) then
       write(6, "(' CUTEst error, status = ', i0, ', stopping')")  stat
       stop
    end if

  end subroutine SDSQP_evalfg

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_evalcj &
       (prob, mode, m, n, nncon, nnjac, ne, x, fCon, gCon, status)
    class(sdsqp_problem) :: prob
    integer(ip), intent(in)    :: m, n, nncon, nnjac, ne, status
    real(rp),    intent(in)    :: x(n)
    integer(ip), intent(inout) :: mode
    real(rp),    intent(inout) :: fcon(m)
    real(rp),    intent(inout) :: gCon(ne)
    !===========================================================================
    integer(ip) :: j, k, l, neg, leng, stat
    logical     :: needG

    if (m == 0) return

    needG  = mode > 0
    lenG   = ne

    select type(prob)
    type is (cutestprob)
       ! Evaluate the problem constraints.
       ! The Jacobian is stored in sparse format.
       call cutest_ccfsg &
            (stat, n, m, x, fCon, neg, leng, &
             prob%G, prob%jGvar, prob%iGfun, needG)
       if (stat /= 0) go to 910

       if (needG) then
          ! Copy the Jacobian from CSCFG, stored in (G,iGfun,jGvar)
          ! into the SNOPT Jacobian stored by columns in gCon.
          ! locG(j) points to the next available position in column j.
          do l = 1, neg
             j = prob%jGvar(l)
             k       = prob%locG(j)
             gCon(k) = prob%G(l)
             prob%locG(j) = k + 1
          end do

          ! Reset locG.
          do j = n, 2, -1
             prob%locG(j) = prob%locG(j-1)
          end do
          prob%locG(1) = 1
       end if
    end select

    return

910 write(6, "(' CUTEst error, status = ', i0, ', stopping')") stat
    stop

  end subroutine SDSQP_evalcj

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SDSQP_evalh &
       (prob, mode, m, n, nncon, nnjac, nnobj, ne, x, y, H, status)
    class(sdsqp_problem) :: prob
    integer(ip), intent(in)    :: m, n, nncon, nnjac, nnobj, ne, status
    real(rp),    intent(in)    :: x(n), y(m)
    integer(ip), intent(inout) :: mode
    real(rp),    intent(inout) :: H(ne)

    !===========================================================================
    integer(ip) :: stat, nnz, lenh

    lenh = ne

    select type(Hes => prob%H)
    type is (snHij)
       if (m > 0) then
!          if (nnCon > 0) then
             call cutest_csh(stat, n, m, x, -y, nnz, lenH, H, Hes%row, Hes%col)
             !          else
!             call cutest_cish(stat, n, x, 0, nnz, lenH, H, Hes%row, Hes%col)
!          end if

       else
          call cutest_ush(stat, n, x, nnz, lenH, H, Hes%row, Hes%col)
       end if
    end select

  end subroutine SDSQP_evalh

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

  integer(ip),  pointer :: &
       m, n, nnCon, nnjac, nnobj, neJ, neh
  integer(ip) :: &
       nm, nNames, nlc, neG, nejac, neq, iobj, nlocG
  real(rp)      :: &
       fObj, ObjAdd

  integer(ip)   :: &
       i, j, jslack, k, l, lenJ, lenH, &
       nFuns, nObjs, nCons, nfCon1, nfObj1
  real(rp)      :: &
       runTime, cpu(2), calls(7)

  type(cutestprob), target :: &
       myprob
  type(sdsqp_options) :: &
       opts

  character(len=10), allocatable :: &
       Names(:)
  integer(ip), allocatable :: &
       rowJ(:), colJ(:)
  real(rp), allocatable :: &
       c(:), valJ(:)
  logical, allocatable :: &
       equation(:), linear(:)

  integer(ip), parameter :: iCutest = 55,  iOut   = 6, io_buffer = 11
  real(rp),    parameter :: zero    = 0.0, infBnd = 1.0d+20

  !------------------------------------------------------------------------

  m     => myprob%m
  n     => myprob%n
  nncon => myprob%nncon
  nnjac => myprob%nnjac
  nnobj => myprob%nnobj
  nej   => myprob%nej
  neh   => myprob%neh

  myprob%funobj => SDSQP_evalfg
  myprob%funcon => SDSQP_evalcj
  myprob%funhes => SDSQP_evalh

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

  if (m == 0) go to 900


  !-------------------------
  ! Start setting up problem
  !-------------------------
  nm = n + m
  allocate &
       (myprob%bl(nm), myprob%bu(nm), myprob%x(nm), myprob%y(m), myprob%z(nm), &
        Names(nm), c(m), stat=alloc_stat)
  if (alloc_stat /= 0) GO TO 990

  allocate(equation(m), linear(m), stat=alloc_stat)
  call CUTEST_csetup &
       (status, iCutest, iOut, io_buffer, &
        n, m, myprob%x(1:n), myprob%bl(1:n), myprob%bu(1:n), myprob%y(1:m), &
        myprob%bl(n+1:nm), myprob%bu(n+1:nm), &
        equation, linear, 0, 2, 1)

  ! Fix bounds for equality constraints
  do j = 1, m
     if (equation(j)) then
        myprob%bl(n+j) = zero
        myprob%bu(n+j) = zero
     end if
  end do

  if (allocated(equation)) deallocate(equation)
  if (allocated(linear) ) deallocate(linear)
  if (status /= 0) go to 910

  ! Compute number of nonlinear variables,
  ! number of nonlinear constraints, and
  ! linear/equality constraints.
  call CUTEST_cstats (status, nnobj, nnjac, neq, nlc)
  !nnCon = m - nlc

  nncon = m
  nnjac = n
  nnobj = n
  if (status /= 0) go to 910


  ! Compute the objective and constraints at x = 0.
  myprob%z(1:nm) = zero
  call CUTEST_cfn (status, n, m, myprob%z, fObj, c)
  if (status /= 0) go to 910


  ! Get names
  call CUTEST_cnames(status, n, m, pname, Names(1:n), Names(n+1:nm))
  if (status /= 0) go to 910
  nNames = n + m


  ! Construct the structures for Jacobian.
  iobj = 0

  call cutest_cdimsj (status, lenj)
  if (status /= 0) go to 910

  ! Extra workspace for Jacobian
  nlocG  = nnJac + 1
  allocate(myprob%locG(nlocG), myprob%iGfun(lenj), myprob%jGvar(lenj), myprob%G(lenj))
  myprob%locG(:)  = 0
  myprob%iGfun(:) = 0
  myprob%jGvar(:) = 0
  myprob%G(:)     = zero

  ! Allocate Jacobian structure as sparse-by-col type
  allocate(snAj :: myprob%J)
  select type(Jac => myProb%J)
  type is (snAj)
     call Jac%init(m, n, lenj)
     allocate(valJ(lenj), rowJ(lenj), colJ(lenj))

     call cutest_csgr &
          (status, n, m, myprob%x, myprob%y, .false., neJ, lenJ, valJ, colJ, rowJ)
     if (status /= 0) return

     Jac%loc(1:n+1) = 0
     myprob%locG(1:nlocG) = 0

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

           myprob%locG(j) = myprob%locG(j) + 1
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
     myprob%locG(nlocG) = neG + 1

     do  j = n, 2, -1
        Jac%loc(j) = Jac%loc(j-1)
     end do
     Jac%loc(1) = 1

     do  j = nnJac, 2, -1
        myprob%locG(j) = myprob%locG(j+1) - myprob%locG(j)
     end do
     myprob%locG(1) = 1
     neJ = neJac
  end select

  ! Set bounds on linear constraints.
  do j = 1, m
     jslack = n + j
     if (j > nnCon) then
        myprob%bu(jslack) = myprob%bu(jslack) - c(j)
        myprob%bl(jslack) = myprob%bl(jslack) - c(j)
     end if

     ! If possible, set slack variables to be nonbasic at zero.
     myprob%x(jslack) = max(zero, myprob%bl(jslack))
     myprob%x(jslack) = min(myprob%x(jslack), myprob%bu(jslack))
  end do


  ! Get structure of Hessian of Lagrangian
  allocate(snHij :: myprob%H)
  select type(Hes => myprob%H)
  type is (snHij)
     call cutest_cdimsh (status, lenH)
     call Hes%init(n, lenH)
     call cutest_cshp (status, n, neH, lenH, Hes%row, Hes%col)
     myprob%nnH   = maxval(Hes%col(1:neH))
  end select

  !if (nnObj == 0) ObjAdd = Objadd - fObj


  !-----------------------------------------------------------------------------
  ! Ok, we're done with CUTEst stuff.
  !-----------------------------------------------------------------------------
  Start = 'cold'
  ObjAdd = zero

  filename = trim(pname) // '.out'

  call sdsqp_initialize(filename, 'screen', opts)
  call sdsqp_specs(opts, 'SDSQP.SPC', info)

  ! Solve the problem
  call sdsqp(myprob, opts, info)

  ! Problem stats
  call CUTEST_creport(status, calls, cpu)
  write(iOut, 2000) &
       pname, n, m, calls(1), calls(2), calls(5), calls(6), info, &
       fobj, cpu(1), cpu(2)
  call CUTEST_cterminate(status)


  call myprob%trash()
  call sdsqp_end(opts)

  ! Deallocate
  if (allocated(Names)) deallocate(Names)
  if (allocated(c))     deallocate(c)
  if (allocated(rowJ))  deallocate(rowJ)
  if (allocated(colJ))  deallocate(colJ)
  if (allocated(valJ))  deallocate(valJ)

  if (allocated(myprob%iGfun)) deallocate(myprob%iGfun)
  if (allocated(myprob%jGvar)) deallocate(myprob%jGvar)
  if (allocated(myprob%locG))  deallocate(myprob%locG)
  if (allocated(myprob%G))     deallocate(myprob%G)

  if (allocated(myprob%rowH))  deallocate(myprob%rowH)
  if (allocated(myprob%colH))  deallocate(myprob%colH)
  if (allocated(myprob%valH))  deallocate(myprob%valH)

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
