program sqic_main
  use sqic_module

  implicit none

  integer(ip), parameter :: iCutest = 55, iOut = 6, io_buffer = 11

  integer(ip)   :: status, alloc_stat, j, k, i
  real(rp)      :: cpu(2), calls(7)

  real(rp),      allocatable :: zero(:), b(:)
  logical,       allocatable :: equation(:), linear(:)

  ! SQIC variables
  integer(ip)   :: INFO, n, m, nm, ncObj, neA, lenA, neH, lenH
  integer(ip)   :: iObj, nNames
  real(rp)      :: Obj, ObjAdd

  character(10) :: Prob
  character(20) :: filename

  character(10), allocatable :: cNames(:)
  character(8),  allocatable :: Names(:)
  integer(ip),   allocatable :: hs(:), hEtype(:)
  real(rp),      allocatable :: bl(:), bu(:), x(:), pi(:), rc(:), cObj(:)

  type(solver_data) :: probdat
  type(snHij)       :: H    ! coordinate H
  type(snAij)       :: A    ! coordinate A

  !-----------------------------------------------------------------------------
  ! Read CUTEst problem data
  open(iCutest, file = 'OUTSDIF.d', form = 'FORMATTED', status = 'old')
  rewind(iCutest)

  call CUTEST_cdimen(status, iCutest, n, m)  ! Problem dimension
  if (status /= 0) go to 910

  if (m == 0) then
     return
  end if


  nm = n+m
  allocate &
       (bl(nm), bu(nm), cObj(n), x(nm), hs(nm), pi(m), rc(nm), &
        hEtype(nm), b(m), zero(n), equation(m), linear(m),     &
        cNames(nm), Names(nm),                                 &
        stat = alloc_stat)
  if (alloc_stat /= 0) GO TO 990

  hEtype(:) = 0
  hs(:)     = 0
  zero(:)   = 0.0
  x(:)      = 0.0

  ! Get initial x, bounds, constraint types
  call CUTEST_csetup &
       (status, iCutest, iOut, io_buffer, &
        n, m, x(1:n), bl(1:n), bu(1:n),   &
        pi, bl(n+1:nm), bu(n+1:nm), equation, linear, 0, 2, 1)
  if (status /= 0) go to 910

  ! Check that all constraints are linear
  ! ...

  deallocate(equation, linear)


  ! Variable/constraint names
  call CUTEST_cnames(status, n, m, Prob, cNames(1:n), cNames(n+1:n+m))
  if (status /= 0) go to 910

  nNames = n + m
  Names(1:n+m) = cNames(1:n+m)
  deallocate(cNames)


  ! Constraint matrix
  call CUTEST_cdimsj(status, lenA)
  if (status /= 0) go to 910

  lenA = lenA - n

  allocate (A%cval(lenA), A%row(lenA), A%col(lenA), stat = alloc_stat)
  if (alloc_stat /= 0) GO TO 990

  call CUTEST_ccfsg(status, n, m, zero, b, neA, lenA, A%cval, A%col, A%row, .true.)
  if (status /= 0) go to 910

  A%m   = m
  A%n   = n
  A%nnz = neA


  ! Adjust bounds
  bl(n+1:nm) = bl(n+1:nm) - b
  bu(n+1:nm) = bu(n+1:nm) - b
  deallocate(b)


  ! Objective
  ! cofg returns ObjAdd and cObj.
  call CUTEST_cofg(status, n, zero, ObjAdd, cObj, .true.)
  if (status /= 0) go to 910


  ! Hessian of the objective
  call CUTEST_cdimsh(status, lenH)
  if ( status /= 0 ) go to 910

  if (lenH > 0) then
     allocate(H%row(2*lenH), H%col(2*lenH), H%cval(2*lenH), stat = alloc_stat)
     if (alloc_stat /= 0) GO TO 990

     call CUTEST_cish(status, n, zero, 0, neH, lenH, H%cval, H%row, H%col)
     if (status /= 0) go to 910

     i = 0
     do k = 1, neH
        if (H%row(k) < H%col(k)) then
           i = i + 1
           H%row(neH+i)  = H%col(k)
           H%col(neH+i)  = H%row(k)
           H%cval(neH+i) = H%cval(k)
        end if
     end do
     neH = neH + i

     H%n   = maxval(H%col(1:neH))
     H%nnz = neH

  else
     H%n   = 0
     H%nnz = 0
  end if


  !-----------------------------------------------------------------------------
  ! Ok, we're done with CUTEst stuff.
  !-----------------------------------------------------------------------------
  iObj  = 0
  ncObj = n

  filename = trim(Prob)//'.out'
  call sqic_initialize(filename, 'screen', probdat)
  call sqic_specs(probdat, 'SQIC.SPC', info)

  call sqic &
       ('Cold', Prob(1:8), m, n, H, ncObj, cObj, ObjAdd, iObj, &
        A, bl, bu, hs, hEtype, x, pi, rc, nNames, Names, &
        probdat, INFO)
  call sqic_end(probdat)


  call CUTEST_creport(status, CALLS, CPU)
  WRITE (iOut, 2000) Prob, n, m, CALLS(1), CALLS(2), &
                      CALLS(5), CALLS(6), info, Obj, CPU(1), CPU(2)


  ! Deallocate
  deallocate(Names)
  deallocate(hs)
  deallocate(hEtype)
  deallocate(cObj)
  deallocate(bl)
  deallocate(bu)
  deallocate(x)
  deallocate(pi)
  deallocate(rc)

  call A%trash
  call H%trash

  close (iCutest)

  call CUTEST_cterminate(status)
  stop

910 CONTINUE
  WRITE(iOut, "(' CUTEst error, status = ', i0, ', stopping')") status
  stop

990 CONTINUE
  WRITE(iOut, "(' Allocation error, status = ', i0)") status
  stop

2000 FORMAT(/, 24('*'), ' CUTEst statistics ', 24('*') //                    &
          ,' Package used            :  SQIC',    /                            &
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

end program sqic_main
