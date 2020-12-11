program sqopt_main
  implicit none

  integer, parameter :: iCutest = 55, iOut = 6, io_buffer = 11

  integer   :: status, alloc_stat, j, k, i
  double precision      :: cpu(2), calls(7)

  double precision,    allocatable :: zero(:), b(:), Hx(:)
  logical,     allocatable :: equation(:), linear(:)

  ! SQOPT variables
  integer :: &
       INFO, n, m, nm, neA, nNames, ncObj, iObj, lenA, nnH, lenH, &
       nS, nInf, mincw, miniw, minrw, &
       iPrint, iSumm, Errors
  double precision :: &
       ObjQP, fObj, sInf

  integer   :: lencw, leniw, lenrw, lencu, leniu, lenru
  character(len=8), allocatable :: &
       cw0(:), cw(:), cu(:), Names(:)
  integer, allocatable :: &
       eType(:), hs(:), indA(:), locA(:), row(:), col(:), iw0(:), iw(:), iu(:)
  double precision,    allocatable :: &
       bl(:), bu(:), x(:), pi(:), rc(:), cObj(:), &
       val(:), valA(:), rw0(:), rw(:), ru(:)

  character(len=10) :: probname
  character(len=20) :: filename

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
  nNames = n+m
  ncObj = n
  nnH = n

  ! Allocate SQOPT structures
  allocate &
       (eType(nm), hs(nm), bl(nm), bu(nm), cObj(n), x(nm), &
        pi(m), rc(nm), names(n+m), &
        stat=alloc_stat)
  allocate &
       (zero(n), equation(m), linear(m), &
        stat = alloc_stat)
  if (alloc_stat /= 0) GO TO 990

  zero(1:n)   = 0.0
  x(1:n+m) = 0.0

  ! Get initial x, bounds, constraint types
  call CUTEST_csetup &
       (status, iCutest, iOut, io_buffer, &
        n, m, x(1:n), bl(1:n), bu(1:n), &
        pi, bl(n+1:nm), bu(n+1:nm), &
        equation, linear, 0, 2, 1)
  if (status /= 0) go to 910

  ! Check that all constraints are linear
  ! ...

  deallocate(equation, linear)
  close (iCutest)

  ! Variable/constraint names
  call CUTEST_cnames(status, n, m, probname, names(1:n), names(n+1:n+m))
  if (status /= 0) go to 910

  ! Constraint matrix
  call CUTEST_cdimsj(status, lenA)
  if (status /= 0) go to 910

  lenA = lenA - n

  allocate(b(m), stat = alloc_stat)
  if (alloc_stat /= 0) GO TO 990

  allocate(row(lenA), col(lenA), val(lenA))
  call CUTEST_ccfsg &
       (status, n, m, zero, b, neA, lenA, val, col, row, .true.)
  if (status /= 0) go to 910

  ! Convert to sparse format
  allocate(indA(neA), valA(neA), locA(n+1))
  call crd2spr(m, n, neA, row, col, val, indA, locA, valA)

  ! Adjust bounds
  bl(n+1:nm) = bl(n+1:nm) - b(1:m)
  bu(n+1:nm) = bu(n+1:nm) - b(1:m)
  deallocate(b)

  ! Objective
  ! cofg returns ObjAdd and cObj.
  call CUTEST_cofg(status, n, zero, fObj, cObj, .true.)
  if (status /= 0) go to 910

  ! Hessian of the objective
  call CUTEST_cdimsh(status, lenH)
  if ( status /= 0 ) go to 910

  deallocate(zero)

  !-----------------------------------------------------------------------------
  ! Ok, we're done with CUTEst stuff.
  !-----------------------------------------------------------------------------
  lencw = 500
  leniw = 500
  lenrw = 500
  allocate(cw(lencw), iw(leniw), rw(lenrw))

  lencu = 1
  leniu = 1
  lenru = 1
  allocate(cu(lencu), iu(leniu), ru(lenru))

  iObj = 0

  filename = trim(probname)//'.out'
  call sqInitF(filename, 'screen', iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw)
  call sqSpecF('SQOPT.SPC', info, cw, lencw, iw, leniw, rw, lenrw)

  call sqMem &
       (info, m, n, neA, ncObj, nnH, &
        mincw, miniw, minrw, &
        cw, lencw, iw, leniw, rw, lenrw)

  if (leniw <= miniw) then
     leniw = miniw
     allocate(iw0(leniw))
     iw0(1:500) = iw(1:500)
     call move_alloc(from=iw0, to=iw)

     call sqSeti &
          ('Total integer workspace', leniw, 0, 0, Errors, &
           cw, lencw, iw, leniw, rw, lenrw )
  end if

  if (lenrw <= minrw) then
     lenrw = minrw
     allocate(rw0(lenrw))
     rw0(1:500) = rw(1:500)
     call move_alloc(from=rw0, to=rw)

     call sqSeti &
          ('Total real workspace', lenrw, 0, 0, Errors, &
          cw, lencw, iw, leniw, rw, lenrw )
  end if

  call sqopt &
       ('Cold', cutestHx, m, n, neA, nNames, ncObj, nnH, iObj, &
        fObj, probname, valA, indA, locA, bl, bu, cObj, Names, &
        eType, hs, x, pi, rc, INFO, mincw, miniw, minrw, &
        nS, nInf, sInf, ObjQP, &
        cu, lencu, iu, leniu, ru, lenru, &
        cw, lencw, iw, leniw, rw, lenrw)

  call sqEndF(cw, lencw, iw, leniw, rw, lenrw)

  ! Deallocate
  deallocate(bl)
  deallocate(bu)
  deallocate(x)
  deallocate(pi)
  deallocate(rc)
  deallocate(indA)
  deallocate(valA)
  deallocate(locA)
  deallocate(cObj)
  deallocate(Names)
  deallocate(row)
  deallocate(col)
  deallocate(val)

  deallocate(cw, iw, rw)
  deallocate(cu, iu, ru)

  call CUTEST_creport(status, CALLS, CPU)
  WRITE (iOut, 2000) probname, n, m, CALLS(1), CALLS(2), &
                     CALLS(5), CALLS(6), info, ObjQP, CPU(1), CPU(2)

  call CUTEST_cterminate(status)
  stop

910 CONTINUE
  WRITE(iOut, "(' CUTEst error, status = ', i0, ', stopping')") status
  stop

990 CONTINUE
  WRITE(iOut, "(' Allocation error, status = ', i0)") status
  stop

2000 FORMAT(/, 24('*'), ' CUTEst statistics ', 24('*') //               &
          ,' Package used            : SQOPT',     /                    &
          ,' Problem                 :  ', A10,    /                    &
          ,' # variables             =      ', I10 /                    &
          ,' # constraints           =      ', I10 /                    &
          ,' # objective functions   =        ', F8.2 /                 &
          ,' # objective gradients   =        ', F8.2 /                 &
          ,' # constraints functions =        ', F8.2 /                 &
          ,' # constraints gradients =        ', F8.2 /                 &
          ,' Exit code               =      ', I10 /                    &
          ,' Final f                 = ', E15.7 /                       &
          ,' Set up time             =      ', 0P, F10.2, ' seconds' /  &
          ,' Solve time              =      ', 0P, F10.2, ' seconds' // &
           66('*') /)

contains

  subroutine cutestHx &
       (nnH, x, Hx, state, cu, lencu, iu, leniu, ru, lenru)
    integer, intent(in) :: nnH, state, lencu, leniu, lenru
    integer,          intent(inout) :: iu(leniu)
    double precision, intent(inout) :: ru(lenru)
    character(len=8), intent(inout) :: cu(lencu)
    double precision, intent(in)    :: x(nnH)
    double precision, intent(out)   :: Hx(nnH)
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    logical :: gotH
    integer :: status
    double precision :: y(m)

    y(1:m) = 0.0d+0
    if (state == 1) then ! First call
       gotH = .false.
    else
       gotH = .true.
    end if

    call cutest_chprod(status, n, m, gotH, x, y, x, Hx)

  end subroutine cutestHx

  !=============================================================================

  subroutine crd2spr(m, n, ne, row, col, val, indA, locA, valA)
    integer, intent(in) :: m, n, ne, row(ne), col(ne)
    double precision, intent(in) :: val(ne)

    integer, intent(out) :: indA(ne), locA(n+1)
    double precision, intent(out) :: valA(ne)

    !===========================================================================
    ! Converts coordinate form (val, irow, jcol) into sparse-by-column format
    ! in val, ind, loc.
    !
    ! 15 Feb 2010: First version.
    !===========================================================================
    integer :: i, j, k, ii

    ! First count up elements in each column and store in loc.
    locA(1:n+1) = 0

    do k = 1, ne
       j = col(k)
       i = row(k)
       locA(j) = locA(j) + 1
    end do
    locA(n+1) = ne + 1

    ! Set the column pointers.
    do k = n, 1, -1
       locA(k) = locA(k+1) - locA(k)
    end do

    ! Put the row indices and values in and let loc(j) track the position of the
    ! elements for the jth column.
    indA(1:ne) = 0
    valA(1:ne) = zero
    do k = 1, ne
       j  = col(k)
       i  = row(k)
       ii = locA(j)

       valA(ii) = val(k)
       indA(ii) = i
       locA(j)  = ii + 1
    end do

    ! Reset the column pointers.
    locA(2:n) = locA(1:n-1)
    locA(1)   = 1

  end subroutine crd2spr

end program sqopt_main
