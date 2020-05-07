program sqic_main
  use sqic_module

  implicit none

  integer(ip), parameter :: iCutest = 55, iOut = 6, io_buffer = 11

  integer(ip)   :: status, alloc_stat, j, k, i
  real(rp)      :: cpu(2), calls(7)

  integer(ip), allocatable :: Hrow(:), Hcol(:)
  real(rp),    allocatable :: zero(:), b(:), Hval(:)
  logical,     allocatable :: equation(:), linear(:)

  ! SQIC variables
  type(sqic_problem) :: prob
  type(sqic_options) :: options
  type(sqic_stats)   :: stats

  integer(ip)   :: INFO, n, m, nm, neA, lenA, neH, lenH
  real(rp)      :: Obj, ObjAdd

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

  prob%m = m
  prob%n = n

  ! Allocate SQIC structures
  allocate &
       (prob%bl(nm), prob%bu(nm), prob%cObj(n), prob%x(nm), &
        prob%pi(m), prob%rc(nm), prob%varnames(n), prob%connames(m), &
        stat=alloc_stat)
  allocate &
       (b(m), zero(n), equation(m), linear(m), &
        stat = alloc_stat)
  if (alloc_stat /= 0) GO TO 990

  allocate(snAij :: prob%A)
  allocate(snHij :: prob%H)

  zero(:)   = 0.0
  prob%x(:) = 0.0

  ! Get initial x, bounds, constraint types
  call CUTEST_csetup &
       (status, iCutest, iOut, io_buffer, &
        n, m, prob%x(1:n), prob%bl(1:n), prob%bu(1:n), &
        prob%pi, prob%bl(n+1:nm), prob%bu(n+1:nm), &
        equation, linear, 0, 2, 1)
  if (status /= 0) go to 910

  ! Check that all constraints are linear
  ! ...

  deallocate(equation, linear)


  ! Variable/constraint names
  call CUTEST_cnames(status, n, m, probname, prob%varnames(1:n), prob%connames(1:m))
  if (status /= 0) go to 910
  prob%name(1:8) = probname(1:8)


  ! Constraint matrix
  call CUTEST_cdimsj(status, lenA)
  if (status /= 0) go to 910

  lenA = lenA - n

  select type(A => prob%A)
  type is(snAij)
     call A%init(m, n, lenA)
     call CUTEST_ccfsg &
          (status, n, m, zero, b, neA, lenA, A%val, A%col, A%row, .true.)
     if (status /= 0) go to 910

     A%nnz = neA
  end select

  ! Adjust bounds
  prob%bl(n+1:nm) = prob%bl(n+1:nm) - b
  prob%bu(n+1:nm) = prob%bu(n+1:nm) - b
  deallocate(b)

  ! Objective
  ! cofg returns ObjAdd and cObj.
  call CUTEST_cofg(status, n, zero, prob%fObj, prob%cObj, .true.)
  if (status /= 0) go to 910

  ! Hessian of the objective
  call CUTEST_cdimsh(status, lenH)
  if ( status /= 0 ) go to 910

  select type(H => prob%H)
  class is (snHij)
     if (lenH > 0) then
        call H%init(n, lenH)

        allocate(Hrow(lenH), Hcol(lenH), Hval(lenH), stat = alloc_stat)
        if (alloc_stat /= 0) GO TO 990

        call CUTEST_cish(status, n, zero, 0, neH, lenH, Hval, Hrow, Hcol)
        if (status /= 0) go to 910

        i = 0
        do k = 1, neH
           if (Hrow(k) <= Hcol(k)) then
              i = i + 1
              H%row(i)  = Hcol(k)
              H%col(i)  = Hrow(k)
              H%val(i)  = Hval(k)
           end if
        end do

        neH   = i
        H%n   = maxval(Hcol(1:neH))
        H%nnz = neH

        deallocate(HRow, Hcol, Hval)
     else
        H%n   = 0
        H%nnz = 0
     end if
  end select


  !-----------------------------------------------------------------------------
  ! Ok, we're done with CUTEst stuff.
  !-----------------------------------------------------------------------------

  filename = trim(probname)//'.out'
  call sqic_initialize(filename, 'screen', options)
  call sqic_specs(options, 'SQIC.SPC', info)
  call sqic('Cold', prob, options, INFO, stats)

  Obj = stats%objective

  call sqic_end(options)


  call CUTEST_creport(status, CALLS, CPU)
  WRITE (iOut, 2000) probname, n, m, CALLS(1), CALLS(2), &
                     CALLS(5), CALLS(6), info, Obj, CPU(1), CPU(2)
  close (iCutest)

  call CUTEST_cterminate(status)
  stop

910 CONTINUE
  WRITE(iOut, "(' CUTEst error, status = ', i0, ', stopping')") status
  stop

990 CONTINUE
  WRITE(iOut, "(' Allocation error, status = ', i0)") status
  stop

2000 FORMAT(/, 24('*'), ' CUTEst statistics ', 24('*') //               &
          ,' Package used            :  SQIC',    /                     &
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

!!$  subroutine cutestHx(nnH, x, Hx, jcol)
!!$    integer(ip), intent(in)  :: nnH, jcol
!!$    real(rp),    intent(in)  :: x(nnH)
!!$    real(rp),    intent(out) :: Hx(nnH)
!!$    !---------------------------------------------------------------------------
!!$    !---------------------------------------------------------------------------
!!$    logical :: &
!!$         gotH
!!$    integer(ip) :: &
!!$         status
!!$    real(rp) :: &
!!$         xx(nn), Hxx(nn), y(mm)
!!$
!!$    xx(:)     = zero
!!$    y(:)      = zero
!!$    xx(1:nnH) = x(1:nnH)
!!$
!!$    if (state == 0) then
!!$       gotH = .true.
!!$    else
!!$       gotH = .false.
!!$    end if
!!$
!!$    call cutest_chprod(status, nn, mm, gotH, xx, y, xx, Hxx)
!!$    Hx(1:nnH) = Hxx(1:nnH)
!!$
!!$  end subroutine cutestHx

end program sqic_main
