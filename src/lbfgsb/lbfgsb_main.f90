!     ( Last modified on 2 Jan 2013 at 13:40:00 )
program          lbfgsb_main
  !
  !  LBFGSB test driver for problems derived from SIF files.
  !
  !  Nick Gould, for CGT Productions.
  !  September 2004
  !  Revised for CUTEst,     January   2013
  !  Updated by Philip Gill, September 2017

  implicit none
  integer          :: i, n, m, maxit, status, nFuns, nFunGs
  integer          :: iflag, iprint, itns, lwa, liwa, nact, isave(44)
  integer          :: io_buffer = 11

  double precision :: epsmch, f0, f, gnorm
  double precision :: dsave(29)
  double precision :: factr, pgtol, time1, time2, wtime1, wtime2, tlimit
  double precision :: cpu(2), calls(4)

  character(len=10):: pname, spcdat
  character(len=60):: task, csave

  logical          :: lsave(4), NoSolution, OptCon(3)

  integer,          allocatable, dimension(:) :: nbd, iwa
  double precision, allocatable, dimension(:) :: x, xl, xu, g, wa
  character(10),    allocatable, dimension(:) :: xnames

  integer,          parameter :: outfile = 9, stdout = 6, input = 55, inspec = 56
  double precision, parameter :: zero   = 0.0d0, one    = 1.0d0, ten = 10.0d+0
  double precision, parameter :: infty  = 1.0d+19
  double precision, parameter :: xfactr = 0.0d0, xpgtol = 0.0d0

  !  Open the Spec file for the method.

  spcdat = 'LBFGSB.SPC'
  open( inspec, file = spcdat, form = 'formatted', status = 'old' )
  rewind inspec

  !Read input Spec data.
  !
  !  M        : the maximum number of variable metric corrections
  !  MAXIT    : the maximum number of iterations,
  !  IPRINT   : print level (<0,none,=0,one line/iteration,>1,more detail)
  !  FACTR    : the function accuracy tolerance
  !  PGTOL    : the required norm of the projected gradient
  !  TLIMIT   : CPU limit
  !
  ! The iteration will stop when
  !
  !  (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch   AND
  !                    |proj g^k|/(1+|f^k|) <  pgtol.
  ! where epsmch is the machine precision.
  ! Typical values for factr: 1.d+12 for low accuracy;
  !                           1.d+7  for moderate accuracy;
  !                           1.d+1  for extremely high accuracy.

  read(inspec, 1000) m, maxit, iprint, factr, pgtol, tlimit

  !  Close input file.

  close( inspec )

  !  Open the relevant file.

  open( input, file = 'OUTSDIF.d', form = 'formatted', status = 'old' )
  !
  !  Check to see if there is sufficient room
  !
  call CUTEST_udimen( status, input, n )
  if (status /= 0) go to 910

  liwa = 3 * n
  !     L-BFGS-B Version 2.0
  !     lwa  = 2 * m * n + 4 * n + 12 * m *m + 12 * m
  !     L-BFGS-B Version 3.0
  lwa  = 2 * m * n + 5 * n + 11 * m *m + 8 * m

  allocate( nbd(n), iwa(liwa), x(n), xl(n), xu(n), g(n), wa(lwa), xnames(n), &
            stat = status )
  if (status /= 0) go to 990

  !  Set up SIF data.

  call CUTEST_usetup( status, input, stdout, io_buffer, n, x, xl, xu )
  if (status /= 0) go to 910

  !  Set bound constraint status

  do i = 1, n
     if (xl(i) <= - infty) then
        if (xu(i) >= infty) then
           nbd(i) = 0
        else
           nbd(i) = 3
        end if
     else
        if (xu(i) >= infty) then
           nbd(i) = 1
        else
           nbd(i) = 2
        end if
     end if
  end do

  !  Obtain variable names.

  call CUTEST_unames( status, n, pname, xnames )
  if (status /= 0) go to 910

  ! Open output file (if not stdout)
  if ( outfile /= 6 ) then
     open(outfile, file=trim(pname)//'.out', status='new')
  end if

  !  Set up algorithmic input data.

  iflag = 0

  !  Optimization loop

  task  = 'START'

  ! Start the CPU time.

  call  timer(time1)
  call wtimer(wtime1)

  !=========================================================================

  do while (task(1:2) == 'FG' .or. task == 'NEW_X' .or. task == 'START')

     optCon(1:3) = .false.

     !  Call the optimizer

     call setulb( n, m, x, xl, xu, nbd, f, g, xfactr, xpgtol, wa, &
                  iwa, task, iprint, csave, lsave, isave, dsave )

     if (task(1:2) == 'FG') then

        ! The minimization routine  requests the
        ! function f and gradient g values at the current x.
        ! Before evaluating f and g we check the CPU time spent.

        call timer(time2)

        if (time2-time1 > tlimit) then

           task = 'STOP: CPU LIMIT EXCEEDED.'

           ! Note: Assigning task(1:4)='STOP' will terminate the run;
           ! setting task(7:9)='CPU' will restore the information at
           ! the latest iterate generated by the code so that it can
           ! be correctly printed by the driver.

           iflag = 4
           write( outfile, "( ' Time limit ' )" )
        else

           ! The time limit has not been reached and we compute
           ! the function value f.

           call CUTEST_uofg( status, n, x, f, g, .true. )
           if (status /= 0 ) go to 910
        endif

        ! Return to the minimization routine.
     else

        f0 = dsave(2); gnorm = dsave(13); epsmch = dsave(5)

        if (task(1:4) == 'CONV') then
           optCon(1) = gnorm  <= pgtol*(one + abs(f))
           optCon(2) = f0 - f <  factr*epsmch*max(abs(f0),abs(f),one)
           optCon(3) = gnorm  <= sqrt(epsmch)
           iflag = 0 ! not going to happen because xfactr, xpgtol=0

        else if (task(1:4) == 'ABNO' ) then
           iflag = 1
           write( outfile, "( ' Abnormal exit ' )" )
        else if (task(1:5) == 'ERROR' ) then
           iflag = 2
           write( outfile, "( ' Error exit ' )" )

        else if (task(1:5) == 'NEW_X') then

           ! The minimization routine has returned with a new iterate.
           ! The time limit has not been reached.
           ! Test the termination conditions.

           if (isave(30) > maxit) then
              iflag = 3
              write( outfile, "( ' Maximum number of iterations exceeded ' )" )
              task = 'STOP: total no. of iterations exceeds limit'
           end if

           optCon(1) = gnorm  <= pgtol*(one + abs(f))
           optCon(2) = f0 - f <  factr*epsmch*max(abs(f0),abs(f),one)
           optCon(3) = gnorm  <= sqrt(epsmch)
           if (OptCon(1) .and. OptCon(2)  .or.  OptCon(3)) then
              task  = 'STOP: the termination conditions are satisfied'
              iflag = 0
           end if
        endif

!!$           write(6,*) ' gnorm:', gnorm
!!$           write(6,*) ' f0 f:', abs(f0), abs(f)
!!$           write(6,*) ' tol:', pgtol, factr, sqrt(epsmch)
!!$           write(6,*) f0 - f <  factr*epsmch*max(abs(f0),abs(f),one)
!!$           write(6,*) gnorm  <= pgtol*(one + abs(f))
!!$           write(6,*) gnorm  <= sqrt(epsmch)
     end if
  end do

  gnorm = dsave( 13 )

  if (iflag ==0  .and. .not. ((OptCon(1) .and. OptCon(2)) .or. OptCon(3))) iflag = 5
 !if (iflag /=0  .and.                 gNorm <=  pgTol*(one + abs(f))*ten) iflag = 6

  ! Terminal exit.
  call wtimer(wtime2)
  call CUTEST_ureport( status, calls, cpu )
  if (status /= 0) go to 910

  NoSolution = .false.
  if (NoSolution) then
     write( outfile, 2010 ) f, gnorm
  else
     write( outfile, 2015 ) f, gnorm
     do  i = 1, n
        write( outfile, 2020 ) xnames(i), xl(i), x(i), xu(i), g(i)
     end do
  end if

  write(outfile, 2000) pname, n, int(calls(1)), int(calls(2)), iflag, f, cpu(1), cpu(2)

  ! Print results to a file
  ! Added by PEG. Aug 28 2017.

  nFuns  = calls(1)
  nFunGs = calls(2)
  itns   = isave(30)
  nact   = isave(39)

  call lbGetStats &
       (pname(1:10), n, iflag, itns, nact, nFuns, nFunGs, f, g, gnorm, &
        x, xl, xu, cpu(2), wtime2 - wtime1)

  close(input)
  call CUTEST_uterminate(status)

  stop

910 continue
  write( stdout, "( ' CUTEst error, status = ', i0, ', stopping' )") status
  stop

990 continue
  write( stdout, "( ' Allocation error, status = ', I0 )" ) status
  stop
!
!  Non-executable statements.
!
1000 format( 3( I10, / ), 2( D10.3, / ), D10.3 )
2000 format( /, 24('*'), ' CUTEst statistics ', 24('*') //, &
          ' Package used            :  L-BFGS-B',     /,   &
          ' Problem                 :  ', A10,    /,       &
          ' # variables             =      ', I10 /,       &
          ' # objective functions   =      ', I10 /,       &
          ' # objective gradients   =      ', I10 /,       &
          ' Exit code               =      ', I10 /,       &
          ' Final f                 = ', E15.7 /,          &
          ' Set up time             =      ', 0P, F10.2, ' seconds' /  &
          ' Solve time              =      ', 0P, F10.2, ' seconds' // &
           66('*') / )
2010 format( ' Final objective function value   = ', 1P, D12.4, &
             /, ' Final norm of projected gradient = ', 1P, D12.4 )
2015 format( ' Final objective function value   = ', 1P, D12.4,     &
              /, ' Final norm of projected gradient = ', 1P, D12.4, &
              //, '                XL           X        XU',       &
                 '           G ' )
2020 format(1x, a10, 1p, 4d12.4)

end program lbfgsb_main

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine lbGetStats &
     (probName, n, INFO, itns, nact, nFuns, nFunGs, &
      f, g, gNorm, x, bl, bu, cpuTime, wallTime)

  implicit         none
  integer          :: INFO, itns, n, nact, nFuns, nFunGs
  double precision :: f, gNorm, cpuTime, wallTime
  double precision :: g(n), x(n), bl(n), bu(n)
  character        :: ProbName*10

  !==================================================================
  ! lbGetStats prints statistics associated with an LBFGSB run.
  !
  !  CALLS(1)  # objective functions
  !  CALLS(2)  # objective gradients
  !
  !  INFO   0  Optimal Solution Found.
  !  INFO   1  Abnormal exit
  !  INFO   2  Error exit
  !  INFO   3  Maximum number of iterations exceeded
  !
  !==================================================================
  character                   :: flag*3, msg*43,  String*150
  logical                     :: FileExists
  integer                     :: j, nDegen, nFixed
  double precision            :: rnFuns
  double precision, parameter :: zero = 0.0d+0

  integer,          parameter :: iAll = 7, iTeX = 8, iStats = 9, iPP = 11

  itns     = mod( itns  , 1000000 )
  nFuns    = mod( nFuns , 1000000 )

  nDegen   = 0
  nFixed   = 0
  do j = 1, n
     if (x(j) == bl(j) .or. x(j) == bu(j)) then
        if (abs(g(j)) <= 2.2d-15) nDegen = nDegen + 1
        nFixed = nFixed + 1
     end if
  end do

  ! INFO   0  Optimal Solution Found.
  ! INFO   1  Abnormal exit
  ! INFO   2  Error exit
  ! INFO   3  Maximum number of iterations exceeded

  flag = ' '

  if (     INFO  ==  0) then
     flag   = '   '
     msg    = 'Optimal solution found'
  else if (INFO  ==  1) then
     flag   = 'abn'
     msg    = 'Abnormal exit'
  else if (INFO  ==  2) then
     flag   = 'Err'
     msg    = 'Error exit'
  else if (INFO  ==  3) then
     flag   = 'Itr'
     msg    = 'Iteration limit exceeded'
  else if (INFO  ==  4) then
     flag   = 'Cpu'
     msg    = 'Cpu time limit exceeded'
  else if (INFO  ==  5) then
     flag  = 'Cbi'
     msg   = 'Current point cannot be improved'
  else if (INFO  ==  6) then
     flag  = 'acc'
     msg   = 'Optimal, but accuracy not achieved'
     !        1234567890123456789012345678901234
  else
     flag   = '???'
     msg    = ' '
  end if

  !------------------------------------------------------------------
  ! Write the summary of the problem statistics
  !------------------------------------------------------------------
  inquire( FILE='temp.stats', EXIST=FileExists )
  if (FileExists) then
     open (iStats,file='temp.stats',status='old',position='append')
  else
     open (iStats,  file='temp.stats',  status='new')
     write(iStats, 2000)
  end if

  write(String, 2100) probName, n, nFixed, nDegen, gNorm, f
  write(iStats, 7000) trim(String)

  close(iStats)

  !------------------------------------------------------------------
  ! Tex summary
  !------------------------------------------------------------------
  inquire( FILE='temp.tex', EXIST=FileExists )
  if (FileExists) then
     open (iTeX,  file='temp.tex', status='old', position='append')
  else
     open (iTeX,  file='temp.tex', status='new')
     write(iTeX, 3000)
  end if

  write(String, 3100) &
       probName, itns, nFuns, gNorm, f, wallTime, cpuTime, flag
  write(iTeX,7000) trim(String)

  close(iTeX)

  !------------------------------------------------------------------
  ! All data
  !------------------------------------------------------------------
  inquire( FILE='temp.all', EXIST=FileExists )
  if (FileExists) then
     open (iAll, file='temp.all', status='old', position='append')
  else
     open (iAll, file='temp.all', status='new')
     write(iAll, 4000)
  end if

  write(String, 4100) probName, n, itns, nact, nFuns, wallTime, cpuTime, gNorm, f, msg
  write(iAll,7000) trim(String)

  close(iAll)

  !------------------------------------------------------------------
  ! Performance profile
  !------------------------------------------------------------------
  inquire( FILE='temp.pp', EXIST=FileExists )
  if (FileExists) then
     open (iPP, file='temp.pp', status='old', position='append')
  else
     open (iPP, file='temp.pp', status='new')
     write(iPP, 6000)
  end if

  wallTime = max(1.0d-4, wallTime)
  cpuTime  = max(1.0d-4, cpuTime)
  rnFuns   = nFuns
  rnFuns   = max(1.0d-3,  rnFuns )
  write(String, 6100) probName, INFO, wallTime, cpuTime, rnFuns
  if (INFO == 0) then
        !        Relax
  else
     String(17:59) = ' '
     String(28:30) = 'NaN'
     String(42:44) = 'NaN'
     String(57:59) = 'NaN'
  end if

  write(iPP,7000) trim(String)

  close (iPP)

2000 format( 'Problem Name', 3x, 'Variables',              &
          5x, 'Fixed', 2x, 'Dual dgn', 2x, 'Norm g(free)', &
          4x, 'Final objective' )
2100 format( a10, 4x, 2i10, i10, 5x, e9.2, 1x, 1p, e18.6)

3000 format(' \\begin{tabular}{lrrrrrrl}' /                &
            ' Problem Name&      Itns',                    &
            ' &       Fun &', '    norm(g)    ',           &
            ' & Objective  &    Time& CPUTime& Status \\\\[2ex]')
3100 format(a12, 2(' &',i10), ' &$',                       &
            1p, e13.6, ' $&$', e9.2 , ' $&', 0p, f8.2, '&', f8.2, '&{\em ',&
            a3, '} \\' )

!3200 format(' \\end{tabular}'
!    &    // ' \\begin{verbatim}')

4000 format(6x, 'Name', 6x, 'n', 6x, 'Itns', 6x, 'nAct', 6x, &
            'nFun', 1x, 3x, 'wall(s)', 1x, 3x, 'cpu(s)', 2x, &
            'norm gfree', 6x, 'Obj', &
            11x, 'Result')
4100 format( a10, i7, 3i10, 2x, f9.2, 1x, f9.2, 3x, 1p, e9.2, 2x, e16.9, 2x, a)

6000 format('  ProbName', 2x, 'Info', 10x, 'Time', 10x, ' CPU', 11x, 'Funs' )
6100 format(a10, 4x, i2, 2x, f12.4, 2x, f12.4, 3x, f12.3)
7000 format( a )

end  subroutine lbGetStats

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine wtimer(wtime)
  use iso_fortran_env, only : real64
  implicit none
  ! Wall timer

  integer :: count, rate
  double precision :: wtime

  call system_clock(count, rate)
  wtime = real(count,kind=real64) / real(rate,kind=real64)

end subroutine wtimer

subroutine timer( ttime )
  implicit none
  !  CPU timer

  double precision :: ttime

  call cpu_time( ttime )

end subroutine timer
