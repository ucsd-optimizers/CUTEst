!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Print problem stats
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program stats_main
  use CUTEst_interface_double
  implicit none

  logical :: &
       qpobj, lpobj, lincon, uncon, bdcon, fileExists
  character :: &
       code*2, pname*10
  integer   :: &
       status, m, n, nnH, nnObj, nnJac, nEQ, nLC, nnzh, nnzj, j, nfree
  double precision, allocatable :: &
       x(:), bl(:), bu(:), pi(:)
  logical, allocatable :: &
       equation(:), linear(:)
  integer, parameter :: &
       iCutest = 55, iPP = 8, iOut = 6
  double precision, parameter :: &
       infbnd = 1.0d+20

  open  ( iCutest, file = 'OUTSDIF.d', form = 'FORMATTED', status = 'old' )
  rewind( iCutest )

  m     = 0
  nnObj = 0
  nnJac = 0
  nEQ   = 0
  nLC   = 0

  nnH   = 0
  nnzh  = 0
  nnzj  = 0

  code  = '  '

  call CUTEST_cdimen( status, iCutest, n, m )
  allocate (x(n), bl(n+m), bu(n+m), equation(m), linear(m), pi(m))

  if ( m > 0 ) then
     call CUTEST_csetup( status, iCutest, iOut, iOut,     &
                         n, m, x(1:n), bl(1:n), bu(1:n),  &
                         pi, bl(n+1:n+m), bu(n+1:n+m),    &
                         equation, linear, 0, 2, 1 )
     call CUTEST_probname( status, pName )

     call CUTEST_cstats( status, nnObj, nnJac, nEq, nLC )

     call CUTEST_cdimsh( status, nnzh )
     call CUTEST_cdimsj( status, nnzj )
  else
     call CUTEST_usetup( status, iCutest, iOut, iOut, n, x, bl(1:n), bu(1:n) )
     call CUTEST_probname( status, pName )
     call CUTEST_udimsh( status, nnzh )
  end if

  nnH = max(nnJac, nnObj)
  if (m == 0) then
     nfree = 0
     do j = 1, n
        if (bl(j) < -infBnd .and. bu(j) > infBnd) then
           nfree = nfree + 1
        end if
     end do
  end if

  !qpobj  =
  lpobj  = nnH == 0
  lincon = nLC == m
  uncon  = m == 0 .and. nfree == n
  bdcon  = m == 0 .and. nfree <  n

!!$  if (qobj) then
!!$     code(1:1) = 'Q'
  if (lpobj) then
     code(1:1) = 'L'
  else
     code(1:1) = 'N'
  end if

  if (lincon) then
     code(2:2) = 'L'
  else if (uncon) then
     code(2:2) = 'U'
  else if (bdcon) then
     code(2:2) = 'B'
  else
     code(2:2) = 'C'
  end if

  if (allocated( bl )      ) deallocate( bl  )
  if (allocated( bu )      ) deallocate( bu  )
  if (allocated( x  )      ) deallocate( x   )
  if (allocated( pi )      ) deallocate( pi  )

  if (allocated( equation )) deallocate( equation )
  if (allocated( linear )  ) deallocate( linear   )

  !-----------------------------------------------------------------------------
  ! Print stats to file.
  !-----------------------------------------------------------------------------
  inquire( FILE='cutest.stats', EXIST=FileExists )
  if (FileExists) then
     open (iPP, file='cutest.stats', status='old', position='append' )
  else
     open (iPP, file='cutest.stats', status='new' )
     write(iPP,1000) &
          'Name      ', &
          '         m', &
          '         n', &
          '     nnObj', &
          '     nnJac', &
          '    neqCon', &
          '    linCon', &
          '    nlnCon', &
          '      nnzh', &
          '      nnzj'
  end if

  write(iPP,2000) pname, m, n, nnObj, nnJac, nEQ, nLC, m-nLC, nnzh, nnzj, code

  close(iCutest)
  close(iPP)

  if ( m > 0 ) then
     call CUTEST_cterminate(status)
  else
     call CUTEST_uterminate(status)
  end if

  stop

1000 format(10a10)
2000 format(a10, 9(2x, i8), 1x, a2)

end program stats_main
