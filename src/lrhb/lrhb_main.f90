!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CUTEst interface for LRHB
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program lrhb_main
  use CUTEst_interface_double
  use lrhb_module
  implicit none

  character     :: pName*10
  integer(ip)   :: j, n, Errors, INFO, status, alloc_stat

  real(rp)      :: objAdd, fobj, gNorm, cpu( 2 ), calls( 7 )

  real(rp), external :: dnrm2

  character(10), allocatable :: vnames(:)
  real(rp),      allocatable :: bl(:), bu(:), x(:), g(:)

  integer(ip),   parameter   :: iCutest = 55,  iOut   = 6, io_buffer = 11
  real(rp),      parameter   :: zero    = 0.0, infBnd = 1.0d+20

  type(lrhb_options) :: options

  !-----------------------
  ! Open problem data file
  !-----------------------
  open  ( iCutest, file = 'OUTSDIF.d', form = 'FORMATTED', status = 'old' )
  rewind( iCutest )

  !----------------------------------
  ! Get dimensions and allocate space
  !----------------------------------
  call CUTEST_udimen( status, iCutest, n )
  if ( status /= 0 ) go to 910


  !-------------------------
  ! Start setting up problem
  !-------------------------
  allocate( bl(n), bu(n), x(n), g(n), vnames(n), stat=alloc_stat )
  if (alloc_stat /= 0) GO TO 990


  call CUTEST_usetup( status, iCutest, iOut, io_buffer, n, x, bl(1:n), bu(1:n) )
  if (status /= 0) go to 910

  call CUTEST_unames( status, n, pName, vnames )

  !-----------------------------------------------------------------------------
  ! Ok, we're done with CUTEst stuff.
  !-----------------------------------------------------------------------------
  objAdd = zero

  call lrhb_initialize( trim(pName) // '.out', 'screen', options )
  call lrhb_specs( options, 'LRHB.SPC', info )

  call lrhb &
     ( pName(1:10), n, objAdd, x, bl, bu, lrhb_usrfun, fobj, g, options, info )

  call lrhb_end( options )

  call CUTEST_ureport( status, calls, cpu )
  write(iOut, 2010) pName, n, int(calls(1)), int(calls(2)), info, fobj, cpu( 1 ), cpu( 2 )

  if (allocated(bl))     deallocate(bl)
  if (allocated(bu))     deallocate(bu)
  if (allocated(x) )     deallocate(x )
  if (allocated(g) )     deallocate(g )
  if (allocated(vnames)) deallocate(vnames)

  close(iCutest)

  call CUTEST_uterminate(status)

  stop

910 continue
  write( iOut, "( ' CUTEst error, status = ', i0, ', stopping' )") status
  stop

990 continue
  write( iOut, "( ' Allocation error, status = ', i0 )" ) status
  stop

1000 format(a10, 2x, i7, 2x, i7, 2es12.2, 2x, i7, 2x, i4)
2010 format( /, 24('*'), ' CUTEst statistics ', 24('*') //                    &
         ,' Package used            : LRHB',     /                            &
         ,' Problem                 :  ', A10,    /                           &
         ,' # variables             =      ', I10 /                           &
         ,' # objective functions   =      ', I10 /                           &
         ,' # objective gradients   =      ', I10 /                           &
         ,' Exit code               =      ', I10 /                           &
         ,' Final f                 = ', E15.7 /                              &
         ,' Set up time             =      ', 0P, F10.2, ' seconds' /         &
         ,' Solve time              =      ', 0P, F10.2, ' seconds' //        &
         66('*') / )
3000 format(a8, 3x, i2, 2x, f12.3, 3x, i12)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lrhb_usrfun( mode, n, x, f, g, state )
    use  CUTEst_interface_double
    use lrhb_module, only : ip, rp
    implicit none
    integer(ip),  intent(in)    :: n, state
    integer(ip),  intent(inout) :: mode
    real(rp),     intent(in)    :: x(n)
    real(rp),     intent(inout) :: f, g(n)
    !===========================================================================
    logical :: needG

    if (mode == 0) then
       needG = .false.
    else
       needG = .true.
    end if

    call cutest_uofg(status, n, x, f, g, needG)

  end subroutine lrhb_usrfun

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end program lrhb_main
