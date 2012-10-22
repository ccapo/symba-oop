program follow_symba
!--------------------------------------------------------------------------
!			FOLLOW_SYMBA.F90
!--------------------------------------------------------------------------
!
! Extracts information on a particle in a binary file, and writes it to
! an ascii file
!
!--------------------------------------------------------------------------
use swift
use util
use io
implicit none

type(param_t) :: param

logical(lk) :: lfound
integer(ik) :: nbod, narg, command_argument_count, ierr
integer(ik) :: i, id, ifol, ntp, iu = 20

real(rk) :: t, tmax, mstar
real(rk) :: mu, rad, a, e, inc, capom, omega, capm, peri, apo, obar

character(len = 24) :: paramfile, plfile
character(len = 24) :: buffer, followfile

! Get the number of command line arguments
narg = command_argument_count()

! Read in the command line argument
if(narg == 3) then

  call get_command_argument(1, buffer)
  read(buffer,*) paramfile
  write(*,'(2a)')   " Parameter data filename = ", trim(adjustl(paramfile))

  call get_command_argument(2, buffer)
  read(buffer,*) plfile
  write(*,'(2a)')   " Particle data filename = ", trim(adjustl(plfile))

  call get_command_argument(3, buffer)
  read(buffer,*) ifol
  ifol = abs(ifol)
  write(*,'(a,i6)') " Following particle number = ", ifol

else

  write(*, '(a)', advance = 'no') " Input parameter data filename: "
  read(*,*) paramfile

  write(*, '(a)', advance = 'no') " Input particle data filename: "
  read(*,*) plfile

  write(*,'(a)', advance = 'no')  " Input the particle number to follow: "
  read(*,*) ifol
  ifol = abs(ifol)

end if

! Get data for the run and the test particles
call io_init_param(paramfile, param)

! Get planet data
open(unit = 7, file = plfile, status = 'OLD', form = 'FORMATTED', iostat = ierr)

! Read number of bodies
read(7,'(1x,i7,1x,1pe22.15)') nbod, mstar

! Close file
close(unit = 7)

! Print binary format and open binary file
select case(param.outtype)

  case("FXDR4")

    !write(*,'(a,/)') ' Reading a single precision FXDR binary file'
    call io_open_fxdr(param.outfile, "R", .TRUE., iu, ierr)

  case("FXDR8")

    !write(*,'(a,/)') ' Reading a double precision FXDR binary file'
    call io_open_fxdr(param.outfile, "R", .TRUE., iu, ierr)

  case("REAL4")

    !write(*,'(a,/)') ' Reading a single precision binary file'
    call io_open(iu, param.outfile, "OLD", "UNFORMATTED", ierr)

  case("REAL8")

    !write(*,'(a,/)') ' Reading a double precision binary file'
    call io_open(iu, param.outfile, "OLD", "UNFORMATTED", ierr)

  case default

    write(*,'(a,/)') ' ERROR: no binary file format specified in param file'
    call util_exit(FAILURE)

end select

write(followfile,'(a7,i1,a4)') 'follow_', ifol, '.dat'
open(unit = 7, file = followfile, status = 'UNKNOWN')

write(*,'(/,a)') ' Output format:'
write(*,'(a)')   ' 1  2   3  4  5    6      7      8     9     10  11  12  '
write(*,'(a,/)') ' t, id, a, e, inc, capom, omega, capm, peri, apo, M, obar'

tmax = param.t0

do

  ierr = io_read_hdr(iu, t, nbod, ntp, param.ioutform, param.outtype)
  if(ierr /= 0) exit

  lfound = .false.

  do i = 2, nbod

    if(param.lgas) then

      ierr = io_read_line(iu, id, a, e, inc, capom, omega, capm, param.outtype, mass = mu, radius = rad)

    else

      ierr = io_read_line(iu, id, a, e, inc, capom, omega, capm, param.outtype, mass = mu)

    end if

    if(ierr /= 0) exit

    ! Check if the current body is the one we are searching
    if(abs(id) == ifol) then

      lfound = .true.
      tmax = t

      inc = inc*DEGRAD
      obar = capom + omega
      if(obar >= TWOPI) obar = obar - TWOPI
      obar = obar*DEGRAD
      capom = capom*DEGRAD
      omega = omega*DEGRAD
      capm = capm*DEGRAD
      peri = a*(1.0_rk - e)
      apo = a*(1.0_rk + e)
      write(7,1000) t, ifol, a, e, inc, capom, omega, capm, peri, apo, mu, obar

    end if

  end do

  ! Did not find particle this time step
  if(.not. lfound) exit

end do

 1000 format((1x,e15.7,1x,i6,1x,f10.4,1x,f8.5,4(1x,f9.4),2(1x,f10.4),1e13.5,1x,f10.4))

write(*,'(a,1pe11.5)') ' Tmax = ', tmax

end program follow_symba
