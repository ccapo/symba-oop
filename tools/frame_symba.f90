program frame_symba
!----------------------------------------------------------------------------
!				FRAME_SYMBA.F90
!----------------------------------------------------------------------------
!
! Extracts information on all particles in a binary file, for certain times,
! and writes to separate ascii files
!
!----------------------------------------------------------------------------
use swift
use orbel
use io
use util
implicit none

! Global parameters
type(param_t) :: param

! Particle type
type(body_t), allocatable :: pbod(:)

integer(ik) :: narg, command_argument_count, system, status
integer(ik) :: nframe, nplot, nf
integer(ik) :: nbod, nbodm, ntp, ialpha
integer(ik) :: i, j, id, ierr, iu = 20

real(rk) :: gm, mu, rad, a, e, inc, capom, omega, capm
real(rk) :: t, mscale, tframe, dtframe, dtdelay, amin, amax, mpl, fg, fs, logr, rplsml
character(len = 24) :: paramfile, plfile
character(len = 80) :: buffer, cmd, fn, tf, aminf, amaxf, mf, gf, sf, rf

 1000   format(i5,1x,e13.5,1x,6(f10.5,1x))

! Get the number of command line arguments
narg = command_argument_count()

! Read in the command line arguments
if(narg == 10) then

  call get_command_argument(1, buffer)
  read(buffer,*) paramfile

  call get_command_argument(2, buffer)
  read(buffer,*) plfile

  call get_command_argument(3, buffer)
  read(buffer,*) nframe

  call get_command_argument(4, buffer)
  read(buffer,*) dtdelay

  call get_command_argument(5, buffer)
  read(buffer,*) amin

  call get_command_argument(6, buffer)
  read(buffer,*) amax

  call get_command_argument(7, buffer)
  read(buffer,*) mpl

  call get_command_argument(8, buffer)
  read(buffer,*) fg

  call get_command_argument(9, buffer)
  read(buffer,*) fs

  call get_command_argument(10, buffer)
  !read(buffer,*) logr
  read(buffer,*) rplsml

else

  write(*, '(a)', advance = 'no') " Input parameter data filename: "
  read(*,*) paramfile

  write(*, '(a)', advance = 'no') " Input particle data filename: "
  read(*,*) plfile

  write(*, '(a)', advance = 'no') " Input the number of frames: "
  read(*,*) nframe

  write(*, '(a)', advance = 'no') " Input the time delay per frame (sec): "
  read(*,*) dtdelay

  write(*, '(a)', advance = 'no') " Input the Min. semi-major axis (AU): "
  read(*,*) amin

  write(*, '(a)', advance = 'no') " Input the Max. semi-major axis (AU): "
  read(*,*) amax

  write(*, '(a)', advance = 'no') " Input the mass for the embryo (M_Earth): "
  read(*,*) mpl

  write(*, '(a)', advance = 'no') " Input the gas scale factor (MMSN): "
  read(*,*) fg

  write(*, '(a)', advance = 'no') " Input the solid scale factor (MMSN): "
  read(*,*) fs

  !write(*, '(a)', advance = 'no') " Input the logarithm of the planetesimal radius (km): "
  !read(*,*) logr

  write(*, '(a)', advance = 'no') " Input the planetesimal radius (km): "
  read(*,*) rplsml

end if

! Get data for the run and the test particles
call io_init_param(paramfile, param)

! Get planet data
call io_init_pl_symba(plfile, param, nbod, nbodm, pbod)

! Useful definitions
mscale = pbod(1)%mass*3.040432646e-06_rk
tframe = param%t0
nf = -1

if(nframe <= 1) nframe = 2
dtframe = (param%tstop - param%t0)/real(nframe - 1, rk)
if(modulo(dtframe,param%dtout) /= 0) dtframe = param%dtout*(floor(dtframe/param%dtout) + 1)
write(*,'(a,i5,a,/)') 'Creating ', nframe, ' snapshots from t0 to tstop'

!rplsml = 10.0**logr
write(tf,'(1pe9.3)') param%t0
write(aminf,'(f4.1)') amin
write(amaxf,'(f4.1)') amax
write(mf,'(f5.2)') mpl
write(gf,'(f4.1)') fg
write(sf,'(f4.1)') fs
write(rf,'(1pe9.3)') rplsml

open(unit = 7, file = 'plsml.dat', status = 'unknown')
open(unit = 8, file = 'embryo.dat', status = 'unknown')

if(param%lgas) then

  do i = 2, nbod

    gm = pbod(1)%mass + pbod(i)%mass
    call orbel_xv2el(pbod(i)%r, pbod(i)%v, gm, ialpha, a, e, inc, capom, omega, capm)
    if(pbod(i)%mass < param%mtiny) then

      write(7,'(i7,8(1x,1pe22.15))') pbod(i)%id, pbod(i)%mass/mscale, a, e, inc, omega, capom, capm, pbod(i)%rdrag

    else

      write(8,'(i7,8(1x,1pe22.15))') pbod(i)%id, pbod(i)%mass/mscale, a, e, inc, omega, capom, capm, pbod(i)%rdrag

    end if

  end do

else

  do i = 2, nbod

    gm = pbod(1)%mass + pbod(i)%mass
    call orbel_xv2el(pbod(i)%r, pbod(i)%v, gm, ialpha, a, e, inc, capom, omega, capm)
    if(pbod(i)%mass < param%mtiny) then

      write(7,'(i7,7(1x,1pe22.15))') pbod(i)%id, pbod(i)%mass/mscale, a, e, inc, omega, capom, capm

    else

      write(8,'(i7,7(1x,1pe22.15))') pbod(i)%id, pbod(i)%mass/mscale, a, e, inc, omega, capom, capm

    end if

  end do

end if

close(unit = 7)
close(unit = 8)

! Deallocate unused variable
deallocate(pbod)

nf = nf + 1
if((0 <= nf) .and. (nf <= 9))       write(fn, '("000", i1)') nf
if((10 <= nf) .and. (nf <= 99))     write(fn, '("00",  i2)') nf
if((100 <= nf) .and. (nf <= 999))   write(fn, '("0",   i3)') nf
if((1000 <= nf) .and. (nf <= 9999)) write(fn, '(i4)') nf

write(tf, '(1pe9.3)') param%t0

write(cmd,*) './@plot_frame ' // trim(fn) // ' ' // trim(tf) // ' ' // trim(aminf) // ' ' // trim(amaxf) // &
             ' ' // trim(mf) // ' ' // trim(gf) // ' ' // trim(sf) // ' ' // trim(rf)

status = system(trim(cmd))

! Print binary format and open binary file
select case(param%outtype)

  case("FXDR4")

    !write(*,'(a,/)') ' Reading a single precision FXDR binary file'
    call io_open_fxdr(param%outfile, "R", .TRUE., iu, ierr)

  case("FXDR8")

    !write(*,'(a,/)') ' Reading a double precision FXDR binary file'
    call io_open_fxdr(param%outfile, "R", .TRUE., iu, ierr)

  case("REAL4")

    !write(*,'(a,/)') ' Reading a single precision binary file'
    call io_open(iu, param%outfile, "OLD", "UNFORMATTED", ierr)

  case("REAL8")

    !write(*,'(a,/)') ' Reading a double precision binary file'
    call io_open(iu, param%outfile, "OLD", "UNFORMATTED", ierr)

  case default

    write(*,'(a,/)') ' ERROR: no binary file format specified in param file'
    call util_exit(FAILURE)

end select

! Loop through entries in the binary file, locating desired frames
do while(tframe <= param%tstop)

  ! Read header of current frame of binary file
  ierr = io_read_hdr(iu, t, nbod, ntp, param%ioutform, param%outtype)
  if(ierr /= 0) then

    write(*,'(a)') 'Could not reader header of binary!'
    stop

  end if

  ! If the current frame matches the desired frame, write contents and create a plot
  if(t == tframe) then

    ! Open output files from planetesimals and embryos
    open(unit = 7, file = 'plsml.dat', status = 'unknown')
    open(unit = 8, file = 'embryo.dat', status = 'unknown')

    ! If the gas drag flag is set, read and write the drag radius in additional to other standard variables
    if(param%lgas) then

      ! Loop over all the entries in this frame
      do i = 2, nbod

        ! Read entry for body i
        ierr = io_read_line(iu, id, a, e, inc, capom, omega, capm, param%outtype, mass = mu, radius = rad)
        if(ierr /= 0) then

          write(*,'(a)') 'Could not read entry in binary!'
          stop

        end if

        ! Write entry to either output file, depending on its mass
        if(mu < param%mtiny) then

          write(7,'(i7,8(1x,1pe22.15))') id, mu/mscale, a, e, inc, omega, capom, capm, rad

        else

          write(8,'(i7,8(1x,1pe22.15))') id, mu/mscale, a, e, inc, omega, capom, capm, rad

        end if

      end do

    else

      ! Loop over all the entries in this frame
      do i = 2, nbod

        ! Read entry for body i
        ierr = io_read_line(iu, id, a, e, inc, capom, omega, capm, param%outtype, mass = mu)
        if(ierr /= 0) then

          write(*,'(a)') 'Did not find time!'
          stop

        end if

        ! Write entry to either output file, depending on its mass
        if(mu < param%mtiny) then

          write(7,'(i7,7(1x,1pe22.15))') id, mu/mscale, a, e, inc, omega, capom, capm

        else

          write(8,'(i7,7(1x,1pe22.15))') id, mu/mscale, a, e, inc, omega, capom, capm

        end if

      end do

    end if

    ! Close output files
    close(unit = 7)
    close(unit = 8)

    ! Increment output file counter, then store value as a string
    nf = nf + 1
    if((0 <= nf) .and. (nf <= 9))       write(fn, '("000", i1)') nf
    if((10 <= nf) .and. (nf <= 99))     write(fn, '("00",  i2)') nf
    if((100 <= nf) .and. (nf <= 999))   write(fn, '("0",   i3)') nf
    if((1000 <= nf) .and. (nf <= 9999)) write(fn, '(i4)') nf

    ! Store the current value of t as a string
    write(tf, '(1pe9.3)') t

    ! Store command as a string
    write(cmd,*) './@plot_frame ' // trim(fn) // ' ' // trim(tf) // ' ' // trim(aminf) // ' ' // trim(amaxf) // &
                 ' ' // trim(mf) // ' ' // trim(gf) // ' ' // trim(sf) // ' ' // trim(rf)

    ! Send command to command line
    status = system(trim(cmd))

    ! Increment tframe
    tframe = tframe + dtframe

  else

    ! If the current frame does not matche the desired frame, go over entries to the next frame
    if(param%lgas) then

      ! Loop over all the entries in this frame
      do i = 2, nbod

        ! Read entry for body i
        ierr = io_read_line(iu, id, a, e, inc, capom, omega, capm, param%outtype, mass = mu, radius = rad)

      end do

    else

      ! Loop over all the entries in this frame
      do i = 2, nbod

        ! Read entry for body i
        ierr = io_read_line(iu, id, a, e, inc, capom, omega, capm, param%outtype, mass = mu)

      end do

    end if

  end if

end do

if(tframe >= param%tstop) then

  write(*,'(a)') '*** Make MPEG ***'
  status = system("rm plsml.dat embryo.dat")
  !status = system("convert *.png movie.mpg >& /dev/null")
  write(cmd,'(a,f3.1,a)') "ffmpeg -r ", dtdelay, " -i 'frame_%04d.png' -y -r 25 'movie.mpg' >/dev/null 2>&1"
  !write(*,*) trim(cmd)
  status = system(trim(cmd))
  status = system("rm *.png")

end if

write(*,'(a)') '*** Done ***'

end program frame_symba
