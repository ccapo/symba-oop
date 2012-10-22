program swift_symba5
!----------------------------------------------------------------------
!			SWIFT_SYMBA5.F90
!----------------------------------------------------------------------
!
! To run one needs two input files:
!
! a parameter file: param.in
! a particle file: pl.in
!
! NOTE: No test particles in this code
!
! Created by Hal and Martin September, 2006
!
! Authors: Hal Levison & Martin Duncan
! Date: 11/21/96
! Last revision: 12/27/96
!              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
!                             - A number of other modifications, see
!                               individual routines for details
!----------------------------------------------------------------------
use swift
use util
use io
use discard
use symba5
implicit none

! Command line variables
integer(ik) :: narg, command_argument_count
character(LEN = 80) :: buffer

! Global parameters
type(param_t) :: param

! Number of total bodies and massive bodies
integer(ik) :: nbod, nbodm, nbod_old

! Particle type
type(body_t), allocatable :: pbod(:)

! Merger list
integer(ik), pointer :: pmerge(:)

! Current time, output time, dump time and time fraction
real(rk) :: t, tout, tdump, tfraction

! Energy offset
real(rk) :: eoffset = 0.0_rk

! Input and dump file names
character(len = 24) :: paramfile, plfile
character(len = 24) :: dumpparamfile = "dump_param.dat", dumpplfile = "dump_pl.dat"

!-----------------!
! Executable code !
!-----------------!

! Print program header
call util_version

! Get the number of command line arguments
narg = command_argument_count()

! Read in parameter and particle data filenames from the command line
if(narg == 2) then

  write(*,'(a)')    " Input parameter and particle filenames:"
  write(*,'(a)')    " ---------------------------------------"

  call get_command_argument(1, buffer)
  read(buffer,*) paramfile
  write(*,'(2a)')   " Parameter data filename = ", trim(adjustl(paramfile))

  call get_command_argument(2, buffer)
  read(buffer,*) plfile
  write(*,'(2a,/)') " Particle data filename  = ", trim(adjustl(plfile))

else ! Otherwise assume default filenames

  write(*,'(a)')    " Input parameter and particle filenames (Default):"
  write(*,'(a)')    " -------------------------------------------------"

  paramfile = "param.in"
  write(*,'(2a)')   " Parameter data filename = ", trim(adjustl(paramfile))

  plfile = "pl.in"
  write(*,'(2a,/)') " Particle data filename  = ", trim(adjustl(plfile))

end if

! Get parameters for simulation
call io_init_param(paramfile, param)

! Read in particle data, defining the size of each vector/array
call io_init_pl_symba(plfile, param, nbod, nbodm, pbod)

! Set initial time and times for first output and dump
t = param.t0
tout = param.t0 + param.dtout
tdump = param.t0 + param.dtdump

! Do the initial write to binary file
call io_write_frame(param, t, nbod, pbod)

! If we are following the energy and angular momentum of the system,
! then compute and write initial values to energy.dat
if(param.lenergy) call io_write_energy(param, t, nbod, nbodm, pbod, eoffset)

! *** Start the Main Loop ***
write(*,'(/,a,/)') " ************** BEGIN MAIN LOOP ******************"

! Loop until we reach tstop, and while there is more than one body in the system
do while((t <= param.tstop) .and. (nbod > 1))

  ! Perform a single step
  call symba5_step(param, t, nbod, nbodm, pbod, pmerge, eoffset)

  ! Update time
  t = t + param.dt

  ! Check if any bodies need to be removed from the simulation
  if(param.lclose) then

    ! Store current number of bodies
    nbod_old = nbod

    ! Test if any bodies need to be removed
    call discard_massive5(param, t, nbod, nbodm, pbod, pmerge, eoffset)

    ! Update the number of massive bodies, and ensure they are sorted in decreasing mass
    if(nbod_old /= nbod) call util_nbodm(param, nbod, nbodm, pbod)

  end if

  ! If it is time, output position and velocity of all the bodies
  if(t >= tout) then

    ! Write information to binary file
    call io_write_frame(param, t, nbod, pbod)

    ! Update output time
    tout = tout + param.dtout

  end if

  ! If it is time, do a dump
  if(t >= tdump) then

    ! Print current status of the simulation
    tfraction = (t - param.t0)/(param.tstop - param.t0)
    write(*,'(a,1pe12.5,a,0pf6.3,a,i9)') " Time = ", t, "; Fraction complete =", tfraction, "; Number of bodies = ", nbod

    ! Dump position and velocity of all the bodies, along with the simulation parameters
    !call io_dump_pl_symba(dumpplfile, param, nbod, pbod)
    call io_dump_param(dumpparamfile, param, t)

    ! Update dump time
    tdump = tdump + param.dtdump

    ! If we are following the energy and angular momentum of the system,
    ! then compute and append to energy.dat
    if(param.lenergy) call io_write_energy(param, t, nbod, nbodm, pbod, eoffset)

  end if

end do

! *** End of the Main Loop ***
write(*,'(/,a)') " *************** END MAIN LOOP *******************"

! Do a final dump of the position and velocity of all the bodies, and the simulation parameters for possible resumption
t = param.tstop
call io_dump_pl_symba(dumpplfile, param, nbod, pbod)
call io_dump_param(dumpparamfile, param, t)

! Deallocate particle list
deallocate(pbod)

! Deallocate merger list
if(associated(pmerge)) deallocate(pmerge)

! We are done!
call util_exit(SUCCESS)

end program swift_symba5
