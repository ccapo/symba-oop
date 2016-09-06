module io
! Module for io routines
!
! Author:  Hal Levison
! Date:    2/21/94
! Last revision: 2/24/94
use swift
use fxdr
use util
use orbel
use anal
implicit none

! Common logical flag to control creating discard_mass.dat,
! since two routines can create or append to said file.
! Variable is made private, so only the routines in this
! module can modify it.
logical(lk), private, save :: lfirst_discard = .true.

contains
!!
!!
  subroutine io_open(iu, fname, fopenstat, fmt, ierr)
  !-------------------------------------------------------------------------
  !         IO_OPEN.F90
  !-------------------------------------------------------------------------
  ! Open files
  !
  ! Input:  iu        ==> Unit number
  !         fname     ==> File name
  !         fopenstat ==> Open status string
  !         fmt       ==> Format string
  !
  ! Output: ierr      ==> Status flag
  !
  ! Authors: Hal Levison
  ! Date: 03/03/94
  ! Last revision: 01/30/98
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik), intent(in)  :: iu
  integer(ik), intent(out) :: ierr
  character(len = *), intent(in) :: fname, fopenstat, fmt

  !-----------------!
  ! Executable Code !
  !-----------------!

  if(fopenstat == "APPEND") then

    open(unit = iu, file = fname, status = "OLD", position = fopenstat, form = fmt, iostat = ierr)

    if(ierr /= 0) then

      write(*,'(a)')  "SWIFT Warning:: io_open"
      write(*,'(3a)') " Could not open ", fname, " with position = append, will open as status = new"

      open(unit = iu, file = fname, status = "NEW", form = fmt, iostat = ierr)

    end if

  else

    open(unit = iu, file = fname, status = fopenstat, form = fmt, iostat = ierr)

  end if

  return
  end subroutine io_open
!!
  subroutine io_open_fxdr(fname, fopenstat, lflag, iu, ierr)
  !--------------------------------------------------------------------------
  !         IO_OPEN_FXDR.F90
  !--------------------------------------------------------------------------
  ! Open files using the FXDR routines
  !
  ! Input:  fname     ==> File name
  !         fopenstat ==> Open status character
  !         lflag     ==> Return on error flag
  !                       - lflag = .FALSE., then routines halt on I/O error
  !                       - lflag = .TRUE., otherwise
  !
  ! Output: iu        ==> Unit number
  !         ierr      ==> Status flag
  !
  ! Authors: Hal Levison
  ! Date: 11/03/99
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !--------------------------------------------------------------------------
  implicit none

  ! Input variables
  logical(lk), intent(in)  :: lflag
  integer(ik), intent(out) :: iu, ierr
  character(len = *), intent(in) :: fname
  character(len = 1), intent(in) :: fopenstat

  !-----------------!
  ! Executable Code !
  !-----------------!

  iu = initxdr(fname, fopenstat, lflag)

  if(iu > 0) then

    ierr = 0

  else

    ierr = iu

  end if

  return
  end subroutine io_open_fxdr
!!
  subroutine io_init_param(infile, param)
  !-----------------------------------------------------------------------------------
  !         IO_INIT_PARAM.F90
  !-----------------------------------------------------------------------------------
  ! Reads in the parameters for the integration, and stores in global parameter
  ! variable param
  !
  ! Input:  infile ==> Parameter file name
  !
  ! Output: param  ==> Global parameters (See swift module)
  !
  ! Authors: Martin Duncan
  ! Date: 03/02/93
  ! Last revision: 05/10/94 HFL
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Now all parameters are stored in a single variable
  !-----------------------------------------------------------------------------------
  !$ use omp_lib
  implicit none

  ! Input variable
  character(len = *) :: infile

  ! Output variable
  type(param_t), intent(out) :: param

  ! Internal variables
  logical(lk) :: lenergy, loblate, lclose, lgas
  logical(lk) :: lstop = .false., lout = .false., ldump = .false.
  integer(ik) :: i, ierr, iu = 7
  real(rk) :: t0, tstop, dt, dtout, dtdump
  real(rk) :: j2, j4, rmin, rmax, rmaxu, qmin
  real(rk) :: mtiny, rhogas0, gpower, zscale, rgi, rgf, taugas, ca, ce
  character(len = 24) :: outfile, outtype, outform, outstat

  ! Namelists
  namelist / sim_params / t0, tstop, dt, dtout, dtdump, outfile, outtype, outform, outstat, lenergy
  namelist / stellar_params / loblate, j2, j4, lclose, rmin, rmax, rmaxu, qmin
  namelist / symba_params / mtiny
  namelist / gas_params / lgas, rhogas0, gpower, zscale, rgi, rgf, taugas, ca, ce

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Open and read in parameters and store values in param
  call io_open(iu, infile, "OLD", "FORMATTED", ierr)
  read(unit = iu, nml = sim_params)

  ! Ensure integration time, along with output intervals are integer multiples of the timestep
  if(modulo(tstop,dt) /= 0) then

    tstop = dt*(floor(tstop/dt) + 1)
    lstop = .true.

  end if

  if(modulo(dtout,dt) /= 0) then

    dtout = dt*(floor(dtout/dt) + 1)
    lout = .true.

  end if

  if(modulo(dtdump,dt) /= 0) then

    dtdump = dt*(floor(dtdump/dt) + 1)
    ldump = .true.

  end if

  param%paramfile = infile
  param%t0 = t0
  param%tstop = tstop
  param%dt = dt
  param%dtout = dtout
  param%dtdump = dtdump
  param%outfile = outfile
  param%outtype = outtype
  param%outstat = outstat
  param%lenergy = lenergy

  param%lhelio = .true.
  param%lbary = .false.
  param%lcanon = .false.

  read(unit = iu, nml = stellar_params)
  param%loblate = loblate
  param%j2rp2 = j2
  param%j4rp4 = j4
  param%lclose = lclose
  param%laccrete = lclose ! Assumes that if we want to resolve close encounters, we want to permit accretion events as well
  param%rmin = rmin
  param%rmax = rmax
  param%rmaxu = rmaxu
  param%qmin = qmin

  read(unit = iu, nml = symba_params)
  param%mtiny = mtiny*mearth ! Convert from Earth Mass to Solar Mass

  read(unit = iu, nml = gas_params)
  param%lgas = lgas
  param%rhogas0 = rhogas0*msun*1.68314195e6_rk ! Convert from g/cm^3 to Solar Mass/AU^3
  param%gpower = gpower
  param%zscale = zscale
  param%rgi = rgi
  param%rgf = rgf
  param%taugas = taugas
  param%ca = ca
  param%ce = ce
  param%sigma0 = param%rhogas0*param%zscale*sqrt(PI)
  param%spower = param%gpower - 1.25_rk

  close(unit = iu)

  ! Print simulation parameters
  write(*,'(a)')             " Simulation parameters:"
  write(*,'(a)')             " ----------------------"
  write(*,'(a,1pe12.5,a)')   " t0     = ", param%t0,     " years"
  write(*,'(a,1pe12.5,a)')   " tstop  = ", param%tstop,  " years"
  write(*,'(a,1pe12.5,a)')   " dt     = ", param%dt,     " years"
  write(*,'(a,1pe12.5,a)')   " dtout  = ", param%dtout,  " years"
  write(*,'(a,1pe12.5,a,/)') " dtdump = ", param%dtdump, " years"

  if(lstop) write(*,'(a,/)') " *** tstop modified to be an integer multiple of dt ***"
  if(lout)  write(*,'(a,/)') " *** dtout modified to be an integer multiple of dt ***"
  if(ldump) write(*,'(a,/)') " *** dtdump modified to be an integer multiple of dt ***"

  ! Print close encounter information, if relevant
  if(param%lclose) then

    write(*,'(a)')             " Close encounter parameters:"
    write(*,'(a)')             " ---------------------------"
    write(*,'(a,1pe12.5,a)')   " rmin   = ", param%rmin,  " AU"
    write(*,'(a,1pe12.5,a)')   " rmax   = ", param%rmax,  " AU"
    write(*,'(a,1pe12.5,a)')   " rmamu  = ", param%rmaxu, " AU"
    write(*,'(a,1pe12.5,a,/)') " qmin   = ", param%qmin,  " AU"

  else

    param%rmin = -1.0_rk
    param%rmax = -1.0_rk
    param%rmaxu = -1.0_rk
    param%qmin = -1.0_rk

  end if

  ! Print binary output parameters
  write(*,'(a)')  " Binary output parameters:"
  write(*,'(a)')  " -------------------------"
  write(*,'(2a)') " Output filename = ", trim(adjustl(param%outfile))
  write(*,'(2a)') " Output type     = ", trim(adjustl(param%outtype))

  ! Set output format integer, and print choice
  if(outform == "EL") then

    param%ioutform = EL ! Osculating Keplerian orbital elements
    write(*,'(a,/)') " Output format   = Osculating Keplerian orbital elements"

  else

    param%ioutform = XV ! Cartesian positions and velocities
    write(*,'(a,/)') " Output format   = Cartesian positions and velocities"

  end if

  ! Print energy information, if relevant
  if(param%lenergy) then

    write(*,'(a)')   " Energy output parameters:"
    write(*,'(a)')   " -------------------------"
    write(*,'(a,/)') " Output filename = energy.dat"

  end if

  ! Print gas parameters, if relevant
  if(param%lgas) then

    write(*,'(a)')           " Gas disk parameters:"
    write(*,'(a)')           " --------------------"
    write(*,'(a,1pe11.5,a)') " Density of gas @ 1 AU        = ", rhogas0,        " g/cm^3"
    write(*,'(a,1pe11.5)')   " Power law index              = ", param%gpower
    write(*,'(a,1pe11.5,a)') " Gas scale height @ 1 AU      = ", param%zscale, " AU"
    write(*,'(a,1pe11.5,a)') " Inner edge of the gas disk   = ", param%rgi,    " AU"
    write(*,'(a,1pe11.5,a)') " Outer edge of the gas disk   = ", param%rgf,    " AU"
    write(*,'(a,1pe11.5,a)') " Gas decay time scale         = ", param%taugas, " years"
    write(*,'(a,1pe11.5)')   " Type-I drag efficiency [c_a] = ", param%ca
    write(*,'(a,1pe11.5,/)') " Type-I drag efficiency [c_e] = ", param%ce

  end if

  ! Define the maximum number of threads
  nthreads = 1                        ! In the *serial* case
  !$ nthreads = omp_get_max_threads() ! In the *parallel* case
  !$ write(*,'(a)')      " OpenMP parameters:"
  !$ write(*,'(a)')      " ------------------"
  !$ write(*,'(a,i3,/)') " Number of threads  = ", nthreads

  return
  end subroutine io_init_param
!!
  subroutine io_init_pl_symba(infile, param, nbod, nbodm, pbod)
  !-------------------------------------------------------------------------------------
  !         IO_INIT_PL_SYMBA.F90
  !-------------------------------------------------------------------------------------
  ! Reads in the data for the central body and other bodies for SyMBA routines
  !
  ! Input:  infile     ==> Input body file name
  !         param      ==> Global parameters (See swift module & io_init_param)
  !
  ! Output: nbod       ==> Number of bodies
  !         nbodm      ==> Index of last massive body
  !         pbod(nbod) ==> Mass of bodies
  !
  ! Remarks: Based on io_init_pl
  ! Authors: Hal Levison
  ! Date: 11/21/96
  ! Last revision: 01/10/97
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  character(len = *) :: infile
  type(param_t), intent(in) :: param

  ! Output variables
  integer(ik) :: nbod, nbodm
  type(body_t), allocatable :: pbod(:)

  ! Internal variables
  integer(ik) :: i, j, k, ierr, iu = 7
  integer(ik) :: ngap, igap(NPLMAX)
  real(rk) :: mstar, j2rp2, j4rp4, fwgap(NPLMAX)
  !integer(ik) :: ibad
  !real(rk) :: rhill(NTPMAX), rhrat

  ! Namelist
  namelist / gap_params / ngap, igap, fwgap

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Open particle data file
  call io_open(iu, infile, "OLD", "FORMATTED", ierr)

  ! Read number of bodies, and the mass of the central body
  read(iu,'(1x,i7,1x,1pe22.15)') nbod, mstar

  if(nbod > NTPMAX) then

    write(*,'(a)')      "SWIFT Error:: io_init_pl_symba"
    write(*,'(a,i6)')   " The number of bodies ", nbod
    write(*,'(a,i6)')   " Must be less than ", NTPMAX
    call util_exit(FAILURE)

  end if

  write(*,'(a,i6,/)') " Total number of bodies (including the central body) = ", nbod

  ! Define the size of each vector
  allocate(pbod(nbod), stat = ierr)
  if(ierr /= 0) then

    write(*,'(a)')       "SWIFT Error:: io_init_pl_symba"
    write(*,'(2(a,i6))') " Error allocating pbod vector: ", ierr, " nbod = ", nbod
    call util_exit(FAILURE)

  end if

  ! Store the mass, heliocentric position and velocity of the central body
  pbod(1)%id = 1
  pbod(1)%mass = mstar
  pbod(1)%r = 0.0_rk
  pbod(1)%v = 0.0_rk

  ! These quantities are not used, but initialize to zero just in case
  pbod(1)%rhill = 0.0_rk
  pbod(1)%rphy = 0.0_rk
  pbod(1)%rdrag = 0.0_rk
  pbod(1)%rho = 0.0_rk

  ! For each body read mass, heliocentric position and velocity, Hill radius, physical radius, density and drag radius
  ! depending on the values of param%lclose and param%lgas
  if(param%lclose) then

    if(param%lgas) then

      do j = 2, nbod

        read(iu,'(1x,i7)')          pbod(j)%id
        read(iu,'(5(1x,1pe22.15))') pbod(j)%mass, pbod(j)%rhill, pbod(j)%rphy, pbod(j)%rdrag, pbod(j)%rho
        read(iu,'(3(1x,1pe22.15))') pbod(j)%r
        read(iu,'(3(1x,1pe22.15))') pbod(j)%v

      end do

    else

      do j = 2, nbod

        read(iu,'(1x,i7)')          pbod(j)%id
        read(iu,'(3(1x,1pe22.15))') pbod(j)%mass, pbod(j)%rhill, pbod(j)%rphy
        read(iu,'(3(1x,1pe22.15))') pbod(j)%r
        read(iu,'(3(1x,1pe22.15))') pbod(j)%v

      end do

    end if

  else

    if(param%lgas) then

      do j = 2, nbod

        read(iu,'(1x,i7)')          pbod(j)%id
        read(iu,'(4(1x,1pe22.15))') pbod(j)%mass, pbod(j)%rhill, pbod(j)%rdrag, pbod(j)%rho
        read(iu,'(3(1x,1pe22.15))') pbod(j)%r
        read(iu,'(3(1x,1pe22.15))') pbod(j)%v

      end do

    else

      do j = 2, nbod

        read(iu,'(1x,i7)')          pbod(j)%id
        read(iu,'(2(1x,1pe22.15))') pbod(j)%mass, pbod(j)%rhill
        read(iu,'(3(1x,1pe22.15))') pbod(j)%r
        read(iu,'(3(1x,1pe22.15))') pbod(j)%v

      end do

    end if

  end if

  close(unit = iu)

  ! Find the location of the last massive particle
  call util_nbodm(param, nbod, nbodm, pbod)

  ! Open and read in gap parameters
  call io_open(iu, param%paramfile, "OLD", "FORMATTED", ierr)
  read(unit = iu, nml = gap_params)
  close(unit = iu)

  ! Store the gap information for each massive body
  if(ngap > 0) then

    do i = 1, nbodm

      pbod(i)%lgap = .false.

      do j = 1, ngap

        if(igap(j) == i) then

          pbod(i)%lgap = .true.
          pbod(i)%wgap = fwgap(j)

        end if

      end do

    end do

  else

    do i = 1, nbodm

      pbod(i)%lgap = .false.

    end do

  end if

  ! Check to see if the hills spheres are alright
  !call util_hills(nbod, pbod, rhill)
  !ibad = 0
  !do j = 2, nbod
  !
  !  rhrat = pbod(j)%rhill/rhill(j)
  !  if((rhrat > 2.0_rk) .or. (rhrat < 0.5_rk)) ibad = ibad + 1
  !
  !end do
  !
  !if(ibad /= 0) then
  !
  !  write(*,'(a)')        "SWIFT Warning:: io_init_pl_symba"
  !  write(*,'(a,i6,a,/)') " Hill's spheres are not consistent on ", ibad, " objects"
  !
  !end if

  return
  end subroutine io_init_pl_symba
!!
  subroutine io_dump_param(dparfile, param, t)
  !--------------------------------------------------------------------------
  !         IO_DUMP_PARAM.F90
  !--------------------------------------------------------------------------
  ! Dumps out the parameters for the integration
  !
  ! Input:  dparfile ==> Name of file to dump parameters
  !         param    ==> Global parameters (See swift module & io_init_param)
  !         t        ==> Current time
  !
  ! Authors: Martin Duncan
  ! Date: 03/02/93
  ! Last revision: 05/10/94 HFL
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !--------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: t
  character(len = 24) :: dparfile
  type(param_t), intent(inout) :: param

  ! Internal variables
  logical(lk) :: lgas, lenergy, loblate, lclose
  integer(ik) :: ngap, igap(NPLMAX), i, ierr
  real(rk) :: t0, tstop, dt, dtout, dtdump
  real(rk) :: j2, j4, rmin, rmax, rmaxu, qmin
  real(rk) :: mtiny, rhogas0, gpower, zscale, rgi, rgf, taugas, ca, ce, fwgap(NPLMAX)
  character(len = 24) :: outfile, outtype, outform, outstat

  ! Namelists
  namelist / sim_params / t0, tstop, dt, dtout, dtdump, outfile, outtype, outform, outstat, lenergy
  namelist / stellar_params / loblate, j2, j4, lclose, rmin, rmax, rmaxu, qmin
  namelist / symba_params / mtiny
  namelist / gas_params / lgas, rhogas0, gpower, zscale, rgi, rgf, taugas, ca, ce
  namelist / gap_params / ngap, igap, fwgap

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Initialize igap and fwgap
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i) SHARED(igap, fwgap)
  do i = 1, NPLMAX

    igap(i) = 0
    fwgap(i) = 0.0_rk

  end do
  !$OMP END PARALLEL DO

  ! Open parameter file
  open(unit = 7, file = trim(adjustl(param%paramfile)), status = "OLD")
  read(unit = 7, nml = sim_params)
  read(unit = 7, nml = stellar_params)
  read(unit = 7, nml = symba_params)
  read(unit = 7, nml = gas_params)
  read(unit = 7, nml = gap_params)
  close(unit = 7)

  ! Update t0 and outstat
  t0 = t
  outstat = "APPEND"

  ! Ensure integration time, along with output intervals are integer multiples of the timestep
  if(modulo(tstop,dt) /= 0) tstop = dt*(floor(tstop/dt) + 1)
  if(modulo(dtout,dt) /= 0) dtout = dt*(floor(dtout/dt) + 1)
  if(modulo(dtdump,dt) /= 0) dtdump = dt*(floor(dtdump/dt) + 1)

  ! Open parameter data file for the dump
  open(unit = 8, file = dparfile, status = "UNKNOWN")
  write(unit = 8, nml = sim_params)
  write(unit = 8, nml = stellar_params)
  write(unit = 8, nml = symba_params)
  write(unit = 8, nml = gas_params)
  write(unit = 8, nml = gap_params)
  close(unit = 8)

  ! Store updated value of outstat in param
  param%outstat = outstat

  return
  end subroutine io_dump_param
!!
  subroutine io_dump_pl_symba(dplfile, param, nbod, pbod)
  !--------------------------------------------------------------------------------
  !         IO_DUMP_PL_SYMBA.F90
  !--------------------------------------------------------------------------------
  ! Dumps the body information (e.g. position, velocity, mass, radius, etc.)
  ! for all bodies
  !
  ! Input:  dplfile    ==> Name of file to dump body information
  !         param      ==> Global parameters (See swift module & io_init_param)
  !         nbod       ==> Number of bodies
  !         pbod(nbod) ==> Mass of bodies
  !
  ! Remarks: Based on io_dump_pl.f
  ! Authors: Hal Levison
  ! Date: 01/08/97
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !--------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  character(len = 24) :: dplfile
  type(param_t), intent(in) :: param
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: j, ierr, iu = 7

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Open dump particle data file
  call io_open(iu, dplfile, "UNKNOWN", "FORMATTED", ierr)

  ! Write number of bodies, and the mass of the central body
  write(iu,'(1x,i7,1x,1pe22.15)') nbod, pbod(1)%mass

  ! For each body write mass, heliocentric position and velocity, Hill radius, physical radius, density and drag radius
  ! depending on the values of param%lclose and param%lgas
  if(param%lclose) then

    if(param%lgas) then

      do j = 2, nbod

        write(iu,'(1x,i7)')          pbod(j)%id
        write(iu,'(5(1x,1pe22.15))') pbod(j)%mass, pbod(j)%rhill, pbod(j)%rphy, pbod(j)%rdrag, pbod(j)%rho
        write(iu,'(3(1x,1pe22.15))') pbod(j)%r
        write(iu,'(3(1x,1pe22.15))') pbod(j)%v

      end do

    else

      do j = 2, nbod

        write(iu,'(1x,i7)')          pbod(j)%id
        write(iu,'(3(1x,1pe22.15))') pbod(j)%mass, pbod(j)%rhill, pbod(j)%rphy
        write(iu,'(3(1x,1pe22.15))') pbod(j)%r
        write(iu,'(3(1x,1pe22.15))') pbod(j)%v

      end do

    end if

  else

    if(param%lgas) then

      do j = 2, nbod

        write(iu,'(1x,i7)')          pbod(j)%id
        write(iu,'(4(1x,1pe22.15))') pbod(j)%mass, pbod(j)%rhill, pbod(j)%rdrag, pbod(j)%rho
        write(iu,'(3(1x,1pe22.15))') pbod(j)%r
        write(iu,'(3(1x,1pe22.15))') pbod(j)%v

      end do

    else

      do j = 2, nbod

        write(iu,'(1x,i7)')          pbod(j)%id
        write(iu,'(2(1x,1pe22.15))') pbod(j)%mass, pbod(j)%rhill
        write(iu,'(3(1x,1pe22.15))') pbod(j)%r
        write(iu,'(3(1x,1pe22.15))') pbod(j)%v

      end do

    end if

  end if

  close(unit = iu)

  return
  end subroutine io_dump_pl_symba
!!
  function io_read_hdr(iu, t, npl, ntp, iout_form, out_type) result(ierr)
  !--------------------------------------------------------------------------------
  !         IO_READ_HDR.F90
  !--------------------------------------------------------------------------------
  ! Reads the current time, the number of massive bodies, the number of test
  ! particles and the output form
  !
  ! Now chooses between single or double precision ordinary binary or FXDR binary
  ! output files based on out_type
  !
  ! Input:  iu        ==> Unit number of binary file
  !         out_type  ==> Output type
  !                       - out_type = [ "REAL4" | "REAL8" | "FXDR4" | "FXDR8" ]
  !
  ! Output: t         ==> Current time
  !         npl       ==> Number of massive bodies
  !         ntp       ==> Number of test particles
  !         iout_form ==> Output form
  !                       - iout_form = [ EL | XV ]
  !         ierr      ==> Status flag
  !                       - Zero if read succeeded, non-zero otherwise
  !
  ! Authors: Hal Levison
  ! Date: 11/02/99
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Modelled on SWIFTER's io routines
  !--------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik), intent(in) :: iu
  character(len = *), intent(in) :: out_type

  ! Output variables
  integer(ik) :: npl, ntp, iout_form, ierr
  real(rk) :: t

  ! Internal variables
  integer(ik) :: nvec(NDIM)
  real(real4) :: t4

  !-----------------!
  ! Executable code !
  !-----------------!

  select case(out_type)

    case("REAL4")

      read(iu, iostat = ierr) t4, npl, ntp, iout_form

      if(ierr /= 0) return
      t = real(t4, rk)

    case("REAL8")

      read(iu, iostat = ierr) t, npl, ntp, iout_form
      if(ierr /= 0) return

    case("FXDR4")

      ierr = ixdrreal(iu, t4)
      if(ierr /= 0) return
      t = real(t4, rk)

      ierr = ixdrimat(iu, NDIM, nvec)
      if(ierr /= 0) return
      npl = nvec(1)
      ntp = nvec(2)
      iout_form = nvec(3)

    case("FXDR8")

      ierr = ixdrdouble(iu, t)
      if(ierr /= 0) return

      ierr = ixdrimat(iu, NDIM, nvec)
      if(ierr /= 0) return
      npl = nvec(1)
      ntp = nvec(2)
      iout_form = nvec(3)

    case default

      write(*,'(a)')  "SWIFT Error:: io_read_hdr"
      write(*,'(2a)') " Unknown output type: ", out_type
      call util_exit(FAILURE)

  end select

  return
  end function io_read_hdr
!!
  function io_read_line(iu, id, d1, d2, d3, d4, d5, d6, out_type, mass, radius) result(ierr)
  !-----------------------------------------------------------------------------------------
  !         IO_READ_LINE.F90
  !-----------------------------------------------------------------------------------------
  ! Reads a body identifier, six elements and optionally a mass and radius for a body
  !
  ! Now chooses between single or double precision ordinary binary or FXDR binary output
  ! files based on out_type
  !
  ! Input:  iu                ==> Unit number of binary file
  !         out_type          ==> Output type
  !                               - out = [ "REAL4" | "REAL8" | "FXDR4" | "FXDR8" ]
  !
  ! Output: id                ==> Body idenifier
  !         d1,d2,d3,d4,d5,d6 ==> Body elements (Keplerian or Cartesian)
  !         mass              ==> Mass of body (Optional)
  !         radius            ==> Radius of body (Optional)
  !         ierr              ==> Status flag
  !                               - Zero if read succeeded, non-zero otherwise
  !
  ! Authors: Hal Levison
  ! Date: 11/02/99
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Modelled on SWIFTER's io routines
  !-----------------------------------------------------------------------------------------
  implicit none

  ! Number of elements
  integer(ik), parameter :: NEL = 6

  ! Input variables
  integer(ik), intent(in) :: iu
  character(len = *), intent(in) :: out_type

  ! Output variables
  integer(ik) :: id, ierr
  real(rk) :: d1, d2, d3, d4, d5, d6
  real(rk), optional :: mass, radius

  ! Internal variables
  logical(lk) :: lmass, lradius
  real(real4) :: smass, sradius, svec(NEL)
  real(rk) :: dvec(NEL)

  !-----------------!
  ! Executable code !
  !-----------------!

  lmass = present(mass)
  lradius = present(radius)

  select case(out_type)

    case("REAL4")

      if(lmass .and. lradius) then

        read(iu, iostat = ierr) id, smass, sradius, svec

      else if(lmass .and. (.not. lradius)) then

        read(iu, iostat = ierr) id, smass, svec

      else

        read(iu, iostat = ierr) id, svec

      end if

      if(ierr /= 0) return
      if(lmass) mass = real(smass, rk)
      if(lradius) radius = real(sradius, rk)

      d1 = real(svec(1), rk)
      d2 = real(svec(2), rk)
      d3 = real(svec(3), rk)
      d4 = real(svec(4), rk)
      d5 = real(svec(5), rk)
      d6 = real(svec(6), rk)

    case("REAL8")

      if(lmass .and. lradius) then

        read(iu, iostat = ierr) id, mass, radius, dvec

      else if(lmass .and. (.not. lradius)) then

        read(iu, iostat = ierr) id, mass, dvec

      else

        read(iu, iostat = ierr) id, dvec

      end if

      if(ierr /= 0) return

      d1 = dvec(1)
      d2 = dvec(2)
      d3 = dvec(3)
      d4 = dvec(4)
      d5 = dvec(5)
      d6 = dvec(6)

    case("FXDR4")

      ierr = ixdrint(iu, id)
      if(ierr /= 0) return

      if(lmass) then

        ierr = ixdrreal(iu, smass)
        if(ierr /= 0) return
        mass = real(smass, rk)

        if(lradius) then

          ierr = ixdrreal(iu, sradius)
          if(ierr /= 0) return
          radius = real(sradius, rk)

        end if

      end if

      ierr = ixdrrmat(iu, NEL, svec)
      if(ierr /= 0) return

      d1 = real(svec(1), rk)
      d2 = real(svec(2), rk)
      d3 = real(svec(3), rk)
      d4 = real(svec(4), rk)
      d5 = real(svec(5), rk)
      d6 = real(svec(6), rk)

    case("FXDR8")

      ierr = ixdrint(iu, id)
      if(ierr /= 0) return

      if(lmass) then

        ierr = ixdrdouble(iu, mass)
        if(ierr /= 0) return

        if(lradius) then

          ierr = ixdrdouble(iu, radius)
          if(ierr /= 0) return

        end if

      end if

      ierr = ixdrdmat(iu, NEL, dvec)
      if (ierr /= 0) return

      d1 = dvec(1)
      d2 = dvec(2)
      d3 = dvec(3)
      d4 = dvec(4)
      d5 = dvec(5)
      d6 = dvec(6)

    case default

      write(*,'(a)')  "SWIFT Error:: io_read_line"
      write(*,'(2a)') " Unknown output type: ", out_type
      call util_exit(FAILURE)

  end select

  return
  end function io_read_line
!!
  subroutine io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
  !-------------------------------------------------------------------------------
  !         IO_WRITE_HDR.F90
  !-------------------------------------------------------------------------------
  ! Writes the current time, the number of massive bodies, the number of test
  ! particles and the output form
  !
  ! Now chooses between single or double precision ordinary binary or FXDR binary
  ! output files based on out_type
  !
  ! Input:  iu        ==> Unit number of binary file
  !         t         ==> Current time
  !         npl       ==> Number of massive bodies
  !         ntp       ==> Number of test particles
  !         iout_form ==> Output form
  !                       - iout_form = [ EL | XV ]
  !         out_type  ==> Output type
  !                       - out_type = [ "REAL4" | "REAL8" | "FXDR4" | "FXDR8" ]
  !
  ! Authors: Hal Levison
  ! Date: 11/02/99
  ! Last revision: 09/02/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Modelled on SWIFTER's io routines
  !-------------------------------------------------------------------------------
  implicit none

  ! Parameter
  integer(ik), parameter :: NEL = 3

  ! Input variables
  integer(ik), intent(in) :: iu, npl, ntp, iout_form
  real(rk), intent(in) :: t
  character(len = *), intent(in) :: out_type

  ! Internal variables
  integer(ik) :: ierr, nvec(NEL)
  real(real4) :: t4

  !-----------------!
  ! Executable code !
  !-----------------!

  t4 = real(t, real4)
  nvec = (/ npl, ntp, iout_form /)

  select case(out_type)

    case("REAL4")

      write(iu, iostat = ierr) t4, npl, ntp, iout_form

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_hdr"
        write(*,'(a)') " Unable to write binary file header"
        call util_exit(FAILURE)

      end if

    case("REAL8")

      write(iu, iostat = ierr) t, npl, ntp, iout_form

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_hdr"
        write(*,'(a)') " Unable to write binary file header"
        call util_exit(FAILURE)

      end if

    case("FXDR4")

      ierr = ixdrreal(iu, t4)

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_hdr"
        write(*,'(a)') " Unable to write binary file header"
        call util_exit(FAILURE)

      end if

      ierr = ixdrimat(iu, NEL, nvec)

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_hdr"
        write(*,'(a)') " Unable to write binary file header"
        call util_exit(FAILURE)

      end if

    case("FXDR8")

      ierr = ixdrdouble(iu, t)

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_hdr"
        write(*,'(a)') " Unable to write binary file header"
        call util_exit(FAILURE)

      end if

      ierr = ixdrimat(iu, NEL, nvec)

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_hdr"
        write(*,'(a)') " Unable to write binary file header"
        call util_exit(FAILURE)

      end if

    case default

      write(*,'(a)')  "SWIFT Error:: io_write_hdr"
      write(*,'(2a)') " Unknown output type: ", out_type
      call util_exit(FAILURE)

  end select

  return
  end subroutine io_write_hdr
!!
  subroutine io_write_line(iu, id, d1, d2, d3, d4, d5, d6, out_type, mass, radius)
  !--------------------------------------------------------------------------------------
  !         IO_WRITE_LINE.F90
  !--------------------------------------------------------------------------------------
  ! Writes a body identifier, six elements and optionally a mass and radius for a body
  !
  ! Now chooses between single or double precision ordinary binary or FXDR binary output
  ! files based on out_type
  !
  ! Input:  iu                     ==> Unit number of binary file
  !         id                     ==> Body idenifier
  !         d1, d2, d3, d4, d5, d6 ==> Body elements (Keplerian or Cartesian)
  !         out_type               ==> Output type
  !                                    - out_type = [ "REAL4" | "REAL8" | "FXDR4" | "FXDR8" ]
  !         mass                   ==> Mass of body (Optional)
  !         radius                 ==> Radius of body (Optional)
  !
  ! Authors: Hal Levison
  ! Date: 11/02/99
  ! Last revision: 09/02/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Modelled on SWIFTER's io routines
  !--------------------------------------------------------------------------------------
  implicit none

  ! Number of elements
  integer(ik), parameter :: NEL = 6

  ! Input variables
  integer(ik), intent(in) :: iu, id
  real(rk), intent(in) :: d1, d2, d3, d4, d5, d6
  real(rk), optional, intent(in) :: mass, radius
  character(len = *), intent(in) :: out_type

  ! Internal variables
  logical(lk) :: lmass, lradius
  integer(ik) :: ierr
  real(rk) :: dvec(NEL)
  real(real4) :: svec(NEL), smass, sradius

  !-----------------!
  ! Executable code !
  !-----------------!

  dvec = (/ d1, d2, d3, d4, d5, d6 /)
  svec = real((/ d1, d2, d3, d4, d5, d6 /), real4)

  lmass = present(mass)
  lradius = present(radius)

  if(lmass) smass = real(mass, real4)
  if(lradius) sradius = real(radius, real4)

  select case(out_type)

    case("REAL4")

      if(lmass .and. lradius) then

        write(iu, iostat = ierr) id, smass, sradius, svec

      else if(lmass .and. (.not. lradius)) then

        write(iu, iostat = ierr) id, smass, svec

      else

        write(iu, iostat = ierr) id, svec

      end if

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_line"
        write(*,'(a)') " Unable to write binary file record"
        call util_exit(FAILURE)

      end if

    case("REAL8")

      if(lmass .and. lradius) then

        write(iu, iostat = ierr) id, mass, radius, dvec

      else if(lmass .and. (.not. lradius)) then

        write(iu, iostat = ierr) id, mass, dvec

      else

        write(iu, iostat = ierr) id, dvec

      end if

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_line"
        write(*,'(a)') " Unable to write binary file record"
        call util_exit(FAILURE)

      end if

    case("FXDR4")

      ierr = ixdrint(iu, id)

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_line"
        write(*,'(a)') " Unable to write binary file record"
        call util_exit(FAILURE)

      end if

      if(lmass) then

        ierr = ixdrreal(iu, smass)

        if(ierr < 0) then

          write(*,'(a)') "SWIFT Error:: io_write_line"
          write(*,'(a)') " Unable to write binary file record"
          call util_exit(FAILURE)

        end if

        if(lradius) then

          ierr = ixdrreal(iu, sradius)

          if(ierr < 0) then

            write(*,'(a)') "SWIFT Error:: io_write_line"
            write(*,'(a)') " Unable to write binary file record"
            call util_exit(FAILURE)

          end if

        end if

      end if

      ierr = ixdrrmat(iu, NEL, svec)

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_line"
        write(*,'(a)') " Unable to write binary file record"
        call util_exit(FAILURE)

      end if

    case("FXDR8")

      ierr = ixdrint(iu, id)

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_line"
        write(*,'(a)') " Unable to write binary file record"
        call util_exit(FAILURE)

      end if

      if(lmass) then

        ierr = ixdrdouble(iu, mass)

        if(ierr < 0) then

          write(*,'(a)') "SWIFT Error:: io_write_line"
          write(*,'(a)') " Unable to write binary file record"
          call util_exit(FAILURE)

        end if

        if(lradius) then

          ierr = ixdrdouble(iu, radius)

          if(ierr < 0) then

            write(*,'(a)') "SWIFT Error:: io_write_line"
            write(*,'(a)') " Unable to write binary file record"
            call util_exit(FAILURE)

          end if

        end if

      end if

      ierr = ixdrdmat(iu, NEL, dvec)

      if(ierr < 0) then

        write(*,'(a)') "SWIFT Error:: io_write_line"
        write(*,'(a)') " Unable to write binary file record"
        call util_exit(FAILURE)

      end if

    case default

      write(*,'(a)')  "SWIFT Error:: io_write_line"
      write(*,'(2a)') " Unknown output type: ", out_type
      call util_exit(FAILURE)

  end select

  return
  end subroutine io_write_line
!!
  subroutine io_write_frame(param, t, nbod, pbod)
  !---------------------------------------------------------------------------
  !         IO_WRITE_FRAME.F90
  !---------------------------------------------------------------------------
  ! Write a frame (header plus records for each massive/non-massive body)
  ! to output binary file
  !
  ! Now chooses between an ordinary binary or FXDR binary output
  ! files based on out_type
  !
  ! Input:  param      ==> Global parameter (See swift_module & io_init_param)
  !         t          ==> Current time
  !         nbod       ==> Number of bodies
  !         pbod(nbod) ==> Mass of bodies
  !
  ! Authors: Hal Levison
  ! Date: 11/02/99
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Modelled on SWIFTER io routine
  !---------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik), intent(in) :: nbod
  real(rk), intent(in) :: t
  type(body_t), intent(in) :: pbod(nbod)
  type(param_t), intent(in) :: param

  ! Internal variables
  logical(lk), save :: lfirst = .TRUE.
  logical(lk), save :: lfxdr = .FALSE.
  integer(ik) :: i, j, ierr, ntp, ialpha, iu = 20
  real(rk) :: gm, mu, a(nbod), e(nbod), inc(nbod), capom(nbod), omega(nbod), capm(nbod)

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Is this our first time through?
  if(lfirst) then

    ! Set FXDR flag, to select FXDR routines
    lfxdr = ((param%outtype == "FXDR4") .or. (param%outtype == "FXDR8"))

    ! Are we appending to an existing file?
    if(param%outstat == "APPEND") THEN

      ! Use FXDR routine open routine if lfxdr is set
      if(lfxdr) then

        call io_open_fxdr(param%outfile, "A", .TRUE., iu, ierr)

      ! Otherwise use default open routine
      else

        call io_open(iu, param%outfile, param%outstat, "UNFORMATTED", ierr)

      end if

    ! Create a new file, as long as the file does not already exist
    else if(param%outstat == "NEW") THEN

      ! Use FXDR routine open routine if lfxdr is set
      if(lfxdr) then

        call io_open_fxdr(param%outfile, "R", .TRUE., iu, ierr)

        if(ierr == 0) then

          write(*,'(a)') "SWIFT Error:: io_write_frame"
          write(*,'(a)') " Binary output file already exists"
          call util_exit(FAILURE)

        end if

        call io_open_fxdr(param%outfile, "W", .TRUE., iu, ierr)

      ! Otherwise use default open routine
      else

        call io_open(iu, param%outfile, param%outstat, "UNFORMATTED", ierr)

        if(ierr /= 0) then

          write(*,'(a)') "SWIFT Error:: io_write_frame"
          write(*,'(a)') " Binary output file already exists"
          call util_exit(FAILURE)

        end if

      end if

    ! Otherwise create a new file, even if it means deleting an existing file
    else

      ! Use FXDR routine open routine if lfxdr is set
      if(lfxdr) then

        call io_open_fxdr(param%outfile, "W", .TRUE., iu, ierr)

      ! Otherwise use default open routine
      else

        call io_open(iu, param%outfile, "REPLACE", "UNFORMATTED", ierr)

      end if

    end if

    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: io_write_frame"
      write(*,'(a)') " Unable to open binary output file"
      call util_exit(FAILURE)

    end if

    ! Set lfirst so we do not return here
    lfirst = .FALSE.

  else ! Otherwise append to an existing file

    ! Use FXDR routine open routine if lfxdr is set
    if(lfxdr) then

      call io_open_fxdr(param%outfile, "A", .TRUE., iu, ierr)

    ! Otherwise use default open routine
    else

      call io_open(iu, param%outfile, "APPEND", "UNFORMATTED", ierr)

    end if

    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: io_write_frame"
      write(*,'(a)') " Unable to open binary output file for append"
      call util_exit(FAILURE)

    end if

  end if

  ! Write the header for this frame
  ! NB: We set the number of test particles to zero, this is so we can
  !     use the existing SWIFTER import template in SWIFTVIS which expects
  !     a value for the number of test particles.
  !     Also, the SWIFTER template permits reading an optional mass
  !     (and rdrag?) entry.
  ntp = 0
  call io_write_hdr(iu, t, nbod, ntp, param%ioutform, param%outtype)

  ! Select the output form (also borrowed from SWIFTER)
  select case(param%ioutform)

    ! Output osculating Keplerian orbital elements
    case(EL)

      ! Convert Cartesian positions and velocities to osculating Keplerian orbital elements
      mu = pbod(1)%mass
      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, gm) FIRSTPRIVATE(ialpha, mu) &
      !$OMP SHARED(nbod, pbod, a, e, inc, capom, omega, capm)
      do i = 2, nbod

        gm = mu + pbod(i)%mass
        call orbel_xv2el(pbod(i)%r, pbod(i)%v, gm, ialpha, a(i), e(i), inc(i), capom(i), omega(i), capm(i))

      end do
      !$OMP END PARALLEL DO

      ! Write osculating Keplerian orbital elements, mass and drag radius (depending on the value of param%lgas)
      if(param%lgas) then

        do i = 2, nbod

          call io_write_line(iu, -pbod(i)%id, a(i), e(i), inc(i), capom(i), omega(i), capm(i), param%outtype, &
          mass = pbod(i)%mass, radius = pbod(i)%rdrag)

        end do

      else

        do i = 2, nbod

          call io_write_line(iu, -pbod(i)%id, a(i), e(i), inc(i), capom(i), omega(i), capm(i), param%outtype, &
          mass = pbod(i)%mass)

        end do

      end if

    ! Output Cartesian positions and velocities
    case(XV)

      ! Write Cartesian positions, velocities, mass and drag radius (depending on the value of param%lgas)
      if(param%lgas) then

        do i = 2, nbod

          call io_write_line(iu, -pbod(i)%id, pbod(i)%r(1), pbod(i)%r(2), pbod(i)%r(3), pbod(i)%v(1), pbod(i)%v(2), pbod(i)%v(3), &
          param%outtype, mass = pbod(i)%mass, radius = pbod(i)%rdrag)

        end do

      else

        do i = 2, nbod

          call io_write_line(iu, -pbod(i)%id, pbod(i)%r(1), pbod(i)%r(2), pbod(i)%r(3), pbod(i)%v(1), pbod(i)%v(2), pbod(i)%v(3), &
          param%outtype, mass = pbod(i)%mass)

        end do

      end if

  end select

  ! Close binary using FXDR routine if lfxdr is set
  if(lfxdr) then

    ierr = ixdrclose(iu)

  ! Otherwise use standard routine
  else

    close(unit = iu, iostat = ierr)

  end if

  if(ierr /= 0) then

    write(*,'(a)') "SWIFT Error:: io_write_frame"
    write(*,'(a)') " Unable to close binary output file"
    call util_exit(FAILURE)

  end if

  return
  end subroutine io_write_frame
!!
  subroutine io_write_energy(param, t, nbod, nbodm, pbod, eoffset)
  !--------------------------------------------------------------------------------
  !         IO_WRITE_ENERGY.F90
  !--------------------------------------------------------------------------------
  ! Write the total energy, along with components of the total angular momentum of
  ! the system at the current time to energy.dat
  !
  ! Input:  param      ==> Global parameters (See swift module & io_init_param)
  !         t          ==> Current time
  !         nbod       ==> Number of bodies
  !         nbodm      ==> Number of massive bodies
  !         pbod(nbod) ==> Mass of bodies
  !         eoffset    ==> An energy offset that is added to the energy
  !
  ! Authors: Hal Levison
  ! Date: 03/04/93
  ! Last revision: 12/27/96
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !--------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm
  real(rk) :: t, eoffset
  type(body_t) :: pbod(nbod)
  type(param_t) :: param

  ! Internal variables
  logical(lk), save :: lfirst = .true.
  integer(ik) :: ierr, iu = 40
  real(rk) :: energy, l(NDIM), ke, pot
  character(len = 80) :: outfile = "energy.dat"

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Compute the total energy and angular momentum of the system
  call anal_energy(param, nbod, nbodm, pbod, ke, pot, energy, l)

  ! Add any energy offset
  energy = energy + eoffset

  ! Write current values of the total energy and angular momentum to energy.out
  if(lfirst) then

    call io_open(iu, trim(outfile), param%outstat, "FORMATTED", ierr)

    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: io_write_energy"
      write(*,'(a)') " Could not open energy.dat"
      call util_exit(FAILURE)

    end if

    lfirst = .false.

  else

    call io_open(iu, trim(outfile), "APPEND", "FORMATTED", ierr)

    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: io_write_energy"
      write(*,'(a)') " Could not open energy.dat"
      call util_exit(FAILURE)

    end if

  end if

  write(iu,'(1x,1pe12.5,4(1x,1pe22.15))') t, energy, l

  close(unit = iu)

  return
  end subroutine io_write_energy
!!
  subroutine io_discard_mass(time, pbod, iwhy)
  !------------------------------------------------------------------------------------
  !         IO_DISCARD_MASS.F90
  !------------------------------------------------------------------------------------
  ! Write out information about a discarded body to discard_mass.dat
  !
  ! Input:  time ==> Current time
  !         pbod ==> Database entry of a single body
  !         iwhy ==> Reason why body id was discarded
  !
  ! Authors:  Hal Levison
  ! Date: 12/30/96
  ! Last revision: 09/02/09 CCC - Converted Fortran 90/95 syntax
  !------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: iwhy
  real(rk) :: time
  type(body_t) :: pbod

  ! Internal variables
  integer(ik) :: ierr, iu = 30

  !-----------------!
  ! Executable code !
  !-----------------!

  ! If first call
  if(lfirst_discard) then

    ! Create discard file
    call io_open(iu, "discard_mass.dat", "UNKNOWN", "FORMATTED", ierr)

    ! Print error message if there is a problem
    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: io_discard_mass"
      write(*,'(a)') " Could not open discard output file"
      call util_exit(FAILURE)

    end if

    ! Set lfirst_discard so we don't return here
    lfirst_discard = .false.

  else

    ! Open existing discard file
    call io_open(iu, "discard_mass.dat", "APPEND", "FORMATTED", ierr)

    ! Print error message if there is a problem
    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: io_discard_mass"
      write(*,'(a)') " Could not open discard output file"
      call util_exit(FAILURE)

    end if

  end if

  ! Write discard information to discard file
  write(iu,'(1x,1pe22.15,1x,i4)') time, iwhy
  write(iu,'("-1",1x,i6,1x,2(1x,1pe22.15))') pbod%id, pbod%mass, pbod%rphy
  write(iu,'(3(1x,1pe22.15))') pbod%r
  write(iu,'(3(1x,1pe22.15))') pbod%v

  ! Close discard file
  close(unit = iu)

  return
  end subroutine io_discard_mass
!!
  subroutine io_discard_merge(time, pbod1, pbod2, pbodn)
  !-------------------------------------------------------------------------------
  !         IO_DISCARD_MERGE.F90
  !-------------------------------------------------------------------------------
  ! Write out information about a merger to discard_mass.dat
  !
  ! Input:  time  ==> Current time
  !         pbod1 ==> Database entry of mergee
  !         pbod2 ==> Database entry of merger
  !         pbodn ==> Database entry of merger product
  !
  ! Authors: Hal Levison
  ! Date: 12/30/96
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: time
  type(body_t) :: pbod1, pbod2, pbodn

  ! Internal variables
  integer(ik) :: ierr, iu = 30

  !-----------------!
  ! Executable code !
  !-----------------!

  ! If first call
  if(lfirst_discard) then

    ! Create discard file
    call io_open(iu, "discard_mass.dat", "UNKNOWN", "FORMATTED", ierr)

    ! Print error message if there is a problem
    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: io_discard_merge"
      write(*,'(a)') " Could not open discard output file"
      call util_exit(FAILURE)

    end if

    ! Set lfirst_discard so we don't return here
    lfirst_discard = .false.

  else

    ! Open existing discard file
    call io_open(iu, "discard_mass.dat", "APPEND", "FORMATTED", ierr)

    ! Print error message if there is a problem
    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: io_discard_merge"
      write(*,'(a)') " Could not open discard output file"
      call util_exit(FAILURE)

    end if

  end if

  write(iu,'(1x,1pe22.15,"  2")') time

  write(iu,'("-1",1x,i6,1x,2(1x,1pe22.15))') pbod1%id, pbod1%mass, pbod1%rphy
  write(iu,'(3(1x,1pe22.15))') pbod1%r
  write(iu,'(3(1x,1pe22.15))') pbod1%v

  write(iu,'("-1",1x,i6,1x,2(1x,1pe22.15))') pbod2%id, pbod2%mass, pbod2%rphy
  write(iu,'(3(1x,1pe22.15))') pbod2%r
  write(iu,'(3(1x,1pe22.15))') pbod2%v

  write(iu,'("+1",1x,i6,1x,2(1x,1pe22.15))') pbodn%id, pbodn%mass, pbodn%rphy
  write(iu,'(3(1x,1pe22.15))') pbodn%r
  write(iu,'(3(1x,1pe22.15))') pbodn%r

  close(unit = iu)

  return
  end subroutine io_discard_merge
!!
!!
end module io
