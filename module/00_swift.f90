module swift
! Module for SWIFT
!
! Author: Hal Levison
! Date: 2/2/93
! Last revision: 3/7/93
!              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
!                             - Uses new data types to simplify code, and
!                               make it easier to add new particle quantities
!                               or global parameters/flags
implicit none

! Type definitions

! Defines precision for integers and reals
integer, parameter :: integer2 = selected_int_kind(4)
integer, parameter :: integer4 = selected_int_kind(9)
integer, parameter :: real4 = selected_real_kind(6)
integer, parameter :: real8 = selected_real_kind(15)

! The user's choice between single and double precision for reals, and range for integers
integer, parameter :: ik = integer4
integer, parameter :: rk = real8
integer, parameter :: lk = kind(.true.)

! Version of SWIFT
real(rk), parameter :: VER_NUM = 3.0_rk

! Return flags
logical(lk), parameter :: SUCCESS = .true.
logical(lk), parameter :: FAILURE = .false.

! Maximum array size
integer(ik), parameter :: NDIM = 3        ! Number of spatial dimensions (update if string theory proves true)
integer(ik), parameter :: NPLMAX = 1024   ! Max. number of planets, including the Sun (2**7)
integer(ik), parameter :: NTPMAX = 524288 ! Max. number of test particles (2**19)

! Size of the test particle integer status flag
integer(ik), parameter :: NSTATP = 3
integer(ik), parameter :: NSTAT = NSTATP + NPLMAX - 1 ! include one @ planet

! Size of the test particle integer status flag
integer(ik), parameter :: NSTATR = NSTAT              ! io_init_tp assumes nstat == nstatr

! Convergence criteria for danby
real(rk), parameter :: DANBYAC = 1.0e-14_rk, DANBYB = 1.0e-13_rk

! Loop limits in the Laguerre attempts
integer(ik), parameter :: NLAG1 = 50, NLAG2 = 400

! A small number
real(rk), parameter :: TINY_NUMBER = 4.0e-15_rk

! Trignometric stuff
real(rk), parameter :: PI = 3.14159265358979324_rk
real(rk), parameter :: TWOPI = 2.0_rk*PI
real(rk), parameter :: PIBY2 = 0.5_rk*PI
real(rk), parameter :: PI3BY2 = 1.5_rk*PI
real(rk), parameter :: DEGRAD = 180.0_rk/PI

! Simulation Parameters
real(rk), parameter :: MSUN = TWOPI**2
real(rk), parameter :: MEARTH = 3.0e-6_rk*MSUN

! Symbolic names for binary output file contents
integer(ik), parameter :: EL = 1
integer(ik), parameter :: XV = 2

! OpenMP Parameters
integer(ik), save :: nthreads = 1
integer(ik), parameter :: NTHRESHOLD = 1000

! Global parameter type definition
type param_t

  ! Coordinate Parameters (only one is true at any time)
  logical(lk) :: lhelio
  logical(lk) :: lbary
  logical(lk) :: lcanon

  ! Simulation Parameters
  logical(lk) :: lenergy
  integer(ik) :: ioutform ! [ EL | XV ]
  real(rk) :: t0
  real(rk) :: tstop
  real(rk) :: dt
  real(rk) :: dtout
  real(rk) :: dtdump
  character(len = 24) :: paramfile ! Input parmater filename
  character(len = 24) :: outfile   ! Output binary filename
  character(len = 24) :: outtype   ! [ "REAL4" | "REAL8" | "FXDR4" | "FXDR8" ]
  character(len = 24) :: outstat   ! [ "NEW" | "APPEND" | "UNKNOWN" ]

  ! Oblateness Parameters
  logical(lk) :: loblate
  real(rk) :: j2rp2
  real(rk) :: j4rp4

  ! Close Encounter Parameters
  logical(lk) :: lclose
  logical(lk) :: laccrete
  real(rk) :: rmin
  real(rk) :: rmax
  real(rk) :: rmaxu
  real(rk) :: qmin

  ! SyMBA Parameters
  real(rk) :: mtiny

  ! Gas Parameters
  logical(lk) :: lgas
  real(rk) :: rhogas0
  real(rk) :: gpower
  real(rk) :: zscale
  real(rk) :: rgi
  real(rk) :: rgf
  real(rk) :: taugas
  real(rk) :: ca
  real(rk) :: ce
  real(rk) :: sigma0
  real(rk) :: spower

end type param_t

! Particle type definition
type body_t

  ! Body identifier
  integer(ik) :: id

  ! Physical properties
  real(rk) :: mass, rhill, rphy, rdrag, rho

  ! Cartesian position, velocity and acceleration
  real(rk), dimension(NDIM) :: r, v, a

  ! Heliocentric positions, velocities and accelerations
  !real(rk), dimension(NDIM) :: rh, vh, ah

  ! Barycentric positions and velocities
  !real(rk), dimension(NDIM) :: rb, vb

  ! Perhelion integer status flag
  integer(ik) :: iperi

  ! Encounter flag
  logical(lk) :: lenc

  ! Recursion level and max recursion level
  integer(ik) :: rlevel, rlevelmax
  
  ! Gap flag
  logical(lk) :: lgap

  ! Location and width of gap in gas disk
  real(rk) :: rgap, wgap

end type body_t

! Encounter list type definition
type enc_t

  ! Number of encounters
  integer(ik) :: nenc

  ! The encounter list
  integer(ik), pointer :: ienc(:)

  ! The lvdotr list
  logical(lk), pointer :: lvdotr(:)

end type enc_t

end module swift
