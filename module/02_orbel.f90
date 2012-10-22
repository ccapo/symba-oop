module orbel
! Module for routines in orbel directory
use swift
implicit none

contains
!!
!!
  subroutine orbel_scget(angle, sx, cx)
  !----------------------------------------------------------------------
  !         ORBEL_SCGET.F90
  !----------------------------------------------------------------------
  ! Given an angle, efficiently compute sine and cosine
  !
  ! Input:  angle ==> Angle [radians]
  !
  ! Output: sx    ==> Sine of angle
  !         cx    ==> Cosine of angle
  !
  ! REMARKS: The HP 700 series won't return correct answers for sin
  !          and cos if the angle is bigger than 3e7. We first reduce it
  !          to the range [0,2*pi) and use the sqrt rather than cos
  !          since its faster
  !          BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
  ! AUTHOR: M. Duncan.
  ! DATE WRITTEN: May 6, 1992
  ! REVISIONS: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !----------------------------------------------------------------------
  implicit none

  ! Input variable
  real(rk) :: angle

  ! Output variables
  real(rk) :: sx, cx

  ! Internal variables
  integer(ik) :: nper
  real(rk) :: x

  !-----------------!
  ! Executable code !
  !-----------------!

  nper = angle/TWOPI
  x = angle - nper*TWOPI
  if(x < 0.0_rk) x = x + TWOPI
  sx = sin(x)
  cx = sqrt(1.0_rk - sx*sx)
  if((x > PIBY2) .and. (x < PI3BY2)) cx = -cx

  return
  end subroutine orbel_scget
!!
  subroutine orbel_schget(angle, shx, chx)
  !----------------------------------------------------------------------
  !         ORBEL_SCHGET.F90
  !----------------------------------------------------------------------
  ! Given an angle, efficiently compute sinh and cosh
  !
  ! Input:  angle ==> Angle [radians]
  !
  ! Output: shx   ==>  sinh(angle)
  !         chx   ==>  cosh(angle)
  !
  ! REMARKS: Based on the routine SCGET for sine's and cosine's
  !          We use the sqrt rather than cosh since its faster
  !          BE SURE THE ANGLE IS IN RADIANS AND IT CAN'T BE LARGER
  !          THAN 300 OR OVERFLOWS WILL OCCUR!
  ! AUTHOR: M. Duncan.
  ! DATE WRITTEN: May 6, 1992
  ! REVISIONS: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !----------------------------------------------------------------------
  implicit none

  ! Input variable
  real(rk) :: angle

  ! Output variables
  real(rk) :: shx, chx

  !-----------------!
  ! Executable code !
  !-----------------!

  shx = sinh(angle)
  chx = sqrt(1.0_rk + shx*shx)

  return
  end subroutine orbel_schget
!!
  function orbel_esolmd(e, m) result(ea)
  !----------------------------------------------------------------------
  !         ORBEL_ESOLMD.F90
  !----------------------------------------------------------------------
  ! Solves Kepler's equation given the eccentricity and mean anomaly
  ! using Danby's quartic convergence algorithm
  !
  ! Input:  e  ==> Eccentricity
  !         m  ==> Mean anomaly
  !
  ! Output: ea ==> Eccentric anomaly
  !
  ! REMARKS: Only good for small eccentricity since it only iterates once,
  !          and also does not put M or E between 0.0 AND TWOPI
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: May 7, 1992
  ! REVISIONS: 02/26/93 HFL
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !----------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: e, m

  ! Output variable
  real(rk) :: ea

  ! Internal variables
  real(rk) :: sm, cm, sea, cea
  real(rk) :: es, ec, f, fp, fpp, fppp, dea

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Begin with a guess accurate to O(e^3)
  call orbel_scget(m, sm, cm)
  ea = m + e*sm*(1.0_rk + e*(cm + e*(1.0_rk - 1.5_rk*sm*sm)))

  ! Go through one iteration for improved estimate
  call orbel_scget(ea, sea, cea)
  es = e*sea
  ec = e*cea
  f = ea - es  - m
  fp = 1.0_rk - ec
  fpp = es
  fppp = ec
  dea = -f/fp
  dea = -f/(fp + 0.5_rk*fpp*dea)
  dea = -f/(fp + 0.5_rk*fpp*dea + fppp*dea*dea/6.0_rk)
  ea = ea + dea

  return
  end function orbel_esolmd
!!
  function orbel_eget(e, m) result(ea)
  !-------------------------------------------------------------------------
  !         ORBEL_EGET.F90
  !-------------------------------------------------------------------------
  ! Solves Kepler's equation given the eccentricity and mean anomaly
  ! using Danby's quartic convergence algorithm
  !
  ! Input:  e  ==> Eccentricity
  !         m  ==> Mean anomaly
  !
  ! Output: ea ==> Eccentric anomaly
  !
  ! REMARKS: For results very near roundoff, give it M between 0 and 2*pi
  !        : If e < 0.18 use ORBEL_ESOLMD for speed and sufficient accuracy
  !        : If e > 0.8 use ORBEL_EHIE which may not converge fast enough
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: May 7, 1992.
  ! REVISIONS: May 21, 1992 - Now have it go through EXACTLY two
  !                           iterations with the premise that it will
  !                           only be called if we have an ellipse with
  !                           0.18 <= e <= 0.8
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: e, m

  ! Output variable
  real(rk) :: ea

  ! Internal variables
  real(rk) :: sm, cm, sea, cea
  real(rk) :: es, ec, f, fp, fpp, fppp, dea

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Begin with a guess accurate to O(e^3)
  call orbel_scget(m, sm, cm)
  ea = m + e*sm*(1.0_rk + e*(cm + e*(1.0_rk - 1.5_rk*sm**2)))

  ! Go through one iteration for improved estimate
  call orbel_scget(ea, sea, cea)
  es = e*sea
  ec = e*cea
  f = ea - es  - m
  fp = 1.0_rk - ec
  fpp = es
  fppp = ec
  dea = -f/fp
  dea = -f/(fp + 0.5_rk*fpp*dea)
  dea = -f/(fp + 0.5_rk*fpp*dea + fppp*(dea**2)/6.0_rk)
  ea = ea + dea

  ! Do another iteration
  call orbel_scget(ea, sea, cea)
  es = e*sea
  ec = e*cea
  f = ea - es  - m
  fp = 1.0_rk - ec
  fpp = es
  fppp = ec
  dea = -f/fp
  dea = -f/(fp + 0.5_rk*fpp*dea)
  dea = -f/(fp + 0.5_rk*fpp*dea + fppp*(dea**2)/6.0_rk)
  ea = ea + dea

  return
  end function orbel_eget
!!
  function orbel_ehie(e, m) result(ea)
  !----------------------------------------------------------------------
  !         ORBEL_EHIE.F90
  !----------------------------------------------------------------------
  ! Solves Kepler's equation given the eccentricity and mean anomaly
  ! using three iterations of Danby's quartic convergence algorithm
  !
  ! N.B. The equation is f(x) = x - e*sin(x + M), where E = x + M
  !
  ! Input:  e  ==> Eccentricity
  !         m  ==> Mean anomaly
  !
  ! Output: ea ==> Eccentric anomaly
  !
  ! REMARKS: First guess is very good for e near 1.0
  !        : Modifies M so that both E and M are in range (0, TWOPI)
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: May 25, 1992
  ! REVISIONS: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !----------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: e, m

  ! Output variable
  real(rk) :: ea

  ! Internal variables
  integer(ik), parameter :: NMAX = 3
  logical(lk) :: lflag
  integer(ik) :: nper, niter
  real(rk) :: dx, x, sa, ca, esa, eca, f, fp

  !-----------------!
  ! Executable code !
  !-----------------!

  ! In this section, bring M into the range (0, TWOPI) and if
  ! the result is greater than PI, solve for (TWOPI - M)
  lflag = .false.
  nper = m/TWOPI
  m = m - nper*TWOPI
  if(m < 0.0_rk) m = m + TWOPI

  if(m > PI) then

    m = TWOPI - m
    lflag = .true.

  end if

  ! Make a first guess that works well for e near 1.0
  x = (6.0_rk*m)**(1.0_rk/3.0_rk) - m

  ! Iteration loop
  do niter = 1, NMAX

    call orbel_scget(x + m, sa, ca)
    esa = e*sa
    eca = e*ca
    f = x - esa
    fp = 1.0_rk - eca
    dx = -f/fp
    dx = -f/(fp + 0.5_rk*dx*esa)
    dx = -f/(fp + 0.5_rk*dx*(esa + eca*dx/3.0_rk))
    x = x + dx

  end do

  ea = x + m

  if(lflag) then

    ea = TWOPI - ea
    m = TWOPI - m

  end if

  return
  end function orbel_ehie
!!
  function orbel_ehybrid(e, m) result(ea)
  !----------------------------------------------------------------------
  !         ORBEL_EHYBRID.F90
  !----------------------------------------------------------------------
  ! Solves Kepler's equation given the eccentricity and mean anomaly
  !
  ! N.B.: If e < 0.18 uses fast routine ORBEL_ESOLMD
  !       If 0.18 <= e <= 0.8 uses ORBEL_EGET
  !       If e > 0.8 uses ORBEL_EHIE
  !
  ! Input:  e  ==> Eccentricity
  !         m  ==> Mean anomaly
  !
  ! Output: ea ==> Eccentric anomaly
  !
  ! REMARKS: Only ORBEL_EHIE brings M and E into range (0,TWOPI)
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: May 25, 1992
  ! REVISIONS: 02/26/93 HFL
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !----------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: e, m

  ! Output variable
  real(rk) :: ea

  !-----------------!
  ! Executable code !
  !-----------------!

  if(e < 0.18_rk) then

    ea = orbel_esolmd(e, m)

  else

    if(e <= 0.8_rk) then

      ea = orbel_eget(e, m)

    else

      ea = orbel_ehie(e, m)

    end if

  end if

  return
  end function orbel_ehybrid
!!
  function orbel_flon(e, capn) result(ea)
  !--------------------------------------------------------------------------
  !         ORBEL_FLON.F90
  !--------------------------------------------------------------------------
  ! Solves Kepler's equation for hyperbola using a power series for N in
  ! terms of F and Newton's method
  !
  ! N.B. Only good for low values of N (e.g. N < 0.636*e - 0.6)
  !
  ! Input:  e    ==> Eccentricity
  !         capn ==> Hyperbola mean anomaly
  !
  ! Output: ea   ==> Eccentric anomaly
  !
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: May 26, 1992
  ! REVISIONS: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !--------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: e, capn

  ! Output variable
  real(rk) :: ea

  ! Internal variables
  integer(ik), parameter :: imax = 10
  real(rk), parameter :: a11 = 156.0_rk, a9 = 17160.0_rk, a7 = 1235520.0_rk
  real(rk), parameter :: a5 = 51891840.0_rk, a3 = 1037836800.0_rk
  real(rk), parameter :: b11 = 11.0_rk*a11, b9 = 9.0_rk*a9, b7 = 7.0_rk*a7
  real(rk), parameter :: b5 = 5.0_rk*a5, b3 = 3.0_rk*a3
  logical(lk) :: lflag
  integer(ik) :: i
  real(rk) :: a, b, sq, biga, bigb
  real(rk) :: x, x2
  real(rk) :: f, fp, dx
  real(rk) :: diff, a0, a1, b1

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Function to solve Kepler's equation for F (here called ea) for given e and CAPN.
  ! Only good for smallish CAPN
  lflag = .false.
  if(capn < 0.0_rk) then

    capn = -capn
    lflag = .true.

  end if

  a1 = 6227020800.0_rk*(1.0_rk - 1.0_rk/e)
  a0 = -6227020800.0_rk*capn/e
  b1 = a1

  ! Set lflag if capn < 0., in which case solve for -capn
  ! and change the sign of the final answer for F.
  ! Begin with a reasonable guess based on solving the cubic for small F
  a = 6.0_rk*(e - 1.0_rk)/e
  b = -6.0_rk*capn/e
  sq = sqrt(0.25_rk*b*b + a*a*a/27.0_rk)
  biga = (sq - 0.5*b)**(1.0_rk/3.0_rk)
  bigb = -(sq + 0.5*b)**(1.0_rk/3.0_rk)
  x = biga + bigb
  ea = x

  ! If capn is tiny (or zero) no need to go further than cubic even for e = 1
  i = 1
  if(capn >= TINY_NUMBER) then

    do i = 1, IMAX

      x2 = x*x
      f = a0 + x*(a1 + x2*(a3 + x2*(a5 + x2*(a7 + x2*(a9 + x2*(a11 + x2))))))
      fp = b1 + x2*(b3 + x2*(b5 + x2*(b7 + x2*(b9 + x2*(b11 + 13.0_rk*x2)))))
      dx = -f/fp
      ea = x + dx

      ! If we have converged here there's no point in going on
      if(abs(dx) <= TINY_NUMBER) exit
      x = ea

    end do

  end if

  ! Check if capn was originally negative
  if(lflag) then

    ea = -ea
    capn = -capn

  end if

  ! Abnormal return here - we've gone thru the loop IMAX times without convergence
  if(i >= imax) then

    diff = e*sinh(ea) - ea - capn
    write(*,'(a)') "SWIFT Warning:: orbel_flon"
    write(*,'(a)') " Returning without complete convergence"
    write(*,'(a,3(1pe14.6))') "N, F, ecc*sinh(F) - F - N: ", capn, ea, diff

  end if

  return
  end function orbel_flon
!!
  function orbel_fget(e, capn) result(ea)
  !------------------------------------------------------------------------
  !         ORBEL_FGET.F90
  !------------------------------------------------------------------------
  ! Solves Kepler's equation for hyperbola using hybrid approach
  !
  ! Input:  e    ==> Eccentricity
  !         capn ==> Hyperbola mean anomaly
  !
  ! Output: ea   ==> Eccentric anomaly
  !
  ! ALGORITHM: Based on pp. 70 - 72 of Fitzpatrick's book "Principles of
  !            Celestical Mechanics"
  !          : Quartic convergence from Danby's book
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: May 11, 1992
  ! REVISIONS: 02/26/93 HFL
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: e, capn

  ! Output variable
  real(rk) :: ea

  ! Internal variables
  integer(ik), parameter :: imax = 10
  logical(lk) :: lflag
  integer(ik) :: i
  real(rk) :: tmp, x, shx, chx
  real(rk) :: esh, ech, f, fp, fpp, fppp, dx

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Begin with a guess proposed by Danby
  if(capn < 0.0_rk) then

    tmp = -2.0_rk*capn/e + 1.8_rk
    x = -log(tmp)

  else

    tmp = 2.0_rk*capn/e + 1.8_rk
    x = log(tmp)

  end if

  ea = x

  do i = 1, IMAX

    call orbel_schget(x, shx, chx)
    esh = e*shx
    ech = e*chx
    f = esh - x - capn
    !write(*,*) 'i, x, f: ', i, x, f
    fp = ech - 1.0_rk
    fpp = esh
    fppp = ech
    dx = -f/fp
    dx = -f/(fp + 0.5_rk*dx*fpp)
    dx = -f/(fp + 0.5_rk*dx*fpp + dx*dx*fppp/6.0_rk)
    ea = x + dx

    ! If we have converged here there's no point in going on
    if(abs(dx) <= TINY_NUMBER) then

      lflag = .false.
      exit

    else

      lflag = .true.

    end if

    x = ea

  end do

  if(lflag) then

    write(*,'(a)') "SWIFT Warning:: orbel_fget"
    write(*,'(a)') " Returning without complete convergence"

  end if

  return
  end function orbel_fget
!!
  function orbel_fhybrid(e, capn) result(ea)
  !------------------------------------------------------------------------
  !         ORBEL_FHYBRID.F90
  !------------------------------------------------------------------------
  ! Solves Kepler's equation for hyperbola using hybrid approach
  !
  ! N.B. If abs(N) < 0.636*e - 0.6, uses ORBEL_FLON
  !      Otherwise uses ORBEL_FGET
  !
  ! Input:  e    ==> Eccentricity
  !         capn ==> Hyperbola mean anomaly
  !
  ! Output: ea   ==> Eccentric anomaly
  !
  !     ALGORITHM: For abs(N) < 0.636*ecc - 0.6, use FLON
  !                For larger N, uses FGET
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: May 26, 1992
  ! REVISIONS: 02/26/93 HFL
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: e, capn

  ! Output variable
  real(rk) :: ea

  ! Internal variable
  real(rk) :: abn

  !-----------------!
  ! Executable code !
  !-----------------!

  abn = capn
  if(capn < 0.0_rk) abn = -abn

  if(abn < 0.636_rk*e - 0.6_rk) then

    ea = orbel_flon(e, capn)

  else

    ea = orbel_fget(e, capn)

  endif

  return
  end function orbel_fhybrid
!!
  function orbel_zget(q) result(ea)
  !----------------------------------------------------------------------
  !         ORBEL_ZGET.F90
  !----------------------------------------------------------------------
  ! Solves the equivalent of Kepler's equation for a parabola given Q
  ! (Fitzpatrick notation)
  !
  ! Input:  q  ==> Parabola mean anomaly
  !
  ! Output: ea ==> Eccentric anomaly
  !
  ! ALGORITHM: pp. 70 - 72 of Fitzpatrick's book "Princ. of Cel. Mech."
  ! REMARKS: For a parabola we can solve analytically
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: May 11, 1992
  ! REVISIONS: May 27, 1992(?) - corrected it for negative Q, and use
  !                              power series for small Q
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !----------------------------------------------------------------------
  implicit none

  ! Input variable
  real(rk) :: q

  ! Output variable
  real(rk) :: ea

  ! Internal variables
  logical(lk) :: lflag
  real(rk) :: x, tmp

  !-----------------!
  ! Executable code !
  !-----------------!

  lflag = .false.
  if(q < 0.0_rk) then

    q = -q
    lflag = .true.

  endif

  if(q < 1.0e-3_rk) then

    ea = q*(1.0_rk - q*q*(1.0_rk - q*q)/3.0_rk)

  else

    x = 0.5_rk*(3.0_rk*q + sqrt(9.0_rk*q*q + 4.0_rk))
    tmp = x**(1.0_rk/3.0_rk)
    ea = tmp - 1.0_rk/tmp

  end if

  if(lflag) then

    ea = -ea
    q = -q

  end if

  return
  end function orbel_zget
!!
  subroutine orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm, r, v)
  !--------------------------------------------------------------------------------------
  !         ORBEL_EL2XV.F90
  !--------------------------------------------------------------------------------------
  ! Computes the Cartesian position and velocity for a body given the central mass,
  ! conic section type, and orbital elements (See Fitzpatrick "Principles of Cel. Mech.")
  !
  ! Input:  gm      ==> G times central mass
  !         ialpha  ==> Conic section type
  !                     - ialpha = +1 (Hyperbola)
  !                     - ialpha =  0 (Parabola)
  !                     - ialpha = -1 (Ellipse)
  !         a       ==> Semi-major axis or pericentric distance if a parabola [AU]
  !         e       ==> Eccentricity
  !         inc     ==> Inclination [radians]
  !         capom   ==> Longitude of ascending node [radians]
  !         omega   ==> Argument of perihelion [radians]
  !         capm    ==> Mean anomoly [radians]
  !
  ! Output: r(NDIM) ==> Cartesian position of body
  !         v(NDIM) ==> Cartesian velocity of body
  !
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: May 11, 1992
  ! REVISIONS: May 26, 1992(?) - Uses better Kepler solver for ellipses (ORBEL_EHYBRID)
  !                              and hyperbolae (ORBEL_FHYBRID)
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !--------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: ialpha
  real(rk) :: gm, a, e, inc, capom, omega, capm

  ! Output variables
  real(rk) :: r(NDIM), v(NDIM)

  ! Internal variables
  logical(lk) :: lnotcircle, lnotellipse, lnothyperbola
  real(rk) :: cape, capf, zpara
  real(rk) :: sp, cp, so, co, si, ci
  real(rk) :: d11, d12, d13, d21, d22, d23
  real(rk) :: scap, ccap, shcap, chcap
  real(rk) :: sqe, sqgma, xfac1, xfac2, ri, vfac1, vfac2

  !-----------------!
  ! Executable code !
  !-----------------!

  if(e < 0.0_rk) then

    write(*,'(a)') "SWIFT Warning:: orbel_el2xv"
    write(*,'(a)') " Found e < 0.0, forcing e = 0.0"
    e = 0.0_rk

  end if

  ! Check for inconsistencies between ialpha and e
  lnotcircle = (ialpha == 0) .and. (abs(e - 1.0_rk) > TINY_NUMBER)
  lnotellipse = (ialpha < 0) .and. (e > 1.0_rk)
  lnothyperbola = (ialpha > 0) .and. (e < 1.0_rk)
  if(lnotcircle .or. lnotellipse .or. lnothyperbola)  then

    write(*,'(a)')              "SWIFT Error:: orbel_el2xv"
    write(*,'(a,i6,a,1pe14.6)') " ialpha and e are inconsistent: ialpha = ", ialpha, " e = ", e
    write(*,'(/,a,f3.1,a)')     "Terminating SWIFT (Version: ", VER_NUM, ") due to ERROR!!!"
    write(*,'(a)')              "------------------------------------------------"

  end if

  ! Generate rotation matrices (See p. 42 from Fitzpatrick)
  call orbel_scget(omega, sp, cp)
  call orbel_scget(capom, so, co)
  call orbel_scget(inc, si, ci)
  d11 = cp*co - sp*so*ci
  d12 = cp*so + sp*co*ci
  d13 = sp*si
  d21 = -sp*co - cp*so*ci
  d22 = -sp*so + cp*co*ci
  d23 = cp*si

  ! Get the other quantities depending on orbit type (i.e. IALPHA)
  select case(ialpha)

    case(-1) ! Elliptical Orbit

      cape = orbel_ehybrid(e, capm)
      call orbel_scget(cape, scap, ccap)
      sqe = sqrt(1.0_rk - e*e)
      sqgma = sqrt(gm*a)
      xfac1 = a*(ccap - e)
      xfac2 = a*sqe*scap
      ri = 1.0_rk/(a*(1.0_rk - e*ccap))
      vfac1 = -ri*sqgma*scap
      vfac2 = ri*sqgma*sqe*ccap

    case(1) ! Hyperbolic Orbit

      capf = orbel_fhybrid(e, capm)
      call orbel_schget(capf, shcap, chcap)
      sqe = sqrt(e*e - 1.0_rk)
      sqgma = sqrt(gm*a)
      xfac1 = a*(e - chcap)
      xfac2 = a*sqe*shcap
      ri = 1.0_rk/(a*(e*chcap - 1.0_rk))
      vfac1 = -ri*sqgma*shcap
      vfac2 = ri*sqgma*sqe*chcap

    case(0) ! Parabolic Orbit

      zpara = orbel_zget(capm)
      sqgma = sqrt(2.0_rk*gm*a)
      xfac1 = a*(1.0_rk - zpara*zpara)
      xfac2 = 2.0_rk*a*zpara
      ri = 1.0_rk/(a*(1.0_rk + zpara*zpara))
      vfac1 = -ri*sqgma*zpara
      vfac2 = ri*sqgma

    case default

      write(*,'(a)')          "SWIFT Error:: orbel_el2xv"
      write(*,'(a,i6)')       " Unrecognized value of ialpha: ", ialpha
      write(*,'(/,a,f4.1,a)') "Terminating SWIFT (Version:", VER_NUM, ") due to ERROR!!!"
      write(*,'(a)')          "------------------------------------------------"
      stop

  end select

  r(1) = d11*xfac1 + d21*xfac2
  r(2) = d12*xfac1 + d22*xfac2
  r(3) = d13*xfac1 + d23*xfac2
  v(1) = d11*vfac1 + d21*vfac2
  v(2) = d12*vfac1 + d22*vfac2
  v(3) = d13*vfac1 + d23*vfac2

  return
  end subroutine orbel_el2xv
!!
  subroutine orbel_xv2el(r, v, gmsum, ialpha, a, e, inc, capom, omega, capm)
  !---------------------------------------------------------------------------------
  !         ORBEL_XV2EL.F90
  !---------------------------------------------------------------------------------
  ! Given the Cartesian position and velocity of a body, computes the osculating
  ! orbital elements
  !
  ! Input:  r(NDIM) ==> Cartesian position of body
  !         v(NDIM) ==> Cartesian velocity of body
  !         gmsum   ==> G times sum of masses of the body and central body
  !
  ! Output: ialpha  ==> Conic section type
  !                     - ialpha = +1 (Hyperbola)
  !                     - ialpha =  0 (Parabola)
  !                     - ialpha = -1 (Ellipse)
  !         a       ==> semi-major axis or pericentric distance if a parabola [AU]
  !         e       ==> Eccentricity
  !         inc     ==> Inclination [radians]
  !         capom   ==> Longitude of ascending node [radians]
  !         omega   ==> Argument of perihelion [radians]
  !         capm    ==> Mean anomoly [radians]
  !
  ! ALGORITHM: See  p. 70 of Fitzpatrick's "Priciples of Cel. Mech."
  ! REMARKS: If the inclination INC is less than TINY_NUMBER, we arbitrarily set
  !          the longitude of the ascending node CAPOM to 0.0, so the ascending
  !          node is then along the X axis
  !        : If the eccentricity E is less than sqrt(TINY_NUMBER), we arbitrarily
  !          set the argument of perihelion OMEGA to 0.0
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: May 8, 1992
  ! REVISIONS: 02/04/00
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !---------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: gmsum, r(NDIM), v(NDIM)

  ! Output variables
  integer(ik) :: ialpha
  real(rk) :: a, e, inc, capom, omega, capm

  ! Internal variables
  real(rk) :: h(NDIM), h2, habs, rdist, v2, vdotr, energy, fac, face, cape, capf, tmpf
  real(rk) :: cw, sw, w, u

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Compute the angular momentum H, and thereby the inclination INC
  h(1) = r(2)*v(3) - r(3)*v(2)
  h(2) = r(3)*v(1) - r(1)*v(3)
  h(3) = r(1)*v(2) - r(2)*v(1)
  h2 = sum(h**2)
  habs = sqrt(h2)
  if(h(3) > habs) h(3) = habs ! Hal's fix
  inc = acos(h(3)/habs)

  ! Compute longitude of ascending node CAPOM and the argument of latitude u
  fac = sqrt(sum(h(1:2)**2))/habs

  if(fac < TINY_NUMBER) then

    capom = 0.0_rk
    u = atan2(r(2), r(1))
    if(abs(inc - PI) < 10.0_rk*TINY_NUMBER) u = -u

  else

    capom = atan2(h(1), -h(2))
    u = atan2(r(3)/sin(inc), r(1)*cos(capom) + r(2)*sin(capom))

  end if

  if(capom < 0.0_rk) capom = capom + TWOPI
  if(u < 0.0_rk) u = u + TWOPI

  ! Compute the distance rdist and velocity squared V2, and the dot
  ! product RDOTV, the energy per unit mass ENERGY
  rdist = sqrt(sum(r**2))
  v2 = sum(v**2)
  vdotr = dot_product(r, v)
  energy = 0.5_rk*v2 - gmsum/rdist

  ! Determine type of conic section and label it via IALPHA
  if(abs(energy*rdist/gmsum) < sqrt(TINY_NUMBER)) then

    ialpha = 0

  else

    if(energy < 0.0_rk) ialpha = -1
    if(energy > 0.0_rk) ialpha = 1

  end if

  ! Depending on the conic type, determine the remaining elements
  select case(ialpha)

    case(-1) ! Elliptical Orbit

      a = -0.5_rk*gmsum/energy
      fac = 1.0_rk - h2/(gmsum*a)

      if(fac > TINY_NUMBER) then

        e = sqrt(fac)
        face = (a - rdist)/(a*e)

        ! Apr. 16/93: watch for case where face is slightly outside unity
        if(face > 1.0_rk) then

          cape = 0.0_rk

        else

          if(face > -1.0_rk) then

            cape = acos(face)

          else

            cape = PI

          end if

        end if

        if(vdotr < 0.0_rk) cape = TWOPI - cape
        cw = (cos(cape) - e)/(1.0_rk - e*cos(cape))
        sw = sqrt(1.0_rk - e*e)*sin(cape)/(1.0_rk - e*cos(cape))
        w = atan2(sw, cw)
        if(w < 0.0_rk) w = w + TWOPI

      else

        e = 0.0_rk
        w = u
        cape = u

      end if

      capm = cape - e*sin(cape)
      omega = u - w
      if(omega < 0.0_rk) omega = omega + TWOPI
      omega = omega - TWOPI*int(omega/TWOPI)

    case(1) ! Hyperbolic Orbit

      a = 0.5_rk*gmsum/energy
      fac = h2/(gmsum*a)

      if(fac > TINY_NUMBER) then

        e = sqrt(1.0_rk + fac)
        tmpf = (a + rdist)/(a*e)
        if(tmpf < 1.0_rk) tmpf = 1.0_rk
        capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))
        if(vdotr < 0.0_rk) capf = -capf
        cw = (e - cosh(capf))/(e*cosh(capf) - 1.0_rk)
        sw = sqrt(e*e - 1.0_rk)*sinh(capf)/(e*cosh(capf) - 1.0_rk)
        w = atan2(sw, cw)
        if(w < 0.0_rk) w = w + TWOPI

      else

        ! We only get here if a hyperbola is essentially a parabola
        ! so we calculate e and w accordingly to avoid singularities
        e = 1.0_rk
        tmpf = 0.5_rk*h2/gmsum
        w = acos(2.0_rk*tmpf/rdist - 1.0_rk)
        if(vdotr < 0.0_rk) w = TWOPI - w
        tmpf = (a + rdist)/(a*e)
        capf = log(tmpf + sqrt(tmpf**2 - 1.0_rk))

      end if

      capm = e*sinh(capf) - capf
      omega = u - w
      if(omega < 0.0_rk) omega = omega + TWOPI
      omega = omega - TWOPI*int(omega/TWOPI)

    case(0) ! Parabolic Orbit
            ! N.B. In this case we use "a" to mean pericentric distance

      a =  0.5_rk*h2/gmsum
      e = 1.0_rk
      w = acos(2.0*a/rdist - 1.0_rk)
      if(vdotr < 0.0_rk) w = TWOPI - w
      tmpf = tan(0.5_rk*w)
      capm = tmpf*(1.0_rk + tmpf*tmpf/3.0_rk)
      omega = u - w
      if(omega < 0.0_rk) omega = omega + TWOPI
      omega = omega - TWOPI*int(omega/TWOPI)

    case default

      write(*,'(a)')          "SWIFT Error:: orbel_xv2el"
      write(*,'(a,i6)')       " Unrecognized value of ialpha: ", ialpha
      write(*,'(/,a,f4.1,a)') "Terminating SWIFT (Version:", VER_NUM, ") due to ERROR!!!"
      write(*,'(a)')          "------------------------------------------------"
      stop

  end select

  return
  end subroutine orbel_xv2el
!!
  subroutine orbel_xv2aeq(r, v, gmsum, ialpha, a, e, q)
  !-------------------------------------------------------------------------------
  !         ORBEL_XV2AEQ.F90
  !-------------------------------------------------------------------------------
  ! Given the Cartesian position and velocity of a body, compute the osculating
  ! orbital elements a, e, and q *only*
  !
  ! Input:  r(NDIM) ==> Cartesian position of body
  !         v(NDIM) ==> Cartesian velocity of body
  !         gmsum   ==> G times sum of masses of body and central body
  !
  ! Output: ialpha  ==> Conic section type
  !                     - ialpha = +1 (Hyperbola)
  !                     - ialpha =  0 (Parabola)
  !                     - ialpha = -1 (Ellipse)
  !         a       ==> Semi-major axis or pericentric distance if a parabola [AU]
  !         e       ==> Eccentricity
  !         q       ==> Perihelion distance (i.e. q = a(1 - e) ) [AU]
  !
  ! ALGORITHM: See p. 70 of Fitzpatrick's "Priciples of Cel. Mech."
  ! AUTHOR: L. Dones
  ! DATE WRITTEN: February 24, 1994
  ! REVISIONS: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: gmsum, r(NDIM), v(NDIM)

  ! Output variables
  integer(ik) :: ialpha
  real(rk) :: a, e, q

  ! Internal variables
  real(rk) :: h(NDIM), h2, rdist, v2, energy, fac

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Compute the angular momentum H, and thereby the inclination INC

  h(1) = r(2)*v(3) - r(3)*v(2)
  h(2) = r(3)*v(1) - r(1)*v(3)
  h(3) = r(1)*v(2) - r(2)*v(1)
  h2 = sum(h**2)

  ! Compute the distance rdist and velocity squared V2, and the dot
  ! product RDOTV, the energy per unit mass ENERGY
  rdist = sqrt(sum(r**2))
  v2 = sum(v**2)
  energy = 0.5_rk*v2 - gmsum/rdist

  ! Determine type of conic section and label it via IALPHA
  if(abs(energy*rdist/gmsum) < sqrt(TINY_NUMBER)) then

    ialpha = 0

  else

    if(energy < 0.0_rk) ialpha = -1
    if(energy > 0.0_rk) ialpha = 1

  end if

  ! Depending on the conic type, determine the remaining elements
  select case(ialpha)

    case(-1) ! Elliptical Orbit

      a = -0.5_rk*gmsum/energy
      fac = 1.0_rk - h2/(gmsum*a)

      if(fac > TINY_NUMBER) then

        e = sqrt(fac)

      else

        e = 0.0_rk

      end if

      q = a*(1.0_rk - e)

    case(1) ! Hyperbolic Orbit

      a = 0.5_rk*gmsum/energy
      fac = h2/(gmsum*a)

      if(fac > TINY_NUMBER) then

        e = sqrt(1.0_rk + fac)

        ! Have to insert minus sign in expression for q because this code
        ! takes a > 0, even for a hyperbola
        q = -a*(1.0_rk - e)

      else

        ! We only get here if a hyperbola is essentially a parabola
        ! so we calculate e accordingly to avoid singularities
        e = 1.0_rk
        q = 0.5_rk*h2/gmsum

      end if

    case(0) ! Parabolic Orbit
            ! N.B. - In this case "a", which is formally infinite,
            !        is arbitrarily set equal to the pericentric distance q

      a = 0.5_rk*h2/gmsum
      e = 1.0_rk
      q = a

    case default

      write(*,'(a)')          "SWIFT Error:: orbel_xv2aeq"
      write(*,'(a,i6)')       " Unrecognized value of ialpha: ", ialpha
      write(*,'(/,a,f4.1,a)') "Terminating SWIFT (Version:", VER_NUM, ") due to ERROR!!!"
      write(*,'(a)')          "------------------------------------------------"
      stop

  end select

  return
  end subroutine orbel_xv2aeq
!!
  subroutine orbel_xv2aqt(r, v, gmsum, ialpha, a, q, capm, tperi)
  !---------------------------------------------------------------------------------
  !         ORBEL_XV2AQT.F90
  !---------------------------------------------------------------------------------
  ! Given the Cartesian position and velocity of a body, compute the osculating
  ! orbital elements a, q, capm and time of perihelion passage *only*
  !
  ! Input:  r(NDIM) ==> Cartesian position of body
  !         v(NDIM) ==> Cartesian velocity of body
  !         gmsum   ==> G times sum of masses of the body and central body
  !
  ! Output: ialpha  ==> Conic section type
  !                     - ialpha = +1 (Hyperbola)
  !                     - ialpha =  0 (Parabola)
  !                     - ialpha = -1 (Ellipse)
  !         a       ==> Semi-major axis or pericentric distance if a parabola [AU]
  !         q       ==> Perihelion distance (e.g. q = a(1 - e) ) [AU]
  !         capm    ==> Mean anomoly [radians]
  !         tperi   ==> Time to next or last perihelion, which ever is less [year]
  !
  ! ALGORITHM: See  p. 70 of Fitzpatrick's "Priciples of Cel. Mech."
  ! AUTHOR: Hal Levison
  ! DATE WRITTEN: 08/07/01
  ! REMARKS: The tperi may not be correct for parabolic orbits.
  !          I think it is OK, but beware!
  ! REVISIONS: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !---------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: gmsum, r(NDIM), v(NDIM)

  ! Output variables
  integer(ik) :: ialpha
  real(rk) :: a, q, capm, tperi

  ! Internal variables
  real(rk) :: h(NDIM), h2, rdist, v2, energy, fac, vdotr, cape, e
  real(rk) :: capf, tmpf, meanmo, face, w

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Compute the angular momentum H, and thereby the inclination INC
  h(1) = r(2)*v(3) - r(3)*v(2)
  h(2) = r(3)*v(1) - r(1)*v(3)
  h(3) = r(1)*v(2) - r(2)*v(1)
  h2 = sum(h**2)

  ! Compute the distance rdist and velocity squared V2, and the dot
  ! product RDOTV, the energy per unit mass ENERGY
  rdist = sqrt(sum(r**2))
  v2 = sum(v**2)
  energy = 0.5_rk*v2 - gmsum/rdist
  vdotr = dot_product(r, v)

  ! Determine type of conic section and label it via IALPHA
  if(abs(energy*rdist/gmsum) < sqrt(TINY_NUMBER)) then

    ialpha = 0

  else

    if(energy < 0.0_rk) ialpha = -1
    if(energy > 0.0_rk) ialpha = 1

  end if

  ! Depending on the conic type, determine the remaining elements
  select case(ialpha)

    case(-1) ! Elliptical Orbits

      a = -0.5_rk*gmsum/energy
      fac = 1.0_rk - h2/(gmsum*a)

      if(fac > TINY_NUMBER) then

        e = sqrt(fac)
        face = (a - rdist)/(a*e)

        ! Apr. 16/93: watch for case where face is slightly outside unity
        if(face > 1.0_rk) then

          cape = 0.0_rk

        else

          if(face > -1.0_rk) then

            cape = acos(face)

          else

            cape = PI

          end if

        end if

        if(vdotr < 0.0_rk) cape = TWOPI - cape

      else

        e = 0.0_rk
        cape = 0.0_rk

      end if

      capm = cape - e*sin(cape)
      q = a*(1.0_rk - e)

    case(1) ! Hyperbolic Orbit

      a = 0.5_rk*gmsum/energy
      fac = h2/(gmsum*a)

      if(fac > TINY_NUMBER) then

        e = sqrt(1.0_rk + fac)

        ! Have to insert minus sign in expression for q because this code
        ! takes a > 0, even for a hyperbola
        q = -a*(1.0_rk - e)

        tmpf = (a + rdist)/(a*e)
        if(tmpf < 1.0) tmpf = 1.0_rk
        capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))
        if(vdotr < 0.0_rk) capf = -capf

      else

        ! We only get here if a hyperbola is essentially a parabola
        ! so we calculate e accordingly to avoid singularities
        e = 1.0_rk
        q = 0.5_rk*h2/gmsum

        tmpf = (a + rdist)/(a*e)
        capf = log(tmpf + sqrt(tmpf*tmpf - 1.0_rk))

      end if

      capm = e*sinh(capf) - capf

    case(0) ! Parabolic Orbit
            ! N.B. - In this case "a", which is formally infinite,
            !        is arbitrarily set equal to the pericentric distance q

      a = 0.5_rk*h2/gmsum
      e = 1.0_rk
      q = a
      w = acos(2.0_rk*a/rdist - 1.0_rk)
      if(vdotr < 0.0_rk) w = TWOPI - w
      tmpf = tan(0.5_rk*w)
      capm = tmpf*(1.0_rk + tmpf*tmpf/3.0_rk)

    case default

      write(*,'(a)')          "SWIFT Error:: orbel_xv2aqt"
      write(*,'(a,i6)')       " Unrecognized value of ialpha: ", ialpha
      write(*,'(/,a,f4.1,a)') "Terminating SWIFT (Version:", VER_NUM, ") due to ERROR!!!"
      write(*,'(a)')          "------------------------------------------------"
      stop

  end select

  meanmo = sqrt(gmsum/(a*a*a))

  if((capm < PI) .or. (ialpha >= 0)) then

    tperi = -capm/meanmo

  else

    tperi = -1.0_rk*(capm - TWOPI)/meanmo

  end if

  return
  end subroutine orbel_xv2aqt
!!
!!
end module orbel
