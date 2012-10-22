module drift
! Module for routines in drift directory
use swift
use orbel, only: orbel_scget
implicit none

contains
!!
!!
  subroutine drift_kepu_stumpff(x, c0, c1, c2, c3)
  !-------------------------------------------------------------------------
  !         DRIFT_KEPU_STUMPFF.F90
  !-------------------------------------------------------------------------
  ! Calculates the Stumpff functions (see Danby p. 172, equations 6.9.15)
  !
  ! Input:  x              ==> Argument
  !
  ! Output: c0, c1, c2, c3 ==> c's from pp. 171 - 172 of Danby
  !
  ! Author: Hal Levison
  ! Date: 02/03/93
  ! Last revision: 02/03/93
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------
  implicit none

  ! Input variable
  real(rk) :: x

  ! Output variables
  real(rk) :: c0, c1, c2, c3

  ! Internal variables
  real(rk), parameter :: xm = 0.1_rk
  integer(ik) :: n, i

  !-----------------!
  ! Executable code !
  !-----------------!

  n = 0
  do while(abs(x) >= xm)

    n = n + 1
    x = x/4.0_rk

  end do

  c2 = (1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x/182.0_rk)/132.0_rk)/90.0_rk)/56.0_rk)/30.0_rk) &
       /12.0_rk)/2.0_rk
  c3 = (1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x*(1.0_rk - x/210.0_rk)/156.0_rk)/110.0_rk)/72.0_rk)/42.0_rk) &
       /20.0_rk)/6.0_rk
  c1 = 1.0_rk - x*c3
  c0 = 1.0_rk - x*c2

  if(n /= 0) then

    do i = n, 1, -1

      c3 = (c2 + c0*c3)/4.0_rk
      c2 = (c1**2)/2.0_rk
      c1 = c0*c1
      c0 = 2.0*c0**2 - 1.0_rk
      x = 4.0_rk*x

    end do

  end if

  return
  end subroutine drift_kepu_stumpff
!!
  subroutine drift_kepu_p3solve(dtau0, r0, mu, alpha, u, s, iflag)
  !----------------------------------------------------------------------------------------
  !         DRIFT_KEPU_P3SOLVE.F90
  !----------------------------------------------------------------------------------------
  ! Returns the real root of cubic often found in solving Kepler's problem in universal
  ! variables.
  !
  ! Input:  dtau0 ==> Time step
  !         r0    ==> Distance between central body and paritcle
  !         mu    ==> Reduced mass of system
  !         alpha ==> Twice the binding energy
  !         u     ==> Dot product of velocity and radial vectors
  !
  ! Output: s     ==> Solution of cubic equation for the universal variable
  !         iflag ==> Status flag
  !                   - Zero if real root found, non-zero otherwise
  !
  ! Author: Martin Duncan
  ! Date: March 12/93
  ! Last revision: March 12/93
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !----------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: dtau0, r0, mu, alpha, u

  ! Output variables
  integer(ik) :: iflag
  real(rk) :: s

  ! Internal variables
  real(rk) :: denom, a0, a1, a2, q, r, sq2, sq, p1, p2

  !-----------------!
  ! Executable code !
  !-----------------!

  denom = (mu - alpha*r0)/6.0_rk
  a2 = 0.5_rk*u/denom
  a1 = r0/denom
  a0 = -dtau0/denom

  q = (a1 - (a2**2)/3.0_rk)/3.0_rk
  r = (a1*a2 - 3.0_rk*a0)/6.0_rk - (a2**3)/27.0_rk
  sq2 = q**3 + r**2

  if(sq2 >= 0.0_rk) then

    sq = sqrt(sq2)

    if((r + sq) <= 0.0_rk) then

      p1 =  -(-(r + sq))**(1.0_rk/3.0_rk)

    else

      p1 = (r + sq)**(1.0_rk/3.0_rk)

    end if

    if((r - sq) <= 0.0_rk) then

      p2 =  -(-(r - sq))**(1.0_rk/3.0_rk)

    else

      p2 = (r - sq)**(1.0_rk/3.0_rk)

    end if

    iflag = 0
    s = p1 + p2 - a2/3.0_rk

  else

    iflag = 1
    s = 0.0_rk

  end if

  return
  end subroutine drift_kepu_p3solve
!!
  subroutine drift_kepmd(dm, es, ec, x, s, c)
  !-------------------------------------------------------------------------------
  !         DRIFT_KEPMD.F90
  !-------------------------------------------------------------------------------
  ! Solves Kepler's equation in difference form for an ellipse, given SMALL dm
  ! and SMALL eccentricity.  See DRIFT_DAN.F90 for the criteria.
  !
  ! ** WARNING **
  ! BUILT FOR SPEED: Does not check how well the original equation is solved!
  ! Can do that in the calling routine by checking how close the quantity is to
  ! zero: (x - ec*s + es*(1.0 - c) - dm)
  !
  ! Input:  dm     ==> Increment in mean anomaly M
  !         es, ec ==> Eccentricity times sine and cosine of E_0
  !
  ! Output: x      ==> Solution to Kepler's difference equation
  !         s, c   ==> Sine and cosine of x
  !
  ! Updated: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: dm, es, ec

  ! Output variables
  real(rk) :: x, s, c

  ! Internal variables
  real(rk), parameter :: A0 = 39916800.0_rk, A1 = 6652800.0_rk, A2 = 332640.0_rk, A3 = 7920.0_rk, A4 = 110.0_rk
  real(rk) :: dx
  real(rk) :: fac1, fac2, q, y
  real(rk) :: f, fp, fpp, fppp

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Calculate the initial guess for the root
  fac1 = 1.0_rk/(1.0_rk - ec)
  q = fac1*dm
  fac2 = fac1*es**2 - ec/3.0_rk
  x = q*(1.0_rk - 0.5_rk*fac1*q*(es - q*fac2))

  ! Excellent approximation to sine and cosine of x for small x
  y = x**2
  s = x*(A0 - y*(A1 - y*(A2 - y*(A3 - y*(A4 - y)))))/A0
  c = sqrt(1.0_rk - s**2)

  ! Compute better value for the root using quartic Newton method
  f = x - ec*s + es*(1.0_rk - c) - dm
  fp = 1.0_rk - ec*c + es*s
  fpp = ec*s + es*c
  fppp = ec*c - es*s
  dx = -f/fp
  dx = -f/(fp + 0.5_rk*fpp*dx)
  dx = -f/(fp + 0.5_rk*fpp*dx + fppp*(dx**2)/6.0_rk)
  x = x + dx

  ! Excellent approximation to sine and cosine of x for small x
  y = x**2
  s = x*(A0 - y*(A1 - y*(A2 - y*(A3 - y*(A4 - y)))))/A0
  c = sqrt(1.0_rk - s**2)

  return
  end subroutine drift_kepmd
!!
  function drift_kepu_guess(dtau0, r0, mu, alpha, u) result(s)
  !------------------------------------------------------------------------------------
  !         DRIFT_KEPU_GUESS.F90
  !------------------------------------------------------------------------------------
  ! Initial guess for solving Kepler's equation using universal variables
  !
  ! Input:  dtau0 ==> Time step
  !         r0    ==> Distance between central body and paritcle
  !         mu    ==> Reduced mass of system
  !         alpha ==> Twice the binding energy
  !         u     ==> Angular momentum
  !
  ! Output: s     ==> Initial guess for the value of universal variable
  !
  ! Author: Hal Levison & Martin Duncan
  ! Date: 03/12/93
  ! Last revision: April 6/93
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: dtau0, r0, mu, alpha, u

  ! Output variable
  real(rk) :: s

  ! Internal variables
  integer(ik) :: iflag
  real(rk) :: y, sy, cy, sigma, es
  real(rk) :: x, a
  real(rk) :: en, ec, e

  !-----------------!
  ! Executable code !
  !-----------------!

  if(alpha > 0.0_rk) then

    ! Find initial guess for elliptic motion
    if(dtau0/r0 <= 0.4_rk)  then

      s = dtau0/r0 - (u*dtau0**2)/(2.0_rk*r0**3)

    else

      a = mu/alpha
      en = sqrt(mu/(a**3))
      ec = 1.0_rk - r0/a
      es = u/(en*a**2)
      e = sqrt(ec**2 + es**2)
      y = en*dtau0 - es
      call orbel_scget(y, sy, cy)
      sigma = sign(1.0_rk, (es*cy + ec*sy))
      x = y + 0.85_rk*sigma*e
      s = x/sqrt(alpha)

    end if

  else

    ! Find initial guess for hyperbolic motion.
    call drift_kepu_p3solve(dtau0, r0, mu, alpha, u, s, iflag)

    if(iflag /= 0) s = dtau0/r0

  end if

  return
  end function drift_kepu_guess
!!
  subroutine drift_kepu_new(s, dtau0, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
  !----------------------------------------------------------------------------------
  !         DRIFT_KEPU_NEW.F90
  !----------------------------------------------------------------------------------
  ! Solving Kepler's equation in universal variables using NEWTON'S METHOD
  !
  ! Input:  s          ==> Inital value of universal variable
  !         dtau0      ==> Time step
  !         r0         ==> Distance between central body and paritcle
  !         mu         ==> Reduced mass of system
  !         alpha      ==> Twice binding energy
  !         u          ==> Angular momentum
  !
  ! Output: s          ==> Final value of universal variable
  !         fp         ==> f' from p. 170 of Danby
  !         c1, c2, c3 ==> c's from pp. 171 - 172 of Danby
  !         iflag      ==> Status flag
  !                        - Zero if converged, non-zero otherwise
  !
  ! Author: Hal Levison
  ! Date: 02/03/93
  ! Last revision: 04/21/93
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !----------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: s, dtau0, r0, mu, alpha, u

  ! Output variables
  integer(ik) :: iflag
  real(rk) :: fp, c1, c2, c3

  ! Internal variables
  integer(ik), parameter :: nmax = 7
  integer(ik) :: i
  real(rk) :: x, c0, ds
  real(rk) :: f, fpp, fppp, fdt

  !-----------------!
  ! Executable code !
  !-----------------!

  do i = 1, nmax

    x = alpha*s**2
    call drift_kepu_stumpff(x, c0, c1, c2, c3)
    c1 = c1*s
    c2 = c2*s**2
    c3 = c3*s**3
    f = r0*c1 + u*c2 + mu*c3 - dtau0
    fp = r0*c0 + u*c1 + mu*c2
    fpp = (mu - r0*alpha)*c1 + u*c0
    fppp = (mu - r0*alpha)*c0 - u*alpha*c1
    ds = -f/fp
    ds = -f/(fp + 0.5_rk*fpp*ds)
    ds = -f/(fp + 0.5_rk*fpp*ds + fppp*(ds**2)/6.0_rk)
    s = s + ds
    fdt = f/dtau0

    ! Quartic convergence test
    if(fdt**2 < DANBYB**2) then

      iflag = 0 ! Newton's method succeeded
      exit

    else

      iflag = 1 ! Newton's method failed

    end if

  end do

  return
  end subroutine drift_kepu_new
!!
  function drift_kepu_fchk(dtau0, r0, mu, alpha, u, s) result(f)
  !-----------------------------------------------------------------------------
  !         DRIFT_KEPU_FCHK.F90
  !-----------------------------------------------------------------------------
  ! Returns the value of the function f of which we are trying to find the root
  ! in universal variables.  So if f is close to zero, then the root is correct.
  !
  ! Input:  dtau0 ==> Time step
  !         r0    ==> Distance between central body and particle
  !         mu    ==> Reduced mass of system
  !         alpha ==> Twice the binding energy
  !         u     ==> Dot product of velocity and radial vectors
  !         s     ==> Approximate root of f
  !
  ! Output: f     ==> Function value
  !
  ! Author:  Martin Duncan
  ! Date:    March 12/93
  ! Last revision: March 12/93
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-----------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: dtau0, r0, mu, alpha, u, s

  ! Output variable
  real(rk) :: f

  ! Internal variables
  real(rk) :: x, c0, c1, c2, c3

  !-----------------!
  ! Executable code !
  !-----------------!

  x = alpha*s**2
  call drift_kepu_stumpff(x, c0, c1, c2, c3)
  c1 = c1*s
  c2 = c2*s**2
  c3 = c3*s**3
  f = r0*c1 + u*c2 + mu*c3 - dtau0

  return
  end function drift_kepu_fchk
!!
  subroutine drift_kepu_lag(s, dtau0, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
  !-----------------------------------------------------------------------------------
  !         DRIFT_KEPU_LAG.F90
  !-----------------------------------------------------------------------------------
  ! Solves Kepler's equation in universal variables using LAGUERRE'S METHOD
  !
  ! Input:  s          ==> Inital value of universal variable
  !         dtau0      ==> Time step
  !         r0         ==> Distance between central body and paritcle
  !         mu         ==> Reduced mass of system
  !         alpha      ==> Twice the binding ergery energy
  !         u          ==> Angular momentum
  !
  ! Output: s          ==> Final value of universal variable
  !         fp         ==> f' from p. 170 of Danby
  !         c1, c2, c3 ==> c's from pp. 171 - 172 of Danby
  !         iflag      ==> Status flag
  !                        - Zero if converged, non-zero otherwise
  !
  ! Author: Hal Levison
  ! Date: 02/03/93
  ! Last revision: 04/21/93
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-----------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: s, dtau0, r0, mu, alpha, u

  ! Output variables
  integer(ik) :: iflag
  real(rk) :: fp, c1, c2, c3

  ! Internal variables
  real(rk), parameter :: ln = 5.0_rk
  integer(ik) :: i, nmax
  real(rk) :: x, fpp, ds, c0, f
  real(rk) :: fdt

  !-----------------!
  ! Executable code !
  !-----------------!

  ! To get the close approach, we needed to take a large number of iterations if alpha < 0
  if(alpha < 0.0_rk) then

    nmax = NLAG2 + 1

  else

    nmax = NLAG1 + 1

  end if

  ! Start Laguerre's method
  do i = 1, nmax

    x = alpha*s**2
    call drift_kepu_stumpff(x, c0, c1, c2, c3)
    c1 = c1*s
    c2 = c2*s**2
    c3 = c3*s**3
    f = r0*c1 + u*c2 + mu*c3 - dtau0
    fp = r0*c0 + u*c1 + mu*c2
    fpp = (mu - r0*alpha)*c1 + u*c0
    ds = -ln*f/(fp + sign(1.0_rk, fp)*sqrt(abs(((ln - 1.0_rk)*fp)**2 - (ln - 1.0_rk)*ln*f*fpp)))
    s = s + ds
    fdt = f/dtau0

    ! Quartic convergence test
    if(fdt**2 < DANBYB**2) then

      iflag = 0 ! Laguerre's method succeeded
      exit

    else

      iflag = 2 ! Laguerre's method failed

    end if

  end do

  return
  end subroutine drift_kepu_lag
!!
  subroutine drift_kepu(dtau0, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
  !---------------------------------------------------------------------------------
  !         DRIFT_KEPU.F90
  !---------------------------------------------------------------------------------
  ! Solves Kepler's equation using universal variables.
  !
  ! Input:  dtau0      ==> Time step
  !         r0         ==> Distance between central body and paritcle
  !         mu         ==> Reduced mass of system
  !         alpha      ==> Twice the binding energy
  !         u          ==> Angular momentun
  !
  ! Output: fp         ==> f' from p. 170 of Danby
  !         c1, c2, c3 ==> c's from pp. 171 - 172 of Danby
  !         iflag      ==> Status flag
  !                        - Zero if converged, non-zero otherwise
  !
  ! Author:  Hal Levison
  ! Date:    2/3/93
  ! Last revision: 02/03/93
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !---------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: dtau0, r0, mu, alpha, u

  ! Output variables
  integer(ik) :: iflag
  real(rk) :: fp, c1, c2, c3

  ! Internal variables
  real(rk) :: s, st, fo, fn

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Compute an initial guess for the solution to Kepler's equation
  s = drift_kepu_guess(dtau0, r0, mu, alpha, u)

  ! Store initial guess for possible use later in
  ! Laguerre's method, in case Newton's method fails
  st = s

  call drift_kepu_new(s, dtau0, r0, mu, alpha, u, fp, c1, c2, c3, iflag)

  if(iflag /= 0) then

    fo = drift_kepu_fchk(dtau0, r0, mu, alpha, u, st)
    fn = drift_kepu_fchk(dtau0, r0, mu, alpha, u, s)

    if(abs(fo) < abs(fn)) s = st

    call drift_kepu_lag(s, dtau0, r0, mu, alpha, u, fp, c1, c2, c3, iflag)

  end if

  return
  end subroutine drift_kepu
!!
  subroutine drift_dan(mu, r0, v0, dtau0, iflag)
  !----------------------------------------------------------------------------
  !         DRIFT_DAN.F90
  !----------------------------------------------------------------------------
  ! Drift a body along its orbit using the Danby algorithm, and determines
  ! which variables to use
  !
  ! Input:  mu       ==> Mass of body
  !         r0(NDIM) ==> Initial jacobi position of body
  !         v0(NDIM) ==> Initial jacobi velocity of body
  !         dtau0    ==> Time step
  !
  ! Output: r0(NDIM) ==> Final jacobi position of body
  !         v0(NDIM) ==> Final jacobi position of body
  !         iflag    ==> Status flag
  !                      - Zero if converged, non-zero otherwise
  !
  ! Authors: Hal Levison & Martin Duncan
  ! Date: 02/10/93
  ! Last revision: April 6/93 - MD adds dtau and keeps dtau0 unchanged
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !----------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk), intent(in) :: mu, dtau0

  ! Input and Output variables
  real(rk), intent(inout) :: r0(NDIM), v0(NDIM)

  ! Output variable
  integer(ik), intent(out) :: iflag

  ! Internal variables
  real(rk) :: rtmp(NDIM), vtmp(NDIM), dtau
  real(rk) :: f, g, fdot, gdot
  real(rk) :: c1, c2, c3
  real(rk) :: u, alpha, fp, r, v2
  real(rk) :: a, a2, en
  real(rk) :: dm, ec, es, e2, xkep
  real(rk) :: fchk, s, c

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Set dtau = dtau0 to be sure timestep is not altered while solving for new coords.
  iflag = 0
  dtau = dtau0
  r = sqrt(sum(r0**2))
  v2 = sum(v0**2)
  u = dot_product(r0, v0)
  alpha = 2.0_rk*mu/r - v2

  if(alpha > 0.0_rk) then

    a = mu/alpha
    a2 = a**2
    en = sqrt(mu/(a*a2))
    ec = 1.0_rk - r/a
    es = u/(en*a2)
    e2 = ec**2 + es**2
    dm = en*dtau - TWOPI*int((en*dtau)/TWOPI)
    dtau = dm/en
	
    if((e2 <= 0.36_rk) .and. (dm**2 <= 0.16_rk) .and. (e2*dm**2 < 0.0016_rk)) then

      call drift_kepmd(dm, es, ec, xkep, s, c)
      fchk = xkep - ec*s + es*(1.0_rk - c) - dm

      ! Check for convergence
      if(fchk**2 <= DANBYB**2) then

        fp = 1.0_rk - ec*c + es*s
        f = (a/r)*(c - 1.0_rk) + 1.0_rk
        g = dtau + (s - xkep)/en
        fdot = -(a/(r*fp))*en*s
        gdot = (c - 1.0_rk)/fp + 1.0_rk

        rtmp = r0*f + v0*g
        vtmp = r0*fdot + v0*gdot

        r0 = rtmp
        v0 = vtmp

        iflag = 0 ! Success

      else

        iflag = 1 ! Failure

      end if
	  
    else

      call drift_kepu(dtau, r, mu, alpha, u, fp, c1, c2, c3, iflag)

      if(iflag == 0) then

        f = 1.0_rk - (mu/r)*c2
        g = dtau - mu*c3
        fdot = -(mu/(fp*r))*c1
        gdot = 1.0_rk - (mu/fp)*c2

        rtmp = r0*f + v0*g
        vtmp = r0*fdot + v0*gdot

        r0 = rtmp
        v0 = vtmp

      end if

    end if
	
  else

    call drift_kepu(dtau, r, mu, alpha, u, fp, c1, c2, c3, iflag)

    if(iflag == 0) then

      f = 1.0_rk - (mu/r)*c2
      g = dtau - mu*c3
      fdot = -(mu/(fp*r))*c1
      gdot = 1.0_rk - (mu/fp)*c2

      rtmp = r0*f + v0*g
      vtmp = r0*fdot + v0*gdot

      r0 = rtmp
      v0 = vtmp

    end if

  end if

  return
  end subroutine drift_dan
!!
  subroutine drift_one(mu, r, v, dtau, iflag)
  !--------------------------------------------------------------------------
  !         DRIFT_ONE.F90
  !--------------------------------------------------------------------------
  ! Drifts a body along its orbit using Danby alogorithm, and determines
  ! which variables to use.  If the accuracy is too poor, will redo the
  ! drift with a smaller step size.
  !
  ! Input:  mu      ==> Mass of central body
  !         r(NDIM) ==> Initial jacobi position of body
  !         v(NDIM) ==> Initial jacobi velocity of body
  !         dtau    ==> Time step
  !
  ! Output: r(NDIM) ==> Final jacobi position of body
  !         v(NDIM) ==> Final jacobi velocity of body
  !         iflag   ==> Status flag
  !                     - Zero if converged, non-zero otherwise
  !
  ! Authors: Hal Levison & Martin Duncan
  ! Date: 02/10/93
  ! Last revision: 02/10/93
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !--------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: mu, dtau

  ! Input and Output variables
  real(rk) :: r(NDIM), v(NDIM)

  ! Output variable
  integer(ik) :: iflag

  ! Internal variables
  integer(ik) :: i
  real(rk) :: dtau_temp

  !-----------------!
  ! Executable code !
  !-----------------!

  call drift_dan(mu, r, v, dtau, iflag)

  if(iflag /= 0) then

    dtau_temp = dtau/10.0_rk

    do i = 1, 10

      call drift_dan(mu, r, v, dtau_temp, iflag)
      if(iflag /= 0) exit ! Abandon all hope ye who enter here

    end do

  end if

  return
  end subroutine drift_one
!!
!!
end module drift
