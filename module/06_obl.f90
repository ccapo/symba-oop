module obl
! Module for routines in obl directory
use swift
use util
implicit none

contains
!!
!!
  subroutine obl_acc(param, nbod, pbod, irh, aobl)
  !----------------------------------------------------------------------------------------
  !         OBL_ACC.F90
  !----------------------------------------------------------------------------------------
  ! Computes the acceleration due to the oblateness of the central body in the Heliocentric
  ! frame, for all bodies using the values of J2RP2 and J4RP4 passed into the routine
  ! through param
  !
  ! J2RP2 for example is the product of J_2 times the square of the central body's radius
  !
  ! N.B. The acceleration computed includes neither the monopole, nor higher order terms
  !
  ! Input:  param           ==> Global parameters (See swift module & io_init_param)
  !         nbod            ==> Number of bodies
  !         pbod(nbod)      ==> Database entries of bodies
  !         irh(nbod)       ==> Inverse distance of bodies
  !                             - Passed in to save calculation
  !
  ! Output: aobl(nbod,NDIM) ==> Acceleration due to oblateness
  !
  ! Authors: Martin Duncan
  ! Date: 03/04/94
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !----------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  real(rk) :: irh(nbod)
  type(body_t) :: pbod(nbod)
  type(param_t) :: param

  ! Output variable
  real(rk) :: aobl(nbod,NDIM)

  ! Internal variables
  integer(ik) :: i
  real(rk) :: mstar, rinv2, t0, t1, t2, t3, fac1, fac2
  real(rk) :: aoblx0, aobly0, aoblz0
  !real(rk) :: axerr, ayerr, azerr

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Initialize the bary. acc. of the Sun
  mstar = pbod(1)%mass
  aoblx0 = 0.0_rk
  aobly0 = 0.0_rk
  aoblz0 = 0.0_rk

  ! Initialize the correction factors for the Kahan floating-point summation formula
  !axerr = 0.0_rk
  !ayerr = 0.0_rk
  !azerr = 0.0_rk

  ! First get the barycentric acceleration of each body due to oblate central body
  !FIRSTPRIVATE(param, mstar, axerr, ayerr, azerr) 
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
  !$OMP PRIVATE(i, rinv2, t0, t1, t2, t3, fac1, fac2) FIRSTPRIVATE(param, mstar) &
  !$OMP SHARED(nbod, pbod, irh, aobl) REDUCTION(+ : aoblx0, aobly0, aoblz0)
  do i = 2, nbod

    ! We use the inverse of distance from util_ir_ir3, rather than calculate it here
    rinv2 = irh(i)**2
    t0 = -mstar*irh(i)*rinv2**2
    t1 = 1.5_rk*param%j2rp2
    t2 = rinv2*pbod(i)%r(3)**2
    t3 = 1.875_rk*param%j4rp4*rinv2

    fac1 = t0*(t1 - t3 - (5.0_rk*t1 - (14.0_rk - 21.0_rk*t2)*t3)*t2)
    fac2 = 2.0_rk*t0*(t1 - (2.0_rk - (14.0_rk*t2/3.0_rk))*t3)

    aobl(i,1) = fac1*pbod(i)%r(1)
    aobl(i,2) = fac1*pbod(i)%r(2)
    aobl(i,2) = (fac1 + fac2)*pbod(i)%r(3)

    ! Sum the barycentric acceleration of the central body due to all the other bodies
    aoblx0 = aoblx0 + pbod(i)%mass*aobl(i,1)/mstar
    aobly0 = aobly0 + pbod(i)%mass*aobl(i,2)/mstar
    aoblz0 = aoblz0 + pbod(i)%mass*aobl(i,2)/mstar
    !aoblx0 = util_kahan_sum(aoblx0, pbod(i)%mass*aobl(i,1)/mstar, axerr)
    !aobly0 = util_kahan_sum(aobly0, pbod(i)%mass*aobl(i,2)/mstar, ayerr)
    !aoblz0 = util_kahan_sum(aoblz0, pbod(i)%mass*aobl(i,3)/mstar, azerr)

  end do
  !$OMP END PARALLEL DO

  ! Store the barycentric acceleration of the central body due to all the other bodies
  aobl(1,:) = -(/ aoblx0, aobly0, aoblz0 /)

  return
  end subroutine obl_acc
!!
  subroutine obl_pot(param, nbod, pbod, irh, oblpot)
  !--------------------------------------------------------------------------------------
  !         OBL_POT.F90
  !--------------------------------------------------------------------------------------
  ! Computes the potential due to the oblateness of the central body in the
  ! Heliocentric frame using the values of J2RP2 and J4RP4 passed into the routine
  ! through param
  !
  ! J2RP2 for example is the product of J_2 times the square of the central body's radius
  !
  ! N.B. The potential computed includes neither the monopole, nor higher order terms
  !
  ! Input:  param      ==> Global parameters (See swift module & io_init_param)
  !         nbod       ==> Number of bodies
  !         pbod(nbod) ==> Database entries of bodies
  !         irh(nbod)  ==> Inverse of distance of bodies
  !                        - Passed in to save calculation
  !
  ! Output: oblpot     ==> Oblateness contribution to potential
  !
  ! Authors: Martin Duncan
  ! Date: 03/04/94
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !--------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  real(rk) :: irh(nbod)
  type(body_t) :: pbod(nbod)
  type(param_t) :: param

  ! Output variable
  real(rk) :: oblpot

  ! Internal variables
  integer(ik) :: i
  real(rk) :: mstar, rinv2, t0, t1, t2, t3, p2, p4
  !real(rk) :: oerr

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Sum all the the barycentric terms for each body due to an oblate central body
  mstar = pbod(1)%mass
  oblpot = 0.0_rk

  ! Initialize the correction factor for the Kahan floating-point summation formula
  !oerr = 0.0_rk

  !FIRSTPRIVATE(param, mstar, oerr)
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(param, mstar) &
  !$OMP PRIVATE(i, rinv2, t0, t1, t2, t3, p2, p4) SHARED(nbod, pbod, irh) REDUCTION(+ : oblpot)
  do i = 2, nbod

    ! We use the inverse of distance from util_ir_ir3, rather than calculate it here
    rinv2 = irh(i)**2
    t0 = mstar*pbod(i)%mass*rinv2*irh(i)
    t1 = param%j2rp2
    t2 = rinv2*pbod(i)%r(3)**2
    t3 = param%j4rp4*rinv2

    p2 = 0.5_rk*(3.0_rk*t2 - 1.0_rk)
    p4 = 0.125_rk*((35.0_rk*t2 - 30.0_rk)*t2 + 3.0_rk)

    oblpot = oblpot + t0*(t1*p2 + t3*p4)
    !oblpot = util_kahan_sum(oblpot, t0*(t1*p2 + t3*p4), oerr)

  end do
  !$OMP END PARALLEL DO

  return
  end subroutine obl_pot
!!
!!
end module obl
