module helio
! Module for routines in the helio directory
use swift
use util
use drift
implicit none

contains
!!
!!
  subroutine helio_lindrift(nbod, pbod, dtau)
  !--------------------------------------------------------------------------------------
  !         HELIO_LINDRIFT.F90
  !--------------------------------------------------------------------------------------
  ! This subroutine takes a linear drift due to mometum of the central body
  !
  ! Input:  nbod       ==> Number of massive bodies
  !         pbod(nbod) ==> Current database entries of bodies
  !         dtau       ==> Time step
  !
  ! Output: pbod(nbod) ==> Updated database entries of bodies
  !
  ! Remarks: Bases on Martin's code h2.f
  ! Authors: Hal Levison
  ! Date: 11/14/96
  ! Last revision: 01/08/97
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !--------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  real(rk) :: dtau

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i
  real(rk) :: mstar, px, py, pz
  !real(rk) :: pxerr, pyerr, pzerr
  real(rk) :: ptmp(NDIM)

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Store the inverse of the mass of the central body, and initialize the momentum
  mstar = pbod(1).mass
  px = 0.0_rk
  py = 0.0_rk
  pz = 0.0_rk

  ! Initialize the correction factor for the Kahan floating-point summation formula
  !pxerr = 0.0_rk
  !pyerr = 0.0_rk
  !pzerr = 0.0_rk

  !FIRSTPRIVATE(mstar, dtau, pxerr, pyerr, pzerr)
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, ptmp) FIRSTPRIVATE(mstar, dtau) &
  !$OMP SHARED(nbod, pbod, px, py, pz)

  !$OMP DO SCHEDULE(STATIC) REDUCTION(+ : px, py, pz)
  do i = 2, nbod

    px = px + pbod(i).mass*pbod(i).v(1)
    py = py + pbod(i).mass*pbod(i).v(2)
    pz = pz + pbod(i).mass*pbod(i).v(3)
    !px = util_kahan_sum(px, pbod(i).mass*pbod(i).v(1), pxerr)
    !py = util_kahan_sum(py, pbod(i).mass*pbod(i).v(2), pyerr)
    !pz = util_kahan_sum(pz, pbod(i).mass*pbod(i).v(3), pzerr)

  end do
  !$OMP END DO

  ! Make a private copy of the momentum vector
  ptmp = (/ px, py, pz /)/mstar

  !$OMP DO SCHEDULE(STATIC)
  do i = 2, nbod

    if(pbod(i).mass /= 0.0_rk) pbod(i).r = pbod(i).r + ptmp*dtau

  end do
  !$OMP END DO NOWAIT

  !$OMP END PARALLEL

  return
  end subroutine helio_lindrift
!!
  subroutine helio_drift(nbod, pbod, dtau)
  !---------------------------------------------------------------------------------------
  !         HELIO_DRIFT.F90
  !---------------------------------------------------------------------------------------
  ! Drifts all bodies along their orbit using Danby routine
  !
  ! Input:  nbod       ==> Number of bodies
  !         pbod(nbod) ==> Current database entries of bodies
  !         dtau       ==> Time step
  !
  ! Output: pbod(nbod) ==> Updated database entries of bodies
  !
  ! Remarks: Based on drift.f
  ! Authors: Hal Levison
  ! Date: 11/14/96
  ! Last revision: 01/08/97 HFL(?) - For SyMBA
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !---------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  real(rk) :: dtau

  ! Input and Output variables
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i, iflag(nbod)
  real(rk) :: mstar

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Store the mass of the central body, and set its flag (which should not be used)
  mstar = pbod(1).mass
  iflag(1) = 0

  ! Take a drift forward dtau
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(mstar, dtau) SHARED(nbod, pbod, iflag)
  do i = 2, nbod

    iflag(i) = 0
    if(pbod(i).mass /= 0.0_rk) call drift_one(mstar, pbod(i).r, pbod(i).v, dtau, iflag(i))

  end do
  !$OMP END PARALLEL DO

  ! If any drifts were flagged as unsuccessful, print them all
  if(any(iflag /= 0)) then

    do i = 2, nbod

      if(iflag(i) /= 0) then

        write(*,'(a)')              "SWIFT Error:: helio_drift"
        write(*,'(a,i7,a)')         "Particle:    ", i, " is lost!!!!!!!!!"
        write(*,'(a,2(1pe14.6))')   "Mass, dt:    ", mstar, dtau
        write(*,'(a,3(1pe14.6))')   "Helio. pos.: ", pbod(i).r
        write(*,'(a,3(1pe14.6),/)') "Bary. vel.:  ", pbod(i).v

      end if

    end do

    call util_exit(FAILURE)

  end if

  return
  end subroutine helio_drift
!!
  subroutine helio_kickvb(nbod, pbod, dtau)
  !---------------------------------------------------------------------------------------
  !         HELIO_KICKVB.F90
  !---------------------------------------------------------------------------------------
  ! Kicks the barcentric velocity of all bodies
  !
  ! Input:  nbod       ==> Number of bodies
  !         pbod(nbod) ==> Current database entries of bodies
  !         dtau       ==> Time step
  !
  ! Output: pbod(nbod) ==> Updated database entries of bodies
  !
  ! AUTHOR: M. Duncan
  ! DATE WRITTEN: Feb. 2, 1993
  ! REVISIONS: 02/18/93 HFL
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                         - Parallelized code using OpenMP
  !---------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  real(rk) :: dtau

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variable
  integer(ik) :: i

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Apply the kick to the velocities
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(dtau) SHARED(nbod, pbod)
  do i = 2, nbod

    pbod(i).v = pbod(i).v + pbod(i).a*dtau

  end do
  !$OMP END PARALLEL DO

  return
  end subroutine helio_kickvb
!!
!!
end module helio
