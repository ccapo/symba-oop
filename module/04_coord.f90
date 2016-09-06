module coord
! Module for the routines contained in coord
use swift
use util
implicit none

contains
!!
!!
  subroutine coord_b2h(param, nbod, pbod)
  !---------------------------------------------------------------------------------------
  !         COORD_B2H.F90
  !---------------------------------------------------------------------------------------
  ! Converts from barycentric to heliocentric coordinates
  !
  ! Input:  nbod       ==> Number of bodies
  !         pbod(nbod) ==> Current database entries for bodies (See swift module)
  !
  ! Output: pbod(nbod) ==> Updated database entries for bodies (See swift module)
  !
  ! REMARKS: Can of course use this to get coordinates relative to any body,
  !          simply by changing the coordinates being subtracted from all bodies
  !
  ! Authors: Martin Duncan
  ! WRITTEN: Jan 27/93
  ! REVISIONS: 02/17/95 HFL
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                         - Parallelized code using OpenMP
  !---------------------------------------------------------------------------------------
  implicit none

  ! Input variable
  integer(ik) :: nbod
  type(param_t) :: param

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i
  real(rk) :: rtmp(NDIM), vtmp(NDIM)

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Convert from barycentric => heliocentric
  if(param%lbary) then

    ! Store the barycentric position and velocity of the central body
    rtmp = pbod(1)%r
    vtmp = pbod(1)%v

    ! Compute the heliocentric position and velocity of all the bodies
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(rtmp, vtmp) &
    !$OMP SHARED(nbod, pbod)
    do i = 1, nbod

      pbod(i)%r = pbod(i)%r - rtmp
      pbod(i)%v = pbod(i)%v - vtmp

    end do
    !$OMP END PARALLEL DO

    ! Set coordinate flags
    param%lhelio = .true.
    param%lbary = .false.
    param%lcanon = .false.

  else

    write(*,'(a)') "SWIFT Warning:: coord_b2h"
    write(*,'(a)') " Coordinates already converted to heliocentric"

  end if

  return
  end subroutine coord_b2h
!!
  subroutine coord_h2b(param, nbod, pbod, msys)
  !---------------------------------------------------------------------------------------
  !         COORD_H2B.F90
  !---------------------------------------------------------------------------------------
  ! Converts from heliocentric to barycentric coordinates
  !
  ! Input:  nbod       ==> Number of bodies
  !         pbod(nbod) ==> Current database entries for bodies (See swift module)
  !
  ! Output: pbod(nbod) ==> Updated database entries for bodies (See swift module)
  !         msys       ==> Mass of the system
  !
  ! Authors: Martin Duncan
  ! WRITTEN: Jan 27/93
  ! REVISIONS: 02/22/94 HFL
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                         - Parallelized code using OpenMP
  !---------------------------------------------------------------------------------------
  implicit none

  ! Input variable
  integer(ik) :: nbod
  type(param_t) :: param

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Output variable
  real(rk) :: msys

  ! Internal variables
  integer(ik) :: i
  real(rk) :: x, y, z, vx, vy, vz
  !real(rk) :: msyserr, xerr, yerr, zerr, vxerr, vyerr, vzerr
  real(rk) :: rtmp(NDIM), vtmp(NDIM)

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Initialize the total system mass
  msys = pbod(1)%mass

  ! Convert from heliocentric => barycentric
  if(param%lhelio) then

    ! Initialize the temporary variables: x, y, z, vx, vy and vz
    msys = pbod(1)%mass
    x = 0.0_rk
    y = 0.0_rk
    z = 0.0_rk
    vx = 0.0_rk
    vy = 0.0_rk
    vz = 0.0_rk

    ! Initialize the correction factors for the Kahan floating-point summation formula
    !msyserr = 0.0_rk
    !xerr = 0.0_rk
    !yerr = 0.0_rk
    !zerr = 0.0_rk
    !vxerr = 0.0_rk
    !vyerr = 0.0_rk
    !vzerr = 0.0_rk

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, rtmp, vtmp) &
    !$OMP SHARED(nbod, pbod, msys, x, y, z, vx, vy, vz)
    !FIRSTPRIVATE(msyserr, xerr, yerr, zerr, vxerr, vyerr, vzerr)

    ! Compute the total mass of the system, along with the mass-weighted sum of positions and velocities
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+ : msys, x, y, z, vx, vy, vz)
    do i = 2, nbod

      msys = msys + pbod(i)%mass
      x = x + pbod(i)%mass*pbod(i)%r(1)
      y = y + pbod(i)%mass*pbod(i)%r(2)
      z = z + pbod(i)%mass*pbod(i)%r(3)
      vx = vx + pbod(i)%mass*pbod(i)%v(1)
      vy = vy + pbod(i)%mass*pbod(i)%v(2)
      vz = vx + pbod(i)%mass*pbod(i)%v(3)
      !msys = util_kahan_sum(msys, pbod(i)%mass, msyserr)
      !x = util_kahan_sum(x, pbod(i)%mass*pbod(i)%r(1), xerr)
      !y = util_kahan_sum(y, pbod(i)%mass*pbod(i)%r(2), yerr)
      !z = util_kahan_sum(z, pbod(i)%mass*pbod(i)%r(3), zerr)
      !vx = util_kahan_sum(vx, pbod(i)%mass*pbod(i)%v(1), vxerr)
      !vy = util_kahan_sum(vy, pbod(i)%mass*pbod(i)%v(2), vyerr)
      !vz = util_kahan_sum(vx, pbod(i)%mass*pbod(i)%v(3), vzerr)

    end do
    !$OMP END DO

    rtmp = [ x, y, z ]/msys
    vtmp = [ vx, vy, vz ]/msys

    !$OMP DO SCHEDULE(STATIC)
    do i = 1, nbod

      pbod(i)%r = pbod(i)%r - rtmp
      pbod(i)%v = pbod(i)%v - vtmp

    end do
    !$OMP END DO NOWAIT

    !$OMP END PARALLEL

    ! Set coordinate flags
    param%lhelio = .false.
    param%lbary = .true.
    param%lcanon = .false.

  else if(param%lcanon) then ! Convert from heliocentric => barycentric for positions *only*

    ! Initialize the temporary variables: x, y, z, vx, vy and vz
    msys = pbod(1)%mass
    x = 0.0_rk
    y = 0.0_rk
    z = 0.0_rk

    ! Initialize the correction factors for the Kahan floating-point summation formula
    !msyserr = 0.0_rk
    !xerr = 0.0_rk
    !yerr = 0.0_rk
    !zerr = 0.0_rk

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, rtmp) &
    !$OMP SHARED(nbod, pbod, msys, x, y, z)
    !FIRSTPRIVATE(msyserr, xerr, yerr, zerr)

    ! Compute the total mass of the system, along with the mass-weighted sum of positions and velocities
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+ : msys, x, y, z)
    do i = 2, nbod

      msys = msys + pbod(i)%mass
      x = x + pbod(i)%mass*pbod(i)%r(1)
      y = y + pbod(i)%mass*pbod(i)%r(2)
      z = z + pbod(i)%mass*pbod(i)%r(3)
      !msys = util_kahan_sum(msys, pbod(i)%mass, msyserr)
      !x = util_kahan_sum(x, pbod(i)%mass*pbod(i)%r(1), xerr)
      !y = util_kahan_sum(y, pbod(i)%mass*pbod(i)%r(2), yerr)
      !z = util_kahan_sum(z, pbod(i)%mass*pbod(i)%r(3), zerr)

    end do
    !$OMP END DO

    rtmp = [ x, y, z ]/msys

    !$OMP DO SCHEDULE(STATIC)
    do i = 1, nbod

      pbod(i)%r = pbod(i)%r - rtmp

    end do
    !$OMP END DO NOWAIT

    !$OMP END PARALLEL

    ! Set coordinate flags
    param%lhelio = .false.
    param%lbary = .true.
    param%lcanon = .false.

  else

    write(*,'(a)') "SWIFT Warning:: coord_h2b"
    write(*,'(a)') " Coordinates already converted to barycentric"

  end if

  return
  end subroutine coord_h2b
!!
  subroutine coord_vb2vh(param, nbod, pbod)
  !---------------------------------------------------------------------------------------
  !         COORD_VB2VH.F90
  !---------------------------------------------------------------------------------------
  ! Converts from barycentric to heliocentric coordinates for *velocities* only
  !
  ! Input:  nbod       ==> Number of bodies
  !         pbod(nbod) ==> Current database entries for bodies (See swift module)
  !
  ! Output: pbod(nbod) ==> Updated database entries for bodies (See swift module)
  !
  ! Authors: Hal Levison
  ! WRITTEN: 11/14/96
  ! REVISIONS: 11/21/96
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                         - Parallelized using OpenMP
  !---------------------------------------------------------------------------------------
  implicit none

  ! Input variable
  integer(ik) :: nbod
  type(param_t) :: param

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i
  real(rk) :: mstar, vx, vy, vz
  !real(rk) :: vxerr, vyerr, vzerr
  real(rk) :: vtmp(NDIM)

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Convert velocities from barycentric => heliocentric
  if(param%lcanon) then

    ! Initialize temporary variables
    mstar = pbod(1)%mass
    vx = -pbod(2)%mass*pbod(2)%v(1)
    vy = -pbod(2)%mass*pbod(2)%v(2)
    vz = -pbod(2)%mass*pbod(2)%v(3)

    ! Initialize the correction factors for the Kahan floating-point summation formula
    !vxerr = 0.0_rk
    !vyerr = 0.0_rk
    !vzerr = 0.0_rk

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, vtmp) FIRSTPRIVATE(mstar) &
    !$OMP SHARED(nbod, pbod, vx, vy, vz)
    !FIRSTPRIVATE(mstar, vxerr, vyerr, vzerr)

    !$OMP DO SCHEDULE(STATIC) REDUCTION(+ : vx, vy, vz)
    do i = 3, nbod

      vx = vx - pbod(i)%mass*pbod(i)%v(1)
      vy = vy - pbod(i)%mass*pbod(i)%v(2)
      vz = vz - pbod(i)%mass*pbod(i)%v(3)
      !vx = util_kahan_sum(vx, -pbod(i)%mass*pbod(i)%v(1), vxerr)
      !vy = util_kahan_sum(vy, -pbod(i)%mass*pbod(i)%v(2), vyerr)
      !vz = util_kahan_sum(vz, -pbod(i)%mass*pbod(i)%v(3), vzerr)

    end do
    !$OMP END DO

    vtmp = [ vx, vy, vz ]/mstar

    !$OMP DO SCHEDULE(STATIC)
    do i = 2, nbod

      pbod(i)%v = pbod(i)%v - vtmp

    end do
    !$OMP END DO NOWAIT

    !$OMP END PARALLEL

    ! Set heliocentric velocity of the central body
    pbod(1)%v = 0.0_rk

    ! Set coordinate flags
    param%lhelio = .true.
    param%lbary = .false.
    param%lcanon = .false.

  else

    write(*,'(a)') "SWIFT Warning:: coord_vb2vh"
    write(*,'(a)') " Velocities already converted to heliocentric"

  end if

  return
  end subroutine coord_vb2vh
!!
  subroutine coord_vh2vb(param, nbod, pbod)
  !----------------------------------------------------------------------------------------
  !         COORD_VH2VB.F90
  !----------------------------------------------------------------------------------------
  ! Converts from heliocentric to barycentric coordinates for *velocities* only
  !
  ! Input:  nbod       ==> Number of bodies
  !         pbod(nbod) ==> Current database entries for bodies (See swift module)
  !
  ! Output: pbod(nbod) ==> Updated database entries for bodies (See swift module)
  !
  ! Authors: Hal Levison
  ! WRITTEN: 11/14/96
  ! REVISIONS: 11/21/96
  !          : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                         - Parallelized code using OpenMP
  !----------------------------------------------------------------------------------------
  implicit none

  ! Input variable
  integer(ik) :: nbod
  type(param_t) :: param

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i
  real(rk) :: msys, vx, vy, vz
  !real(rk) :: msyserr, vxerr, vyerr, vzerr
  real(rk) :: vtmp(NDIM)

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Convert velocities from heliocentric => barycentric
  if(param%lhelio) then

    ! Initialize the total system mass, and the temporary variable: vx, vy and vz
    msys = pbod(1)%mass
    vx = 0.0_rk
    vy = 0.0_rk
    vz = 0.0_rk

    ! Initialize the correction factors for the Kahan floating-point summation formula
    !msyserr = 0.0_rk
    !vxerr = 0.0_rk
    !vyerr = 0.0_rk
    !vzerr = 0.0_rk

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, vtmp) &
    !$OMP SHARED(nbod, pbod, msys, vx, vy, vz)
    !FIRSTPRIVATE(msyserr, vxerr, vyerr, vzerr)

    !$OMP DO SCHEDULE(STATIC) REDUCTION(+ : msys, vx, vy, vz)
    do i = 2, nbod

      msys = msys + pbod(i)%mass
      vx = vx + pbod(i)%mass*pbod(i)%v(1)
      vy = vy + pbod(i)%mass*pbod(i)%v(2)
      vz = vz + pbod(i)%mass*pbod(i)%v(3)
      !msys = util_kahan_sum(msys, pbod(i)%mass, msyserr)
      !vx = util_kahan_sum(vx, pbod(i)%mass*pbod(i)%v(1), vxerr)
      !vy = util_kahan_sum(vy, pbod(i)%mass*pbod(i)%v(2), vyerr)
      !vz = util_kahan_sum(vz, pbod(i)%mass*pbod(i)%v(3), vzerr)

    end do
    !$OMP END DO

    vtmp = -1.0_rk*[ vx, vy, vz ]/msys

    !$OMP DO SCHEDULE(STATIC)
    do i = 2, nbod

      pbod(i)%v = pbod(i)%v + vtmp

    end do
    !$OMP END DO NOWAIT

    !$OMP END PARALLEL

    pbod(1)%v = -1.0_rk*[ vx, vy, vz ]/msys

    ! Set coordinate flags
    param%lhelio = .false.
    param%lbary = .false.
    param%lcanon = .true.

  else

    write(*,'(a)') "SWIFT Warning:: coord_vh2vb"
    write(*,'(a)') " Velocities already converted to barycentric"

  end if

  return
  end subroutine coord_vh2vb
!!
!!
end module coord
