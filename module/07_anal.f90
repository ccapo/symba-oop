module anal
! Module for routines in the anal directory
use swift
use util
use coord
use obl
implicit none

contains
!!
!!
  subroutine anal_energy(param, nbod, nbodm, pbod, ke, pe, energy, l)
  !---------------------------------------------------------------------------------------
  !         ANAL_ENERGY.F90
  !---------------------------------------------------------------------------------------
  ! Calculates the kinetic, potential and total energy of the system, along with the
  ! components of the total angular momentem of the system
  !
  ! N.B. We assume that the gravitational constant is G = 1.0
  !
  ! Input:  param      ==> Global parameters (See swift module & io_init_param)
  !         nbod       ==> Number of bodies
  !         nbodm      ==> Number of massive bodies
  !         pbod(nbod) ==> Database entries for bodies (See swift module)
  !
  ! Output: ke         ==> Kinetic energy
  !         pe         ==> Potential energy
  !         energy     ==> Total energy
  !         l(NDIM)    ==> Components of total angular momentum
  !
  ! Authors: Martin Duncan
  ! Date: ?
  ! Last revision: 01/24/97 HFL
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !---------------------------------------------------------------------------------------
  !$ use omp_lib
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm
  type(body_t), intent(in) :: pbod(nbod)
  type(param_t) :: param

  ! Output variables
  real(rk) :: ke, pe, energy, l(NDIM)

  ! Internal variables
  integer(ik) :: i, j, k
  real(rk) :: dr(NDIM), dr2, msys, pe_obl, rbi(NDIM), massi
  real(rk) :: lx, ly, lz !, kerr, perr, lxerr, lyerr, lzerr
  type(body_t) :: p(nbod), pi
  type(param_t) :: param_tmp

  !-----------------!
  ! Executable code !
  !-----------------!
  
  ! Make a copy of param and pbod so we don't alter the values stored in either
  param_tmp = param
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, pbod, p)
  do i = 1, nbod

    p(i) = pbod(i)

  end do
  !$OMP END PARALLEL DO

  ! Convert from heliocentric to barycentric coordinates
  call coord_h2b(param_tmp, nbod, p, msys)

  ! Initialize the kinetic, potential and total energy, along with the components of angular momentum
  ke = 0.5_rk*p(nbod)%mass*sum(p(nbod)%v**2)
  pe = 0.0_rk
  lx = p(nbod)%mass*(p(nbod)%r(2)*p(nbod)%v(3) - p(nbod)%r(3)*p(nbod)%v(2))
  ly = p(nbod)%mass*(p(nbod)%r(3)*p(nbod)%v(1) - p(nbod)%r(1)*p(nbod)%v(3))
  lz = p(nbod)%mass*(p(nbod)%r(1)*p(nbod)%v(2) - p(nbod)%r(2)*p(nbod)%v(1))

  ! Initialize the correction factors for the Kahan summation formula
  !kerr = 0.0_rk; perr = 0.0_rk
  !lxerr = 0.0_rk; lyerr = 0.0_rk; lzerr = 0.0_rk

  ! Compute the kinetic and components of angular momentum for the massive bodies, along with the potential energy
  ! amongst massive bodies and between massive and non-massive bodies
  do i = 1, nbodm

    ! Angular momentum components
    lx = lx + p(i)%mass*(p(i)%r(2)*p(i)%v(3) - p(i)%r(3)*p(i)%v(2))
    ly = ly + p(i)%mass*(p(i)%r(3)*p(i)%v(1) - p(i)%r(1)*p(i)%v(3))
    lz = lz + p(i)%mass*(p(i)%r(1)*p(i)%v(2) - p(i)%r(2)*p(i)%v(1))
    !lx = util_kahan_sum(lx, p(i)%mass*(p(i)%r(2)*p(i)%v(3) - p(i)%r(3)*p(i)%v(2)), lxerr)
    !ly = util_kahan_sum(ly, p(i)%mass*(p(i)%r(3)*p(i)%v(1) - p(i)%r(1)*p(i)%v(3)), lyerr)
    !lz = util_kahan_sum(lz, p(i)%mass*(p(i)%r(1)*p(i)%v(2) - p(i)%r(2)*p(i)%v(1)), lzerr)

    ! Kinetic energy
    ke = ke + 0.5_rk*p(i)%mass*sum(p(i)%v**2)
    !ke = util_kahan_sum(ke, 0.5_rk*p(i)%mass*sum(p(i)%v**2), kerr)

    ! Make a copy the information for body i available for each thread
    pi = p(i)

    ! Potential energy
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(j, dr, dr2) FIRSTPRIVATE(i, pi) &
    !$OMP SHARED(nbod, nbodm, p) REDUCTION(+ : pe)
    !FIRSTPRIVATE(i, pi, perr) LASTPRIVATE(perr)
    do j = i + 1, nbod

      if((pi%mass /= 0.0_rk) .and. (p(j)%mass /= 0.0_rk)) then

        dr = p(j)%r - pi%r
        dr2 = sum(dr**2)
        pe = pe - (pi%mass*p(j)%mass)/sqrt(dr2)
        !pe = util_kahan_sum(pe, -(pi%mass*p(j)%mass)/sqrt(dr2), perr)

      end if

    end do
    !$OMP END PARALLEL DO

  end do

  ! Compute the kinetic energy and components of angular momentum for the non-massive bodies
  !FIRSTPRIVATE(lxerr, lyerr, lzerr, kerr)
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) &
  !$OMP SHARED(nbod, nbodm, p) REDUCTION(+ : lx, ly, lz, ke)
  do i = nbodm + 1, nbod - 1

    ! Angular momentum components
    lx = lx + p(i)%mass*(p(i)%r(2)*p(i)%v(3) - p(i)%r(3)*p(i)%v(2))
    ly = ly + p(i)%mass*(p(i)%r(3)*p(i)%v(1) - p(i)%r(1)*p(i)%v(3))
    lz = lz + p(i)%mass*(p(i)%r(1)*p(i)%v(2) - p(i)%r(2)*p(i)%v(1))
    !lx = util_kahan_sum(lx, p(i)%mass*(p(i)%r(2)*p(i)%v(3) - p(i)%r(3)*p(i)%v(2)), lxerr)
    !ly = util_kahan_sum(ly, p(i)%mass*(p(i)%r(3)*p(i)%v(1) - p(i)%r(1)*p(i)%v(3)), lyerr)
    !lz = util_kahan_sum(lz, p(i)%mass*(p(i)%r(1)*p(i)%v(2) - p(i)%r(2)*p(i)%v(1)), lzerr)

    ! Kinetic energy
    ke = ke + 0.5_rk*p(i)%mass*sum(p(i)%v**2)
    !ke = util_kahan_sum(ke, 0.5_rk*p(i)%mass*sum(p(i)%v**2), kerr)

  end do
  !$OMP END PARALLEL DO

  ! Store components of angular momentum in a vector
  l = (/ lx, ly, lz /)

  ! If oblateness terms are needed, add to the potential
  if(param%loblate) then

    pe_obl = anal_energy_obl(param, nbod, p)
    pe = pe + pe_obl

  end if

  ! Total energy of the system
  energy = ke + pe

  return
  end subroutine anal_energy
!!
  function anal_energy_obl(param, nbod, pbod) result(oblpot)
  !--------------------------------------------------------------------------------------
  !         ANAL_ENERGY_OBL.F90
  !--------------------------------------------------------------------------------------
  ! Compute the oblateness contribution to the potential
  !
  ! Input:  param      ==> Global parameters (See swift module & io_init_param)
  !         nbod       ==> Number of bodies
  !         pbod(nbod) ==> Database entries for bodies (See swift module)
  !
  ! Output: oblpot     ==> Oblateness contribution to the potential
  !
  ! By: Chris Capobianco
  ! Date: 03/06/09 CCC - Written in Fortran 90/95 syntax
  !                    - Parallelized code using OpenMP
  !--------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  type(body_t) :: pbod(nbod)
  type(param_t) :: param

  ! Output variable
  real(rk) :: oblpot

  ! Internal variables
  real(rk) :: irh(nbod), ir3h(nbod)

  !-----------------!
  ! Executable code !
  !-----------------!

  !***********************************************************************!
  !*** Consider combining: util_ir_ir3 and obl_pot with the loop below ***!
  !***********************************************************************!

  ! Compute 1/r and 1/r^3
  call util_ir_ir3(nbod, 2, pbod, irh, ir3h)

  ! Compute the oblateness terms contribution to the potential
  call obl_pot(param, nbod, pbod, irh, oblpot)

  return
  end function anal_energy_obl
!!
!!
end module anal
