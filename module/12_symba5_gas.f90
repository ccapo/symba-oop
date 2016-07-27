module symba5_gas
! Module for routines in symba5_gas directory
use swift
use orbel
use symba5
implicit none

contains
!!
!!
  function symba5_gas_drag_coefficient(ma, kn) result(cd)
  !-----------------------------------------------------------------------
  !         SYMBA5_GAS_DRAG_COEFFICIENT.F90
  !-----------------------------------------------------------------------
  !
  ! Computes the drag coefficient (C_D) as a function of Knudsen number
  ! (Kn) and Mach number (Ma), via the Reynolds number (Re)
  !
  ! Input:  ma ==> mach number
  !         kn ==> knudsen number
  !
  ! Output: cd ==> drag coefficient
  !
  ! Author: Martin Duncan
  ! Date: 07/11/07
  !     : 02/09/09 CCC - Coverted to Fortran 90/95 Syntax
  !-----------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: ma, kn

  ! Output variable
  real(rk) :: cd

  ! Internal variable
  real(rk) :: re

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Reynolds number
  re = 4.448823857787100274_rk*ma/kn

  if(kn < 1.0_rk) then

    if(ma >= 1.0_rk) then

      cd = 2.0_rk

    else

      if(re > 1.0e3_rk) then

        cd = 0.44_rk + 1.56_rk*ma**2

      else

        cd = 2.0_rk*ma**2 + 24.0_rk*(1.0_rk - ma**2)*(1.0_rk + 0.15_rk*re**0.687_rk)/re

      end if

    end if

  else

    if(ma < 1.8_rk) then

      cd = 3.6_rk/ma

    else

      cd = 2.0_rk

    end if

  endif

  return
  end function symba5_gas_drag_coefficient
!!
  subroutine symba5_gas_drag(param, time, dtau, nbod, nbodm, pbod, lhill)
  !---------------------------------------------------------------------------------
  !         SYMBA5_GAS_DRAG.F90
  !---------------------------------------------------------------------------------
  ! Computes the acceleration due to the aerodynamic gas drag for the
  ! non-gravitating bodies
  !
  ! Input:  param       ==> Global parameters (See swift module & io_init_param)
  !         time        ==> Current time
  !         dtau        ==> Timestep
  !         nbod        ==> Number of bodies
  !         nbodm       ==> Index of last massive body
  !         pbod(nbod)  ==> Current database entries of bodies
  !         lhill(nbod) ==> Hill sphere flag of bodies
  !
  ! Output: pbod(nbod)  ==> Updated database entreis of bodies
  !
  ! By: Martin Duncan & Hal Levison
  ! Date: ?
  !     : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                    - Parallelized code using OpenMP
  !---------------------------------------------------------------------------------
  implicit none

  ! Input/Output variable
  integer(ik) :: nbod, nbodm
  logical(lk) :: lhill(nbod)
  real(rk) :: time, dtau
  type(body_t) :: pbod(nbod)
  type(param_t) :: param

  ! Internal variable
  integer(ik) :: j
  real(rk) :: exp_factor, s, r, z_0, rhogas, eta, vkep, vfac
  real(rk) :: vrelabs, mach, knudsen, c_d, gdrag, zarg, mstar
  real(rk) :: vgas(NDIM), vrel(NDIM), agas(NDIM)

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Exponential decay factor
  exp_factor = exp(-time/param%taugas)

  ! Mass of central body
  mstar = pbod(1)%mass

  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(exp_factor, mstar, dtau, param) &
  !$OMP PRIVATE(j, s, r, z_0, zarg, rhogas, eta, vkep, vfac, vgas, vrel, vrelabs, mach, knudsen, c_d, gdrag, agas) &
  !$OMP SHARED(nbod, nbodm, pbod, lhill)
  do j = nbodm + 1, nbod

    ! Projected distance in the xy-plane
    s = sqrt(sum(pbod(j)%r(1:2)**2))

    ! If body j is within the hill sphere of a massive body, or its projected distance
    ! in the xy-plane is outside the gas disk region
    if(lhill(j) .or. (s < param%rgi) .or. (s > param%rgf)) then

      ! Set the acceleration due to gas drag to zero
      agas = 0.0_rk

    else

      ! Scale height of the gas
      z_0 = param%zscale*s**1.25_rk ! scale height (Hayashi et al., 1980)

      ! Density follows power law model: s^(-gpower)*exp(-(z/z0)^2)
      zarg = pbod(j)%r(3)/z_0
      rhogas = param%rhogas0*exp(-zarg**2)*s**(-param%gpower)

      ! Pressure of gas parameter eta (e.g. Kokubo and Ida)
      eta = 6.0e-4_rk*(param%gpower + 0.5_rk)*sqrt(s)

      ! Distance and keplerian velocity for body j
      r = sqrt(sum(pbod(j)%r**2))
      vkep = sqrt(mstar/r)

      ! Gas moves slower than keplerian velocity
      vfac = vkep*sqrt(1.0_rk - 2.0_rk*eta)            ! Adachi et al., 1976
      vgas(1) = -vfac*pbod(j)%r(2)/s                   ! Brakes vx
      vgas(2) = vfac*pbod(j)%r(1)/s                    ! Brakes vy
      vgas(3) = 0.0_rk                                 ! No braking in vz

      vrel = pbod(j)%v - vgas                          ! Relative velocity between body j and gas
      vrelabs = sqrt(sum(vrel**2))                     ! Magnitude of relative velocity between body j and gas

      mach = 3.32126045_rk*vrelabs*s**0.25_rk          ! Mach number
      knudsen = 1.1097308e-6_rk/(rhogas*pbod(j)%rdrag) ! Knudsen number

      ! Compute the drag coefficient using the Mach and Knudsen number
      c_d = symba5_gas_drag_coefficient(mach, knudsen)

      ! Compute the inverse of aerodynamic gas drag timescale
      gdrag = 0.844260701751676972_rk*c_d*vrelabs*rhogas*exp_factor/(pbod(j)%rdrag*pbod(j)%rho)

      ! Compute the acceleration due to aerodynamic gas drag
      agas = -gdrag*vrel

      ! Apply the aerodynamic gas drag acceleration
      pbod(j)%v = pbod(j)%v + agas*dtau

    end if

  end do
  !$OMP END PARALLEL DO

  return
  end subroutine symba5_gas_drag
!!
  subroutine symba5_gas_typeI(param, time, dtau, nbod, nbodm, pbod)
  !-------------------------------------------------------------------------------------
  !         SYMBA5_GAS_TYPEI.F90
  !-------------------------------------------------------------------------------------
  ! Computes the acceleration due to Type - I migration for massive gravitating bodies
  !
  ! Input:  param      ==> Global parameters (See swift module & io_init_param)
  !         time       ==> Current time
  !         dtau       ==> Timestep
  !         nbod       ==> Number of bodies
  !         nbodm      ==> Index of last massive body
  !         pbod(nbod) ==> Current database entries of bodies
  !
  ! Output: pbod(nbod) ==> Updated database entreis of bodies
  !
  ! By: Martin Duncan & Hal Levison
  ! Date: ?
  !     : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------------------
  implicit none

  ! Input/Output variables
  integer(ik) :: nbod, nbodm
  real(rk) :: time, dtau
  type(body_t) :: pbod(nbod)
  type(param_t) :: param

  ! Internal variables
  integer(ik) :: i, ialpha
  real(rk) :: s2, r2, gm, apl, epl, qpl, sigma
  real(rk) :: hovera, omega_min, ta, te, ti, eratio, mratio
  real(rk) :: vdotr, faca, mdisk, fgap, atypeI(NDIM)

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Exponential decay of surface gas density
  sigma = param%sigma0*exp(-time/param%taugas)

  do i = 2, nbodm

    ! Projected distance squared in the xy-plane
    s2 = sum(pbod(i)%r(1:2)**2)

    ! If body i is outside the gas disk region
    if((s2 < param%rgi**2) .or. (s2 > param%rgf**2)) then

      ! Set the acceleration due to Type-I to zero
      atypeI = 0.0_rk

    else

      ! Get the keplerian orbital elements for body i
      gm = pbod(1)%mass + pbod(i)%mass
      call orbel_xv2aeq(pbod(i)%r, pbod(i)%v, gm, ialpha, apl, epl, qpl)

      ! Compute the semi-major axis and eccentricity damping timescales
      mratio = pbod(i)%mass/pbod(1)%mass
      hovera = param%zscale*apl**0.25_rk
      omega_min = sqrt(apl**3/gm)
      mdisk = PI*sigma*(apl**(2.0_rk - param%spower))/pbod(1)%mass
      ta = (omega_min*hovera**2)/(param%ca*mratio*mdisk)
      te = (omega_min*hovera**4)/(param%ce*mratio*mdisk)

      ! Now multiply by the eccentricity correction from Papaloizou & Larwood (2000)
      if(pbod(i)%lgap) then

        fgap = (apl - pbod(i)%rgap)/(abs(apl - pbod(i)%rgap) + apl*pbod(i)%wgap)

      else

        eratio = epl/hovera
        ta = ta*(1.0_rk + (eratio/1.3_rk)**5)/(1.0_rk - (eratio/1.1_rk)**4)
        te = te*(1.0_rk + 0.25_rk*eratio**3)
        fgap = 1.0_rk

      end if

      ! For *now* assume timescale for vertical damping ti equals te even though
      ! Papaloizou and Larwood don't discuss inclination dependent term in te.
      ! What's done here could generate anomolous inclination damping timescale
      ! if inc takes body above scale height and ecc remains less than hovera.
      ! Can we legitimately use parameter like sin(inc)/hovera to modify ti as
      ! we do for te with e/hovera?
      ti = te

      ! Now get the acceleration components as in Papaloizou & Larwood (2000)
      vdotr = dot_product(pbod(i)%r, pbod(i)%v)
      r2 = sum(pbod(i)%r**2)
      faca = 2.0_rk*vdotr/(r2*te)
      !faca = fgap*2%d0*vdotr/(r2*te)  ! ecc and inc damping timescales -> 0 when in a gap

      atypeI(1) = -pbod(i)%v(1)*fgap/ta - faca*pbod(i)%r(1)
      atypeI(2) = -pbod(i)%v(2)*fgap/ta - faca*pbod(i)%r(2)
      atypeI(3) = -pbod(i)%v(3)*fgap/ta - faca*pbod(i)%r(3) - 2.0_rk*pbod(i)%v(3)/ti

      ! Apply the Type-I acceleration
      pbod(i)%v = pbod(i)%v + atypeI*dtau

    end if

  end do

  return
  end subroutine symba5_gas_typeI
!!
  subroutine symba5_gas_drag_kick(param, time, dtau, nbod, nbodm, pbod)
  !-----------------------------------------------------------------------------------
  !         SYMBA5_GAS_DRAG_KICK.F90
  !-----------------------------------------------------------------------------------
  ! Computes the acceleration due to aerodynamic gas drag and Type - I migration,
  ! then applies the kick to the velocities of the appropriate bodies.
  !
  ! Input:  param      ==> Global parameters (See swift module & io_init_param)
  !         time       ==> Current time
  !         dtau       ==> Time step
  !         nbod       ==> Number of bodies
  !         nbodm      ==> Index of last massive body
  !         pbod(nbod) ==> Current database entries of bodies
  !
  ! Output: pbod(nbod) ==> Updated database entreis of bodies
  !
  ! By: Hal Levison & Martin Duncan
  ! Date: ?
  !     : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                    - Parallelized code using OpenMP
  !-----------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm
  real(rk) :: time, dtau
  type(param_t) :: param

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  logical(lk) :: lhill(nbod)
  integer(ik) :: i, j, ialpha
  real(rk) :: gm, apl, epl, qpl, dr(NDIM), dr2
  type(body_t) :: pbodi

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Initialize the Hill radius flag
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i) SHARED(nbod, lhill)
  do i = 1, nbod

    lhill(i) = .false.

  end do
  !$OMP END PARALLEL DO

  ! Check if planetesimals are within the Hill sphere *any* massive body
  do i = 2, nbodm

    ! First update the location of any gaps
    if(pbod(i)%lgap) then

      gm = pbod(1)%mass + pbod(i)%mass
      call orbel_xv2aeq(pbod(i)%r, pbod(i)%v, gm, ialpha, apl, epl, qpl)
      pbod(i)%rgap = apl

    end if

    ! Make a copy of the information for body i available for each thread
    pbodi = pbod(i)

    ! Set Hill radius flag for a planetesimal if within the Hill sphere of *any* massive body
    ! This will turn off the aerodynamic gas drag calculation in the subsequent routine
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(j, dr, dr2) &
    !$OMP FIRSTPRIVATE(pbodi) SHARED(nbod, nbodm, pbod, lhill)
    do j = nbodm + 1, nbod

      dr = pbod(j)%r - pbodi%r
      dr2 = sum(dr**2)
      if(dr2 < pbodi%rhill**2) lhill(j) = .true.

    end do
    !$OMP END PARALLEL DO

  end do

  ! Compute and apply the acceleration due aerodynamic gas drag for planetesimals
  call symba5_gas_drag(param, time, dtau, nbod, nbodm, pbod, lhill)

  ! Compute and apply the acceleration due to Type-I migration for massive bodies
  call symba5_gas_typeI(param, time, dtau, nbod, nbodm, pbod)

  return
  end subroutine symba5_gas_drag_kick
!!
  subroutine symba5_gas_step(param, time, nbod, nbodm, pbod, pmerge, eoffset)
  !---------------------------------------------------------------------------------------------------
  !         SYMBA5_GAS_STEP.F90
  !---------------------------------------------------------------------------------------------------
  ! Takes a step in heliocentric frame by doing a KICK, a DRIFT than a KICK, but also handles
  ! close encounters and mergers.  This version also incorporates analytical prescriptions of
  ! aerodynamic gas drag and Type - I migration.
  !
  ! Input:  param    ==> Global parameters (See swift module & io_init_param)
  !         time       ==> Current time
  !         nbod       ==> Number of bodies
  !         nbodm      ==> Index of the last massive body
  !         pbod(nbod) ==> Current database entries of bodies
  !         eoffset    ==> Current energy offset
  !
  ! Output: pbod(nbod) ==> Updated database entries of bodies (if a merger occurred)
  !         pmerge     ==> Merger list (if a merger occurred)
  !         eoffset    ==> Updated energy offset (if a merger occurred)
  !
  ! Remarks: Based on symba2_step_pl.f
  ! Authors: Hal Levison
  ! Date: 11/27/97
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !---------------------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm
  real(rk) :: time
  type(param_t) :: param

  ! Input and Output variables
  real(rk) :: eoffset
  type(body_t) :: pbod(nbod)

  ! Output variable
  integer(ik), pointer :: pmerge(:)

  ! Internal variables
  logical(lk) :: lenc
  integer(ik) :: i
  real(rk) :: dth
  type(enc_t) :: penc(nbodm)

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Half time step
  dth = 0.5_rk*param%dt
  
  ! Initialize merger list
  if(associated(pmerge)) deallocate(pmerge)

  ! Construct encounter list
  if(param%lclose) call symba5_encounter_list(param, nbod, nbodm, pbod, penc, lenc)
  !lenc = .false.

  ! Apply gas drag and Type-I formulae
  if(param%lgas) call symba5_gas_drag_kick(param, time, dth, nbod, nbodm, pbod)

  ! Do a step
  if(lenc) then

    ! If there are encounters, perform a recursive step
    call symba5_step_interp(param, time, nbod, nbodm, pbod, penc, pmerge, eoffset)

  else

    ! If there are no encounters, perform a simple step
    call symba5_step_helio(param, nbod, nbodm, pbod)

  end if

  ! Apply gas drag and Type-I formulae
  if(param%lgas) call symba5_gas_drag_kick(param, time, dth, nbod, nbodm, pbod)

  ! Print number of encounters and mergers found in this time step
  !write(100,*) time, sum(penc%nenc), size(pmerge)

  ! Deallocate encounter list
  do i = 2, nbodm

    if(associated(penc(i)%ienc)) then

      penc(i)%nenc = 0
      deallocate(penc(i)%ienc)

    end if

  end do

  return
  end subroutine symba5_gas_step
!!
!!
end module symba5_gas
