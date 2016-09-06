module symba5
! Module for routines in symba5 directory
!
! Remarks:  Copied from symba5.inc
! Author:  Hal Levison
! Date:    3/20/97
! Last revision:
use swift
use orbel
use util
use coord
use obl
use drift
use helio
use anal
use discard
implicit none

! Maximum number of encounters
integer(ik), parameter :: NENMAX = 131072 ! must be less than 2**17 = 131072

! Scale factor for Hill's sphere to take shorter time step
real(rk), parameter :: RHSCALE = 6.5_rk

! Ratio of the number of time steps in the adjoining shells
integer(ik), parameter :: NTENC = 3

! Ratio of shell radii squared
!real(rk), parameter :: RSHELL = 0.48075_rk ! rshell ~ ntenc^(-2/3)
!real(rk), parameter :: RSHELL = real(NTENC, rk)**(-2.0_rk/3.0_rk) ! rshell ~ ntenc^(-2/3)
real(rk), parameter :: RSHELL = 0.480749856769136133_rk ! rshell ~ ntenc^(-2/3)

! Common logical flag to control force calculation in symb5_step_helio
! Variable is made private, so it can't accidently be changed
! when loading this module, only these routines can change it.
logical, private, save :: lfirst_force = .true.

contains
!!
!!
  subroutine symba5_check(irec, dtau, pbodi, pbodj, lvdotr, lenc)
  !---------------------------------------------------------------------------------------
  !         symba5_check.F90
  !---------------------------------------------------------------------------------------
  ! Checks if the two bodies are in an encounter region
  !
  ! Input:  irec   ==> Current recursion level
  !         dtau   ==> Time step
  !         pbodi  ==> Heliocentric position of both bodies
  !         pbodj  ==> Heliocentric velocity of both bodies
  !
  ! Output: lvdotr ==> Approach flag
  !                    - lvdotr = .TRUE. if bodies are approaching
  !                    - lvdotr = .FALSE. if bodies are receding
  !         lenc   ==> Encounter flag
  !                    - lenc = .TRUE. if bodies are in an encounter region
  !                    - lenc = .FALSE. if bodies are not in an encounter region
  !
  ! Remarks: Based on plh_chk.f, same as symba_chk.f
  ! Authors: Hal Levison
  ! Date: 03/20/97
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Now only handles *two* bodies at a time
  !---------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: irec
  real(rk) :: dtau
  type(body_t) :: pbodi, pbodj

  ! Output variables
  logical(lk) :: lvdotr, lenc

  ! Internal variables
  real(rk) :: dr(NDIM), dv(NDIM), dr2, dv2, vdotr
  real(rk) :: r2crit, tmin, r2min

  !-----------------!
  ! Executable code !
  !-----------------!

  ! First check if we're already in the encounter region.
  r2crit = ((pbodi%rhill + pbodj%rhill)*RHSCALE*RSHELL**irec)**2

  dr = pbodj%r - pbodi%r
  dr2 = sum(dr**2)

  dv = pbodj%v - pbodi%v
  dv2 = sum(dv**2)

  vdotr = dot_product(dr, dv)
  lvdotr = vdotr < 0.0_rk

  tmin = -vdotr/dv2

  ! If we're heading outward, use dr2 to determine if we are in an encounter region
  if(vdotr > 0.0_rk) then

    if(dr2 >= r2crit) then

      lenc = .false.

    else

      lenc = .true.

    end if

  else

    ! We are converging, so we need to calculate the minimum separation attained in time dtau
    if(tmin < dtau) then

      r2min = dr2 - (vdotr**2)/dv2

    else

      r2min = dr2 + 2.0_rk*vdotr*dtau + dv2*dtau**2

    end if

    r2min = min(r2min, dr2) ! Really make sure

    if(r2min <= r2crit) then

      lenc = .true.

    else

      lenc = .false.

    end if

  end if

  return
  end subroutine symba5_check
!!
  subroutine symba5_encounter_list(param, nbod, nbodm, pbod, penc, lenc)
  !-------------------------------------------------------------------------------------------------
  !         SYMBA5_ENCOUNTER_LIST.F90
  !-------------------------------------------------------------------------------------------------
  ! Searches for potential encounters between the massive gravitating bodies, and non-gravitating
  ! bodies
  !
  ! Input:  param       ==> Global parameters (See swift module & io_init_param)
  !         nbod        ==> Number of bodies
  !         nbodm       ==> Index of last massive gravitating body
  !         pbod(nbod)  ==> Current database entries of bodies
  !
  ! Output: pbod(nbod)  ==> Updated database entries of bodies
  !         penc(nbodm) ==> Encounter list for massive gravitating bodies
  !         lenc        ==> Encounter flag
  !                         - lenc = .TRUE. if any encounters are found
  !                         - lenc = .FALSE. otherwise
  !
  ! By: Chris Capobianco
  ! Date: 10/10/08
  !-------------------------------------------------------------------------------------------------
  !$ use omp_lib
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm
  type(param_t) :: param

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Output variables
  logical(lk) :: lenc
  type(enc_t) :: penc(nbodm)

  ! Internal variables
  logical(lk) :: lvdotr, lenc_local
  integer(ik) :: i, j, k, id, ierr, jstart
  integer(ik) :: istart(nthreads), nenc(nthreads), penc_local(nbod)
  real(rk) :: dr(NDIM), dv(NDIM)
  type(body_t) :: pbodi

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Initialize the recursion level, max recursion level and encounter and merger flag for all bodies
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i) SHARED(nbod, pbod)
  do i = 1, nbod

    pbod(i)%lenc = .false.
    pbod(i)%rlevel = -1
    pbod(i)%rlevelmax = 0

  end do
  !$OMP END PARALLEL DO

  ! Initialize the number of encounters for each massive body
  penc%nenc = 0

  ! Initialize global encounter flag
  lenc = .false.

  ! Check for encounters
  do i = 2, nbodm

    ! Extract the information for body i, and make a copy available for each thread
    pbodi = pbod(i)

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(j, id, lvdotr, lenc_local) FIRSTPRIVATE(i, pbodi, param) &
    !$OMP SHARED(nbod, nbodm, nthreads, pbod, istart, nenc, penc_local)

    id = 1                                  ! Set thread identifier for *serial* case
    !$ id = omp_get_thread_num() + 1        ! Set thread identifier for *parallel* case
    nenc(id) = 0                            ! Initialize encounter counter for each thread
    istart(id) = (id - 1)*nbod/nthreads + 1 ! Set the starting position in penc_local & lvdotr_local for each thread

    !$OMP DO SCHEDULE(STATIC)
    do j = i + 1, nbod

      ! Check if we have an encounter
      call symba5_check(0, param%dt, pbodi, pbod(j), lvdotr, lenc_local)

      ! If we have an encounter, store relevant information
      if(lenc_local) then

        ! Encounter detected
        pbod(j)%lenc = .true.                          ! Set the encounter flag for body j
        pbod(j)%rlevel = 0                             ! Set the recursion level to 0 for body j
        penc_local(istart(id) + nenc(id)) = pbod(j)%id ! Store identifier of body j for this thread
        nenc(id) = nenc(id) + 1                        ! Increment the encounter counter for this thread

      end if

    end do
    !$OMP END DO NOWAIT

    !$OMP END PARALLEL

    ! Sum the number of encounters for body i
    penc(i)%nenc = sum(nenc)

    ! If there any encounters, append to the global encounter list for body i
    if(penc(i)%nenc > 0) then

      ! Set global encounter flag
      lenc = .true.

      ! Set the encounter flag, and the recursion level to 0 for body i
      pbod(i)%lenc = .true.
      pbod(i)%rlevel = 0

      ! Define the size of the global encounter list for body i
      allocate(penc(i)%ienc(penc(i)%nenc), penc(i)%lvdotr(penc(i)%nenc))

      ! Initialize starting location in global encounter list
      jstart = 0

      ! Search through each thread's local encounter list
      do k = 1, nthreads

        ! If an encounter was found for thread k, append nenc(k) entries to penc(i)%ienc and penc(i)%lvdotr
        if(nenc(k) > 0) then

          do j = 1, nenc(k)

            penc(i)%ienc(jstart + j) = penc_local(istart(k) + j - 1)

          end do

          ! Update the starting position in the global encounter list
          jstart = jstart + nenc(k)

        end if

      end do

    end if

  end do

  return
  end subroutine symba5_encounter_list
!!
  subroutine symba5_helio_drift(irec, dtau, nbod, pbod)
  !--------------------------------------------------------------------------
  !         SYMBA5_HELIO_DRIFT.F90
  !--------------------------------------------------------------------------
  ! Drifts all bodies along their orbit using Danby routine, provided that
  ! their recursion level matches the current recursion level
  !
  ! Input:  irec       ==> Current recursion level
  !         dtau       ==> Time step
  !         nbod       ==> Number of bodies
  !         pbod(nbod) ==> Current database entries of bodies
  !
  ! Output: pbod(nbod) ==> Updated database entries of bodies
  !
  ! Remarks: Based on helio_drift.f
  ! Authors: Hal Levison
  ! Date: 01/20/97
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !--------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, irec
  real(rk) :: dtau

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i, iflag(nbod)
  real(rk) :: mstar

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Make a copy of the central body's mass, and set its flag (which should not be used)
  mstar = pbod(1)%mass
  iflag(1) = 0

  ! Take a drift forward dtau if the recursion level of the body matches current recursion level
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(mstar, dtau, irec) PRIVATE(i) &
  !$OMP SHARED(nbod, pbod, iflag)
  do i = 2, nbod

    iflag(i) = 0
    if((pbod(i)%rlevel == irec) .and. (pbod(i)%mass /= 0.0_rk)) call drift_one(mstar, pbod(i)%r, pbod(i)%v, dtau, iflag(i))

  end do
  !$OMP END PARALLEL DO

  ! If any drifts were unsuccessful, print all of them, then abort
  if(any(iflag /= 0)) then

    do i = 2, nbod

      if(iflag(i) /= 0) then

        write(*,'(a)')              "SWIFT Error:: symba5_helio_drift"
        write(*,'(a,i6,a)')         " Particle:    ", i, " is lost!!!!!!!!!"
        write(*,'(a,2(1pe14.6))')   " Mstar, dtau:  ", mstar, dtau
        write(*,'(a,3(1pe14.6))')   " Helio. pos.: ", pbod(i)%r
        write(*,'(a,3(1pe14.6),/)') " Bary. vel.:  ", pbod(i)%v

      end if

    end do

    call util_exit(FAILURE)

  end if

  return
  end subroutine symba5_helio_drift
!!
  subroutine symba5_enc_drift(irec, dtau, nbod, nbodm, pbod, penc)
  !--------------------------------------------------------------------------
  !         SYMBA5_HELIO_DRIFT.F90
  !--------------------------------------------------------------------------
  ! Drifts all bodies along their orbit using Danby routine, provided that
  ! their recursion level matches the current recursion level
  !
  ! Input:  irec       ==> Current recursion level
  !         dtau       ==> Time step
  !         nbod       ==> Number of bodies
  !         pbod(nbod) ==> Current database entries of bodies
  !
  ! Output: pbod(nbod) ==> Updated database entries of bodies
  !
  ! Remarks: Based on helio_drift.f
  ! Authors: Hal Levison
  ! Date: 01/20/97
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !--------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm, irec
  real(rk) :: dtau
  type(enc_t) :: penc(nbodm)

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  logical(lk) :: lflag, ldflag(nbod)
  integer(ik) :: i, j, k, iflag
  real(rk) :: mstar

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Initialize error index and drift flag for each body
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) SHARED(nbod, ldflag)
  do i = 1, nbod

    ldflag(i) = .true.

  end do
  !$OMP END PARALLEL DO

  ! Make a copy of the central body's mass, and set its flag (which should not be used)
  mstar = pbod(1)%mass

  ! Take a drift forward dtau if the recursion level of the body matches current recursion level
  do i = 2, nbodm

    ! If there are any bodies in the encounter region around body i
    if(penc(i)%nenc > 0) then

      ! Perform a drift for body i
      if((pbod(i)%rlevel == irec) .and. (pbod(i)%mass /= 0.0_rk) .and. ldflag(i)) then

        call drift_one(mstar, pbod(i)%r, pbod(i)%v, dtau, iflag)
        ldflag(i) = .false. ! Set so we don't return here

        if(iflag /= 0) then

          write(*,'(a)')              "SWIFT Error:: symba5_enc_drift"
          write(*,'(a,i6,a)')         " Particle:    ", j, " is lost!!!!!!!!!"
          write(*,'(a,2(1pe14.6))')   " Mstar, dtau: ", mstar, dtau
          write(*,'(a,3(1pe14.6))')   " Helio. pos.: ", pbod(j)%r
          write(*,'(a,3(1pe14.6),/)') " Bary. vel.:  ", pbod(j)%v
          call util_exit(FAILURE)

        end if

      end if

      ! Perform a drift for all bodies in the encounter region around body i
      do k = 1, penc(i)%nenc

        j = penc(i)%ienc(k)

        if((pbod(j)%rlevel == irec) .and. (pbod(j)%mass /= 0.0_rk) .and. ldflag(j)) then

          call drift_one(mstar, pbod(j)%r, pbod(j)%v, dtau, iflag)
          ldflag(j) = .false. ! Set so we don't return here

          if(iflag /= 0) then

            write(*,'(a)')              "SWIFT Error:: symba5_enc_drift"
            write(*,'(a,i6,a)')         " Particle:    ", j, " is lost!!!!!!!!!"
            write(*,'(a,2(1pe14.6))')   " Mstar, dtau: ", mstar, dtau
            write(*,'(a,3(1pe14.6))')   " Helio. pos.: ", pbod(j)%r
            write(*,'(a,3(1pe14.6),/)') " Bary. vel.:  ", pbod(j)%v
            call util_exit(FAILURE)

          end if

        end if

      end do

    end if

  end do

  return
  end subroutine symba5_enc_drift
!!
  subroutine symba5_helio_obl(param, nbod, pbod)
  !-------------------------------------------------------------------------------------
  !         SYMBA5_HELIO_OBL.F90
  !-------------------------------------------------------------------------------------
  ! Calculates the acceleration due to the oblateness of the central body on the massive
  ! bodies in the Heliocentric frame
  !
  ! Input:  param    ==> Global parameters (See swift module & io_init_param)
  !         nbod       ==> Number of bodies
  !         pbod(nbod) ==> Current database entries of bodies
  !
  ! Output: pbod(nbod) ==> Updated database entries of bodies
  !
  ! Remarks Based on helio_getacch.f
  ! Author: Hal Levison
  ! Date: 09/12/99
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !-------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  type(param_t) :: param

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i
  real(rk) :: aobl(nbod,NDIM), aobl0(NDIM)
  real(rk) :: ir3h(nbod), irh(nbod)

  !-----------------!
  ! Executable code !
  !-----------------!

  !***********************************************************************!
  !*** Consider combining: util_ir_ir3 and obl_acc with the loop below ***!
  !***********************************************************************!

  ! Now compute the inverse of the distance for each body, and the acceleration
  ! they experience due to the oblateness of the central body
  call util_ir_ir3(nbod, 2, pbod, irh, ir3h)
  call obl_acc(param, nbod, pbod, irh, aobl)

  ! Dummy vector to store the oblateness terms for the central body
  ! A private copy is made for each thread, so they do not have to wait for any of the other threads
  aobl0 = aobl(1,:)

  ! Add contribution of acceleration due to the oblateness of then central body to each body
  !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) FIRSTPRIVATE(aobl0) &
  !$OMP SHARED(nbod, pbod, aobl)
  do i = 2, nbod

    if(pbod(i)%mass /= 0.0_rk) pbod(i)%a = pbod(i)%a + aobl(i,:) - aobl0

  end do
  !$OMP END PARALLEL DO

  return
  end subroutine symba5_helio_obl
!!
  subroutine symba5_helio_getacch(param, nbod, nbodm, pbod)
  !-------------------------------------------------------------------------------------
  !         SYMBA5_HELIO_GETACCH.F90
  !-------------------------------------------------------------------------------------
  ! Calculates the acceleration on the massive bodies in the Heliocentric frame
  !
  ! Input:  param    ==> Global parameters (See swift module & io_init_param)
  !         nbod       ==> Number of bodies
  !         nbodm      ==> Index of last massive gravitating body
  !         pbod(nbod) ==> Current database entries of bodies
  !
  ! Output: pbod(nbod) ==> Updated database entries of bodies
  !
  ! Remarks Based on helio_getacch.f
  ! Author: Hal Levison
  ! Date: 09/12/99
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !-------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm
  type(param_t) :: param

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i, j
  real(rk) :: dr(NDIM), dr2, idr32, faci, facj
  real(rk) :: daxj, dayj, dazj
  !real(rk) :: daxerr, dayerr, dazerr
  type(body_t) :: pbodi

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Initialize the accelerations
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i) SHARED(nbod, pbod)
  do i = 1, nbod

    pbod(i)%a = 0.0_rk

  end do
  !$OMP END PARALLEL DO

  ! Compute the acceleration between gravitating bodies, the acceleration on non-gravitating
  ! bodies due to the gravitating bodies, and their back reaction on the gravitating bodies.
  ! No mutual gravitation interaction is computed for bodies nbodm + 1 to nbod.
  !
  ! N.B. The parallel approach outlined below is not reccomended for a large number of
  ! gravitating bodies (i.e. nbodm >~ 1000), otherwise the overhead for creating a destroying
  ! the parallel region in the inner loop will start to become a factor.
  do i = 2, nbodm

    ! Initialize the sum of the acceleration from bodies j
    daxj = 0.0_rk
    dayj = 0.0_rk
    dazj = 0.0_rk

    ! Initialize the correction factors for the Kahan floating-point summation formula
    !daxerr = 0.0_rk
    !dayerr = 0.0_rk
    !dazerr = 0.0_rk

    ! Extract a copy of the information for particle i, and make a copy available for each thread
    pbodi = pbod(i)

    !FIRSTPRIVATE(i, pbodi, daxerr, dayerr, dazerr)
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(i, pbodi) &
    !$OMP PRIVATE(j, dr, dr2, idr32, faci, facj) SHARED(nbod, pbod) REDUCTION(+ : daxj, dayj, dazj)
    do j = i + 1, nbod

      dr = pbod(j)%r - pbodi%r
      dr2 = sum(dr**2)
      idr32 = 1.0_rk/(dr2*sqrt(dr2))

      ! Acceleration contributed by gravitating body i on body j (gravitating/non-gravitating), and vice versa
      faci = pbodi%mass*idr32
      facj = pbod(j)%mass*idr32

      ! Sum the acceleration for gravitating body i due to body j (gravitating/non-gravitating)
      daxj = daxj + facj*dr(1)
      dayj = dayj + facj*dr(2)
      dazj = dazj + facj*dr(3)
      !daxj = util_kahan_sum(daxj, facj*dr(1), daxerr)
      !dayj = util_kahan_sum(dayj, facj*dr(2), dayerr)
      !dazj = util_kahan_sum(dazj, facj*dr(3), dazerr)

      ! Update the acceleration for body j (gravitating/non-gravitating) due to gravitating body i
      pbod(j)%a = pbod(j)%a - faci*dr

    end do
    !$OMP END PARALLEL DO

    ! Update the acceleration for gravitating body i due to bodies j (gravitating/non-gravitating)
    pbod(i)%a = pbod(i)%a + (/ daxj, dayj, dazj /)

  end do

  ! Compute the acceleration due to the oblateness of the central body, if desired
  if(param%loblate) call symba5_helio_obl(param, nbod, pbod)

  return
  end subroutine symba5_helio_getacch
!!
  subroutine symba5_step_helio(param, nbod, nbodm, pbod)
  !--------------------------------------------------------------------------------
  !         SYMBA5_STEP_HELIO.F90
  !--------------------------------------------------------------------------------
  ! Takes a step in heliocentric frame by doing a KICK, a DRIFT than a KICK
  !
  ! N.B. Does *not* handle close encounters/mergers
  !
  ! Input:  param       ==> Global parameters (See swift module & io_init_param)
  !         nbod          ==> Number of bodies
  !         nbodm         ==> Index of last massive gravitating body
  !         pbod(nbod)    ==> Current database entries of bodies
  !
  ! Output: pbod(nbod)    ==> Updated database entries of bodies
  !
  ! Authors: Hal Levison
  ! Date: 03/20/97
  ! Last revision: 09/06/08 CCC - Converted to Fortran 90/95 syntax
  !--------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm
  type(param_t) :: param

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  real(rk) :: dth

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Half the time step
  dth = 0.5_rk*param%dt

  ! Convert velocities from the heliocentric frame to the barycentric frame
  call coord_vh2vb(param, nbod, pbod)

  ! Do the linear drift due to momentum of the Sun
  call helio_lindrift(nbod, pbod, dth)

  ! If this our first time here, compute the accelerations in heliocentric frame
  if(lfirst_force) call symba5_helio_getacch(param, nbod, nbodm, pbod)

  ! Apply a heliocentric kick for a half step
  call helio_kickvb(nbod, pbod, dth)

  ! Drift in heliocentric frame for the full step
  call helio_drift(nbod, pbod, param%dt)

  ! Compute the accelerations in heliocentric frame
  call symba5_helio_getacch(param, nbod, nbodm, pbod)

  ! Apply a heliocentric kick for a half step
  call helio_kickvb(nbod, pbod, dth)

  ! Do the linear drift due to momentum of the Sun
  call helio_lindrift(nbod, pbod, dth)

  ! Convert velocities from the barycentric frame to the heliocentric frame
  call coord_vb2vh(param, nbod, pbod)

  ! Set symba5 force flag so we can recycle the old values of the accelerations
  lfirst_force = .false.

  return
  end subroutine symba5_step_helio
!!
  subroutine symba5_getacch(param, nbod, nbodm, pbod, penc)
  !-------------------------------------------------------------------------------------------
  !         SYMBA5_GETACCH.F90
  !-------------------------------------------------------------------------------------------
  ! Calculates the acceleration on the massive bodies in the Heliocentric frame, and removes
  ! the calculations for those undergoing an encounter.
  !
  ! Input:  param       ==> Global parameters (See swift module & io_init_param)
  !         nbod        ==> Number of bodies
  !         nbodm       ==> Location of last massive body
  !         pbod(nbod)  ==> Current database entries of bodies
  !         penc(nbodm) ==> List of encounters
  !
  ! Output: pbod(nbod)  ==> Updated database entries of bodies
  !
  ! Remarks: Based on helio_getacch.f, but does not include the forces of
  !          an body B on body A, if body B and A are having an encounter.
  !
  ! Author:  Hal Levison
  ! Date: 03/20/97
  ! Last revision: 11/22/97
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !-------------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm
  type(enc_t) :: penc(nbodm)
  type(param_t) :: param

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i, j, k
  real(rk) :: dr(NDIM), dr2, idr32, faci, facj
  real(rk) :: daxj, dayj, dazj
  !real(rk) :: daxerr, dayerr, dazerr
  type(body_t) :: pbodi

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Initialize the accelerations
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i) SHARED(nbod, pbod)
  do i = 1, nbod

    pbod(i)%a = 0.0_rk

  end do
  !$OMP END PARALLEL DO

  ! Compute the acceleration between gravitating bodies, the acceleration on non-gravitating
  ! bodies due to the gravitating bodies, and their back reaction on the gravitating bodies.
  ! No mutual gravitation interaction is computed for bodies nbodm + 1 to nbod.
  !
  ! N.B. The parallel approach outlined below is not reccomended for a large number of
  ! gravitating bodies (i.e. nbodm >~ 10), otherwise the overhead for creating a destroying
  ! the parallel region in the inner loop will start to become a factor.
  do i = 2, nbodm

    ! Initialize the sum of the acceleration from bodies j
    daxj = 0.0_rk
    dayj = 0.0_rk
    dazj = 0.0_rk

    ! Initialize the correction factors for the Kahan floating-point summation formula
    !daxerr = 0.0_rk
    !dayerr = 0.0_rk
    !dazerr = 0.0_rk

    ! Extract a copy of the information for particle i, and make a copy available for each thread
    pbodi = pbod(i)

    !FIRSTPRIVATE(i, pbodi, daxerr, dayerr, dazerr
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(i, pbodi) &
    !$OMP PRIVATE(j, dr, dr2, idr32, faci, facj) SHARED(nbod, pbod) REDUCTION(+ : daxj, dayj, dazj)
    do j = i + 1, nbod

      dr = pbod(j)%r - pbodi%r
      dr2 = sum(dr**2)
      idr32 = 1.0_rk/(dr2*sqrt(dr2))

      ! Acceleration contributed by gravitating body i on body j (gravitating/non-gravitating), and vice versa
      faci = pbodi%mass*idr32
      facj = pbod(j)%mass*idr32

      ! Sum the acceleration for gravitating body i due to body j (gravitating/non-gravitating)
      daxj = daxj + facj*dr(1)
      dayj = dayj + facj*dr(2)
      dazj = dazj + facj*dr(3)
      !daxj = util_kahan_sum(daxj, facj*dr(1), daxerr)
      !dayj = util_kahan_sum(dayj, facj*dr(2), dayerr)
      !dazj = util_kahan_sum(dazj, facj*dr(3), dazerr)

      ! Update the acceleration for body j (gravitating/non-gravitating) due to gravitating body i
      pbod(j)%a = pbod(j)%a - faci*dr

    end do
    !$OMP END PARALLEL DO

    ! Update the acceleration for gravitating body i due to bodies j (gravitating/non-gravitating)
    pbod(i)%a = pbod(i)%a + (/ daxj, dayj, dazj /)

    ! Now subtract off anyone in an encounter
    !
    ! See above comment on potential problems with this method of parallization

    ! Check if gravitating body i has any encounters
    if(penc(i)%nenc > 0) then

      ! Initialize the sum of the acceleration from bodies j
      daxj = 0.0_rk
      dayj = 0.0_rk
      dazj = 0.0_rk

      ! Initialize the correction factors for the Kahan floating-point summation formula
      !daxerr = 0.0_rk
      !dayerr = 0.0_rk
      !dazerr = 0.0_rk

      ! Extract a copy of the information for particle i, and make a copy available for each thread
      pbodi = pbod(i)

      !if(penc(i)%nenc > NTHRESHOLD) then

        !FIRSTPRIVATE(i, pbodi, daxerr, dayerr, dazerr)
        !OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) FIRSTPRIVATE(i, pbodi) &
        !OMP PRIVATE(j, k, dr, dr2, idr32, faci, facj) SHARED(nbod, pbod, penc) REDUCTION(+ : daxj, dayj, dazj)
        !do k = 1, penc(i)%nenc

          ! Extract the index of the encounter partner
          !j = penc(i)%ienc(k)

          !dr = pbod(j)%r - pbodi%r
          !dr2 = sum(dr**2)
          !idr32 = 1.0_rk/(dr2*sqrt(dr2))

          ! Acceleration contributed by gravitating body i on body j (gravitating/non-gravitating), and vice versa
          !faci = pbodi%mass*idr32
          !facj = pbod(j)%mass*idr32

          ! Sum the acceleration for gravitating body i due to body j (gravitating/non-gravitating)
          !daxj = daxj + facj*dr(1)
          !dayj = dayj + facj*dr(2)
          !dazj = dazj + facj*dr(3)
          !daxj = util_kahan_sum(daxj, facj*dr(1), daxerr)
          !dayj = util_kahan_sum(dayj, facj*dr(2), dayerr)
          !dazj = util_kahan_sum(dazj, facj*dr(3), dazerr)

          ! Remove the acceleration from body j (gravitating/non-gravitating) due to gravitating body i
          !pbod(j)%a = pbod(j)%a + faci*dr

        !end do
        !OMP END PARALLEL DO

      !else

        do k = 1, penc(i)%nenc

          ! Extract the index of the encounter partner
          j = penc(i)%ienc(k)

          dr = pbod(j)%r - pbodi%r
          dr2 = sum(dr**2)
          idr32 = 1.0_rk/(dr2*sqrt(dr2))

          ! Acceleration contributed by gravitating body i on body j (gravitating/non-gravitating), and vice versa
          faci = pbodi%mass*idr32
          facj = pbod(j)%mass*idr32

          ! Sum the acceleration for gravitating body i due to body j (gravitating/non-gravitating)
          daxj = daxj + facj*dr(1)
          dayj = dayj + facj*dr(2)
          dazj = dazj + facj*dr(3)
          !daxj = util_kahan_sum(daxj, facj*dr(1), daxerr)
          !dayj = util_kahan_sum(dayj, facj*dr(2), dayerr)
          !dazj = util_kahan_sum(dazj, facj*dr(3), dazerr)

          ! Remove the acceleration from body j (gravitating/non-gravitating) due to gravitating body i
          pbod(j)%a = pbod(j)%a + faci*dr

        end do

      !end if

      ! Remove the acceleration from body i due to bodies j (gravitating/non-gravitating)
      pbod(i)%a = pbod(i)%a - (/ daxj, dayj, dazj /)

    end if

  end do

  ! Compute the acceleration due to the oblateness of the central body, if desired
  if(param%loblate) call symba5_helio_obl(param, nbod, pbod)

  return
  end subroutine symba5_getacch
!!
  subroutine symba5_merge(param, time, dtau, ireci, nbod, nbodm, pbod, penc, pmerge, eoffset)
  !---------------------------------------------------------------------------------
  !         SYMBA5_MERGE.F90
  !---------------------------------------------------------------------------------
  ! Checks if there are any mergers, and processes mergers
  !
  ! Input:  param       ==> Global parameters (See swift module & io_init_param)
  !         time        ==> Current time
  !         dtau        ==> Time step
  !         ireci       ==> Current recursion level
  !         nbod        ==> Number of bodies
  !         nbodm       ==> Index of last massive body
  !         pbod(nbod)  ==> Current database entries of bodies
  !         penc(nbodm) ==> Encounter list for massive gravitating bodies
  !         pmerge      ==> Current list of mergers
  !         eoffset     ==> Current energy offset
  !
  ! Updated only if a merger occurs
  ! Output: pbod(nbod)  ==> Updated database entries of bodies
  !         pmerge      ==> Updated list of mergers
  !         eoffset     ==> Updated energy offset
  !
  ! Authors:  Hal Levison
  ! Date: 01/02/97
  ! Last revision: 01/24/97
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !---------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm, ireci
  real(rk) :: time, dtau
  type(param_t) :: param

  ! Input and Output variables
  integer(ik), pointer :: pmerge(:)
  real(rk) :: eoffset
  type(body_t) :: pbod(nbod)
  type(enc_t) :: penc(nbodm)

  ! Internal variables
  integer(ik) :: i, j, k, ialpha
  real(rk) :: dr(NDIM), dv(NDIM), dr2, dv2, vdotr, tcross2
  real(rk) :: rsum, rsum2, msum, a, e, q, dt2
  real(rk) :: ke, pot, energy1, energy2, l(NDIM)

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Loop over each massive body, in reverse order to ensure we account for all mergers
  do i = nbodm, 2, -1

    ! Check if body i has any encounters
    if(penc(i)%nenc > 0) then

      ! Loop over encounters for body i
      do k = 1, penc(i)%nenc

        ! Extract index of encounter partner
        j = penc(i)%ienc(k)

        ! Proceed if both bodies are currently at a recursion level
        ! that is greater or equal to the input recursion level
        if((pbod(i)%rlevel >= ireci) .and. (pbod(j)%rlevel >= ireci)) then

          ! Displacement and distance squared
          dr = pbod(j)%r - pbod(i)%r
          dr2 = sum(dr**2)

          ! Sum of the physical radii, and physical radii squared
          rsum = pbod(i)%rphy + pbod(j)%rphy
          rsum2 = rsum**2

          ! Abort if sum of radii is zero
          if(rsum == 0.0_rk) then

            write(*,'(a)') "SWIFT Error:: symba5_merge"
            write(*,'(a)') " Sum of radii for merging particles is zero!"
            call util_exit(FAILURE)

          end if

          ! If the sum of the physical radii squared is greater
          ! than or equal to their distance squared, then merge
          if(rsum2 >= dr2) then

            ! If we are following the energy of the system, compute the total energy *before* merging the bodies
            if(param%lenergy) call anal_energy(param, nbod, nbodm, pbod, ke, pot, energy1, l)

            ! Append body j to the merger list
            pmerge => util_append_index(pmerge, j)

            ! Merge body i and body j, removing body j from the encounter list
            call discard_mass_merge5(nbodm, pbod(1)%mass, time, pbod(i), pbod(j), penc)

            ! If we are following the energy of the system, compute the total energy *after* merging the bodies
            if(param%lenergy) then

              ! Compute the total energy of the system *after* merging the bodies
              call anal_energy(param, nbod, nbodm, pbod, ke, pot, energy2, l)

              ! Update the energy offset
              eoffset = eoffset + energy1 - energy2

            end if

          else

            dv = pbod(j)%v - pbod(i)%v   ! Relative velocity components
            dv2 = sum(dv**2)             ! Magnitude of relative velocity squared
            vdotr = dot_product(dr, dv)  ! Dot product of displacement and relative velocity components

            ! If the value of vdotr has changed from negative (approaching) to positive (receding), then proceed
            if(penc(i)%lvdotr(k) .and. (vdotr > 0.0_rk)) then

              ! Crossing time squared for body i and body j
              tcross2 = dr2/dv2

              ! If the square of the crossing time is less than or equal to the local time step, then proceed
              if(tcross2 <= dtau**2) then

                msum = pbod(i)%mass + pbod(j)%mass               ! Sum of masses
                call orbel_xv2aeq(dr, dv, msum, ialpha, a, e, q) ! Compute closest approach to centre of mass

                ! If the closest approach to centre of mass is less than the sum of physical radii, then merge
                if(q < rsum) then

                  ! If we are following the energy of the system, compute the total energy *before* merging the bodies
                  if(param%lenergy) call anal_energy(param, nbod, nbodm, pbod, ke, pot, energy1, l)

                  ! Append body j to the merger list
                  pmerge => util_append_index(pmerge, j)

                  ! Merge body i and body j, removing body j from the encounter list
                  call discard_mass_merge5(nbodm, pbod(1)%mass, time, pbod(i), pbod(j), penc)

                  ! If we are following the energy of the system, compute the total energy *after* merging the bodies
                  if(param%lenergy) then

                    ! Compute the total energy of the system *after* merging the bodies
                    call anal_energy(param, nbod, nbodm, pbod, ke, pot, energy2, l)

                    ! Update the energy offset
                    eoffset = eoffset + energy1 - energy2

                  end if

                end if

              end if

            end if

          end if

        end if

      end do

    end if

  end do

  return
  end subroutine symba5_merge
!!
  subroutine symba5_kick(nbod, nbodm, irec, pbod, penc, dtau, force_dir)
  !----------------------------------------------------------------------------------------------
  !         SYMBA5_KICK.F90
  !----------------------------------------------------------------------------------------------
  ! Do a symba5 kick
  !
  ! Input:  nbod        ==> Number of bodies
  !         nbodm       ==> Index of the last massive body
  !         irec        ==> Current recursion level
  !         pbod(nbod)  ==> Mass of bodies
  !         penc(nbodm) ==> Encounter list for massive gravitating bodies
  !         dtau        ==> Time step
  !         force_dir   ==> Add or subtract the force
  !                         - Add: +1.0_rk, Subtract: -1.0_rk
  !
  ! Output: pbod(nbod)  ==> Updated barycentric velocity of bodies
  !
  ! Remarks: Uses Man Hoi's force
  ! Authors: Hal Levison
  ! Date: 03/20/97
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !----------------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm, irec
  real(rk) :: dtau, force_dir
  type(enc_t) :: penc(nbodm)

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i, j, k, irec_above, irecl
  real(rk) :: dr(NDIM), dr2, idr32, fac, faci, facj
  real(rk) :: r2crit, r2crit_below, rcrit, rratio
  real(rk) :: daxj, dayj, dazj
  !real(rk) :: daxerr, dayerr, dazerr
  type(body_t) :: pbodi

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Decide the local recursion level based on force_dir
  if(force_dir < 0.0_rk) then

    irecl = irec - 1

  else

    irecl = irec

  end if

  ! Store the recursion level above irec
  irec_above = irec - 1

  ! Calculate the accelerations
  do i = 2, nbodm

    if(penc(i)%nenc > 0) then

      ! Initialize the sum of the acceleration from bodies j
      daxj = 0.0_rk
      dayj = 0.0_rk
      dazj = 0.0_rk

      ! Initialize the correction factors for the Kahan floating-point summation formula
      !daxerr = 0.0_rk
      !dayerr = 0.0_rk
      !dazerr = 0.0_rk

      ! Extract a copy of the information for particle i, and make a copy available for each thread
      pbodi = pbod(i)

      !if(penc(i)%nenc > NTHRESHOLD) then

        !FIRSTPRIVATE(i, irec_above, irecl, pbodi, daxerr, dayerr, dazerr)
        !OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
        !OMP FIRSTPRIVATE(i, irec_above, irecl, pbodi, dtau, force_dir) &
        !OMP PRIVATE(j, k, dr, dr2, idr32, fac, faci, facj, rratio, rcrit, r2crit, r2crit_below) &
        !OMP SHARED(nbod, pbod, penc) REDUCTION(+ : daxj, dayj, dazj)
        !do k = 1, penc(i)%nenc

          !j = penc(i)%ienc(k) ! Extract the index of the collision partner

          !if((pbodi%rlevel >= irec_above) .and. (pbod(j)%rlevel >= irec_above)) then

            !r2crit = ((pbodi%rhill + pbod(j)%rhill)*RHSCALE*RSHELL**irecl)**2
            !r2crit_below = r2crit*RSHELL**2 ! Square of critical radius for recursion level irecl + 1

            !dr = pbod(j)%r - pbodi%r
            !dr2 = sum(dr**2)
            !idr32 = 1.0_rk/(dr2*sqrt(dr2))

            !if(dr2 < r2crit_below) then

              !fac = 0.0_rk

            !else if(dr2 < r2crit) then

              !rcrit = sqrt(r2crit)
              !rratio = (rcrit - sqrt(dr2))/(rcrit*(1.0_rk - RSHELL))
              !fac = (1.0_rk - 3.0_rk*rratio**2 + 2.0_rk*rratio**3)*idr32

            !else

              !fac = idr32

            !end if

            ! Acceleration contributed by gravitating body i on body j (gravitating/non-gravitating), and vice versa
            !faci = pbodi%mass*fac
            !facj = pbod(j)%mass*fac

            ! Sum the acceleration for gravitating body i due to body j (gravitating/non-gravitating)
            !daxj = daxj + facj*dr(1)
            !dayj = dayj + facj*dr(2)
            !dazj = dazj + facj*dr(3)
            !daxj = util_kahan_sum(daxj, facj*dr(1), daxerr)
            !dayj = util_kahan_sum(dayj, facj*dr(2), dayerr)
            !dazj = util_kahan_sum(dazj, facj*dr(3), dazerr)

            ! Apply the acceleration for body j (gravitating/non-gravitating) due to gravitating body i
            !pbod(j)%v = pbod(j)%v - force_dir*faci*dr*dtau

          !end if

        !end do
        !OMP END PARALLEL DO

      !else

        do k = 1, penc(i)%nenc

          j = penc(i)%ienc(k) ! Extract the index of the collision partner

          if((pbodi%rlevel >= irec_above) .and. (pbod(j)%rlevel >= irec_above)) then

            r2crit = ((pbodi%rhill + pbod(j)%rhill)*RHSCALE*RSHELL**irecl)**2
            r2crit_below = r2crit*RSHELL**2 ! Square of critical radius for recursion level irecl + 1

            dr = pbod(j)%r - pbodi%r
            dr2 = sum(dr**2)
            idr32 = 1.0_rk/(dr2*sqrt(dr2))

            if(dr2 < r2crit_below) then

              fac = 0.0_rk

            else if(dr2 < r2crit) then

              rcrit = sqrt(r2crit)
              rratio = (rcrit - sqrt(dr2))/(rcrit*(1.0_rk - RSHELL))
              fac = (1.0_rk - 3.0_rk*rratio**2 + 2.0_rk*rratio**3)*idr32

            else

              fac = idr32

            end if

            ! Acceleration contributed by gravitating body i on body j (gravitating/non-gravitating), and vice versa
            faci = pbodi%mass*fac
            facj = pbod(j)%mass*fac

            ! Sum the acceleration for gravitating body i due to body j (gravitating/non-gravitating)
            daxj = daxj + facj*dr(1)
            dayj = dayj + facj*dr(2)
            dazj = dazj + facj*dr(3)
            !daxj = util_kahan_sum(daxj, facj*dr(1), daxerr)
            !dayj = util_kahan_sum(dayj, facj*dr(2), dayerr)
            !dazj = util_kahan_sum(dazj, facj*dr(3), dazerr)

            ! Apply the acceleration for body j (gravitating/non-gravitating) due to gravitating body i
            pbod(j)%v = pbod(j)%v - force_dir*faci*dr*dtau

          end if

        end do

      !end if

      ! Apply the net acceleration for gravitating body i due to bodies j (gravitating/non-gravitating)
      pbod(j)%v = pbod(j)%v + force_dir*[ daxj, dayj, dazj ]*dtau

    end if

  end do

  return
  end subroutine symba5_kick
!!
  subroutine symba5_rlevel_check(nbod, nbodm, ireci, dtau, pbod, penc, lrecur)
  !----------------------------------------------------------------------------------------------
  !         SYMBA5_RLEVEL_CHECK.F90
  !----------------------------------------------------------------------------------------------
  ! Do a symba5 kick
  !
  ! Input:  nbod        ==> Number of bodies
  !         nbodm       ==> Index of the last massive body
  !         ireci       ==> Current recursion level
  !         dtau        ==> Time step
  !         pbod(nbod)  ==> Entries for
  !         penc(nbodm) ==> Encounter list for massive gravitating bodies
  !
  ! Output: pbod(nbod)  ==> Entries for
  !         penc(nbodm) ==> Encounter list for massive gravitating bodies
  !         lrecur      ==> Recursion flag
  !
  ! Authors: Chris Capobianco
  ! Date: 02/16/10
  !----------------------------------------------------------------------------------------------
  !$ use omp_lib
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm, ireci
  real(rk) :: dtau
  type(enc_t) :: penc(nbodm)

  ! Input and Output variable
  type(body_t) :: pbod(nbod)
  logical(lk) :: lrecur

  ! Internal variables
  logical(lk) :: lenc, lvdotr, lrecur_local(nthreads)
  integer(ik) :: i, j, k, id, irecp
  type(body_t) :: pbodi

  !-----------------!
  ! Executable code !
  !-----------------!

  ! The recursion level below ireci
  irecp = ireci + 1
  
  ! Initialize global and local recursion flags
  lrecur = .false.
  lrecur_local = .false.

  ! Calculate the accelerations
  do i = 2, nbodm

    if(penc(i)%nenc > 0) then

      ! Extract a copy of the information for particle i, and make a copy available for each thread
      pbodi = pbod(i)

      ! Initialize the local recursion flags
      lrecur_local = .false.

      !if(penc(i)%nenc > NTHRESHOLD) then

        !OMP PARALLEL DEFAULT(NONE) FIRSTPRIVATE(i, ireci, irecp, dtau, pbodi) &
        !OMP PRIVATE(j, k, id, lenc) SHARED(nbod, pbod, penc, lrecur_local)

        !id = 1                           ! Set thread identifier for *serial* case
        ! id = omp_get_thread_num() + 1 ! Set thread identifier for *parallel* case

        !OMP DO SCHEDULE(STATIC)
        !do k = 1, penc(i)%nenc

          ! Extract the index of the collision partner
          !j = penc(i)%ienc(k)

          !if((pbodi%rlevel >= ireci) .and. (pbod(j)%rlevel >= ireci)) then

            ! Determine if bodies i and j are still in an encounter region
            !call symba5_check(irecp, dtau, pbodi, pbod(j), penc(i)%lvdotr(k), lenc)

            ! If body j are in an encounter region with body i, then update the
            ! recursion level information for body j and set the local recursion flag
            !if(lenc) then

              !pbod(j)%rlevel = irecp
              !pbod(j)%rlevelmax = max(irecp, pbod(j)%rlevelmax)
              !lrecur_local(id) = .true.

            !end if

          !end if

        !end do
        !OMP END DO NOWAIT

        !OMP END PARALLEL

      !else

        do k = 1, penc(i)%nenc

          ! Extract the index of the collision partner
          j = penc(i)%ienc(k)

          if((pbodi%rlevel >= ireci) .and. (pbod(j)%rlevel >= ireci)) then

            ! Determine if bodies i and j are still in an encounter region
            call symba5_check(irecp, dtau, pbodi, pbod(j), penc(i)%lvdotr(k), lenc)

            ! If body j are in an encounter region with body i, then update the
            ! recursion level information for body j and set the local recursion flag
            if(lenc) then

              pbod(j)%rlevel = irecp
              pbod(j)%rlevelmax = max(irecp, pbod(j)%rlevelmax)
              lrecur_local(1) = .true.

            end if

          end if

        end do

      !end if

    end if

    ! If any of the local recursion flags are set, then update the
    ! recursion level information for body i and set the global recursion flag
    if(any(lrecur_local)) then

      pbod(i)%rlevel = irecp
      pbod(i)%rlevelmax = max(irecp, pbod(i)%rlevelmax)
      lrecur = .true.

    end if

  end do

  return
  end subroutine symba5_rlevel_check
!!
  recursive subroutine symba5_step_recur(param, t, dtau, ireci, nbod, nbodm, pbod, penc, pmerge, eoffset)
  !---------------------------------------------------------------------------------------
  !         SYMBA5_STEP_RECUR.F90
  !---------------------------------------------------------------------------------------
  ! Follow the trajectory for all bodies *undergoing* an encounter, using a
  ! recursively smaller time step to resolve the trajectory as warranted
  !
  ! Input:  param       ==> Global parameters (See swift module & io_init_param)
  !         t           ==> Current time
  !         dtau        ==> Time step
  !         ireci       ==> Input recursion level
  !         nbod        ==> Number of bodies
  !         nbodm       ==> Index of last massive body
  !         pbod(nbod)  ==> Mass of bodies
  !         penc(nbodm) ==> Encounter list for massive gravitating bodies
  !         pmerge      ==> Current list of mergers
  !         eoffset     ==> Current energy offset
  !
  ! Output: pbod(nbod)  ==> Updated mass of bodies (if a merger occurred)
  !         pmerge      ==> List of mergers (if a merger occurred)
  !         eoffset     ==> Updated energy offset (if a merger occurred)
  !
  ! Remarks: If a merger occurs, does not change nbod and puts the mass
  !          of one of the bodies to zero.
  ! Authors: Hal Levison
  ! Date: 03/20/97
  ! Last revision: 05/13/99
  !---------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm, ireci
  real(rk) :: t, dtau
  type(enc_t) :: penc(nbodm)
  type(param_t) :: param

  ! Input and Output variables
  integer(ik), pointer :: pmerge(:)
  real(rk) :: eoffset
  type(body_t) :: pbod(nbod)

  ! Internal variables
  logical(lk) :: lrecur
  integer(ik) :: i, j, irecp
  real(rk) :: dtl, dth

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Compute the local time step, and half local time step
  dtl = dtau/real(NTENC**ireci, rk)
  dth = 0.5_rk*dtl

  if(dtl/dtau <= TINY_NUMBER) then

    write(*,'(a)') "SWIFT Error:: symba_step_recur"
    write(*,'(a)') " Local timestep too small, roundoff will be important!!!!"
    call util_exit(FAILURE)

  end if

  ! Recursion level above input recursion level
  irecp = ireci + 1

  if(ireci == 0) then

    ! For all pairs of bodies undergoing an encounter, test if any are
    ! still in an encounter region for the next recursion level
    call symba5_rlevel_check(nbod, nbodm, ireci, dtl, pbod, penc, lrecur)

    !force_dir = 1.0_rk
    call symba5_kick(nbod, nbodm, irecp, pbod, penc, dth, 1.0_rk)

    call symba5_helio_drift(ireci, dtl, nbod, pbod)

    if(lrecur) call symba5_step_recur(param, t, dtau, irecp, nbod, nbodm, pbod, penc, pmerge, eoffset)

    !force_dir = 1.0_rk
    call symba5_kick(nbod, nbodm, irecp, pbod, penc, dth, 1.0_rk)

    ! Look for mergers from encounter list
    if(param%laccrete) call symba5_merge(param, t, dtl, ireci, nbod, nbodm, pbod, penc, pmerge, eoffset)

    ! Decrease the recursion level for all bodies at irecp to ireci
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i), FIRSTPRIVATE(irecp, ireci) SHARED(nbod, pbod)
    do i = 2, nbod

      if(pbod(i)%rlevel == irecp) pbod(i)%rlevel = ireci

    end do
    !$OMP END PARALLEL DO

  else

    do j = 1, NTENC

      ! For all pairs of bodies undergoing an encounter, test if any are
      ! still in an encounter region for the next recursion level
      call symba5_rlevel_check(nbod, nbodm, ireci, dtl, pbod, penc, lrecur)

      !force_dir = 1.0_rk
      call symba5_kick(nbod, nbodm, irecp, pbod, penc, dth, 1.0_rk)

      !force_dir = -1.0_rk
      call symba5_kick(nbod, nbodm, irecp, pbod, penc, dth, -1.0_rk)

      call symba5_helio_drift(ireci, dtl, nbod, pbod)

      if(lrecur) call symba5_step_recur(param, t, dtau, irecp, nbod, nbodm, pbod, penc, pmerge, eoffset)

      !force_dir = 1.0_rk
      call symba5_kick(nbod, nbodm, irecp, pbod, penc, dth, 1.0_rk)

      !force_dir = -1.0_rk
      call symba5_kick(nbod, nbodm, irecp, pbod, penc, dth, -1.0_rk)

      ! Look for mergers from encounter list
      if(param%laccrete) call symba5_merge(param, t, dtl, ireci, nbod, nbodm, pbod, penc, pmerge, eoffset)

      ! Decrease the recursion level for all bodies at irecp to ireci
      !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i), FIRSTPRIVATE(irecp, ireci) SHARED(nbod, pbod)
      do i = 2, nbod

        if(pbod(i)%rlevel == irecp) pbod(i)%rlevel = ireci

      end do
      !$OMP END PARALLEL DO

    end do

  end if

  return
  end subroutine symba5_step_recur
!!
  subroutine symba5_step_interp(param, time, nbod, nbodm, pbod, penc, pmerge, eoffset)
  !---------------------------------------------------------------------------------------
  !         SYMBA5_STEP_INTERP.F90
  !---------------------------------------------------------------------------------------
  ! Takes a step in heliocentric frame by doing a KICK, a DRIFT than a KICK, but also
  ! handles close encounters and mergers by interpolating between successive smaller
  ! spheres around each massive body.
  !
  ! Input:  param       ==> Global parameters (See swift module & io_init_param)
  !         nbod        ==> Number of bodies
  !         nbodm       ==> Index of last massive body
  !         time        ==> Current time
  !         pbod(nbod)  ==> Mass of bodies
  !         penc(nbodm) ==> Encounter list for massive gravitating bodies
  !         eoffset     ==> Current energy offset
  !
  ! Output: pbod(nbod)  ==> Updated mass of bodies (if a merger occurred)
  !         pmerge      ==> List of mergers (if a merger occurred)
  !         eoffset     ==> Updated energy offset (if a merger occurred)
  !
  ! Authors: Hal Levison
  ! Date: 11/21/96
  ! Last revision: 5/13/99
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !---------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm
  real(rk) :: time
  type(param_t) :: param
  type(enc_t) :: penc(nbodm)

  ! Input and Output variables
  real(rk) :: eoffset
  type(body_t) :: pbod(nbod)

  ! Output variable
  integer(ik), pointer :: pmerge(:)

  ! Internal variables
  integer(ik) :: i, irec
  real(rk) :: dth

  !-----------------!
  ! Executable Code !
  !-----------------!

  ! Half the time step
  dth = 0.5_rk*param%dt

  ! Initialize recursion index
  irec = 0

  ! Convert velocities from heliocentric to barycentric frame
  call coord_vh2vb(param, nbod, pbod)

  ! Do the linear drift due to momentum of the Sun
  call helio_lindrift(nbod, pbod, dth)

  ! Get the accelerations in heliocentric frame
  ! For each object only include those bodies with which they are *not* undergoing an encounter
  call symba5_getacch(param, nbod, nbodm, pbod, penc)

  ! Apply a heliocentric kick for a half dt
  call helio_kickvb(nbod, pbod, dth)

  ! Do a drift for full dt for all bodies *not* undergoing an encounter
  call symba5_helio_drift(-1, param%dt, nbod, pbod)

  ! Follow the trajectory for all bodies *undergoing* an encounter, using a
  ! recursively smaller time step to resolve the trajectory as warranted
  call symba5_step_recur(param, time, param%dt, irec, nbod, nbodm, pbod, penc, pmerge, eoffset)

  ! Get the accelerations in helio frame
  ! For each object only include those bodies with which they are *not* undergoing an encounter
  call symba5_getacch(param, nbod, nbodm, pbod, penc)

  ! Apply a heliocentric kick for a half step
  call helio_kickvb(nbod, pbod, dth)

  ! Do the linear drift due to momentum of the Sun
  call helio_lindrift(nbod, pbod, dth)

  ! Convert velocities from barycentric to heliocentric frame
  call coord_vb2vh(param, nbod, pbod)

  ! Reset symba5 force calculation flag
  lfirst_force = .true.

  return
  end subroutine symba5_step_interp
!!
  subroutine symba5_step(param, time, nbod, nbodm, pbod, pmerge, eoffset)
  !---------------------------------------------------------------------------------------------------
  !         SYMBA5_STEP.F90
  !---------------------------------------------------------------------------------------------------
  ! Takes a step in heliocentric frame by doing a KICK, a DRIFT than a KICK, but also handles
  ! close encounters and mergers.
  !
  ! Input:  param      ==> Global parameters (See swift module & io_init_param)
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

  ! Do a step
  if(lenc) then

    ! If there are encounters, perform a recursive step
    call symba5_step_interp(param, time, nbod, nbodm, pbod, penc, pmerge, eoffset)

  else

    ! If there are no encounters, perform a simple step
    call symba5_step_helio(param, nbod, nbodm, pbod)

  end if

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
  end subroutine symba5_step
!!
!!
end module symba5
