module discard
! Module for routines in discard directory
use swift
use util
use coord
use anal
use io
implicit none

contains
!!
!!
  subroutine discard_mass_reorder5(ip, nbod, pbod)
  !------------------------------------------------------------------------------------
  !         DISCARD_MASS_REORDER5.F90
  !------------------------------------------------------------------------------------
  ! Remove body ip by shifting the values (ip + 1:nbod) up one index, then decreasing
  ! nbod by one
  !
  ! Input:  ip         ==> Index of body to be removed
  !         nbod       ==> Current number of bodies
  !         pbod(nbod) ==> Current database entries for bodies (See swift module)
  !
  ! Output: nbod       ==> Updated number of bodies
  !         pbod(nbod) ==> Updated database entries for bodies (See swift module)
  !
  ! Authors: Hal Levison
  ! Date: 01/02/97
  ! Last revision: 05/13/99
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !------------------------------------------------------------------------------------
  implicit none

  ! Input variable
  integer(ik) :: ip

  ! Input and Output variables
  integer(ik) :: nbod
  type(body_t) :: pbod(nbod)

  ! Internal variable
  integer(ik) :: i

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Shift all the entries up one index
  if(ip < nbod) then

    do i = ip, nbod - 1

      pbod(i) = pbod(i + 1)

    end do

  end if

  ! Decrease the number of bodies by one
  nbod = nbod - 1

  return
  end subroutine discard_mass_reorder5
!!
  subroutine discard_mass_merge5(nbodm, mstar, time, pbodi, pbodj, penc)
  !--------------------------------------------------------------------------------
  !         DISCARD_MASS_MERGE5.F90
  !--------------------------------------------------------------------------------
  ! Merges two bodies, stores the merger product in first index and writes
  ! information of the merger to discard_mass.dat
  !
  ! Input:  nbodm ==> Number of massive bodies
  !         mstar ==> Mass of central star
  !         time  ==> Current time
  !         pbodi ==> Current database entries for both bodies (See swift module)
  !         pbodj ==> Current database entries for both bodies (See swift module)
  !         penc  ==> Current database entries for both bodies (See swift module)
  !
  ! Output: pbodi ==> Updated database entries for both bodies (See swift module)
  !         pbodj ==> Current database entries for both bodies (See swift module)
  !         penc  ==> Current database entries for both bodies (See swift module)
  !
  ! Remarks:
  ! Authors: Hal Levison
  ! Date: 12/30/96
  ! Last revision: 01/30/97
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Now handles *two* bodies only
  !                             - Energy offset updated in symba5_merge
  !--------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbodm
  real(rk) :: mstar, time

  ! Input and Output variables
  type(body_t) :: pbodi, pbodj
  type(enc_t) :: penc(nbodm)

  ! Internal variables
  integer(ik) :: i, j, k
  type(body_t) :: pbodi_org, pbodj_org

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Announce merger event
  write(*,'(2(a,i6),a,1pe12.5)') " Merging bodies: ", pbodi.id, " and ", pbodj.id, " at t = ", time

  ! Make a copy of the relevant quantities for both bodies
  pbodi_org = pbodi
  pbodj_org = pbodj

  ! Note: I am just putting these particles together here, which is clearly wrong.
  !       I should integrate back to the time of close approach, then merge them.
  pbodi.mass = pbodi_org.mass + pbodj_org.mass
  pbodi.rphy = (pbodi_org.rphy**3 + pbodj_org.rphy**3)**(1.0_rk/3.0_rk)
  pbodi.rhill = util_hills_one(mstar, pbodi)
  pbodi.v = (pbodi_org.mass*pbodi_org.v + pbodj_org.mass*pbodj_org.v)/pbodi.mass

  ! Put in zeros for body pbodj
  pbodj.mass = 0.0_rk
  pbodj.rphy = 0.0_rk
  pbodj.rhill = 0.0_rk
  pbodj.rdrag = 0.0_rk
  pbodj.rho = 0.0_rk
  pbodj.v = 0.0_rk
  pbodj.r = 1.0e10_rk*pbodj_org.r ! So Danby does not fail

  ! Write merger information to discard_mass.dat
  call io_discard_merge(time, pbodi_org, pbodj_org, pbodi)

  ! If body pbodj is a massive body, decrease nbodm by one
  !if(pbodj.id <= nbodm) nbodm = nbodm - 1

  ! Remove body pbodj from the encounter list for each massive body
  do i = 2, nbodm

    ! Skip if there are no entries in the encounter list for massive body pbodi
    if(penc(i).nenc > 0) then

      ! Loop over all the entries in the encounter list for massive body pbodi
      inner: do j = 1, penc(i).nenc

        ! Search for entries that match the identifier of body pbodj
        if(penc(i).ienc(j) == pbodj.id) then

          if(j < penc(i).nenc) then

            do k = j, penc(i).nenc - 1

              penc(i).ienc(k) = penc(i).ienc(k + 1)

            end do

          end if

          ! Decrease the number of encounters by one
          penc(i).nenc = penc(i).nenc - 1

          ! Exit the inner do loop now that we have found the entry for body pbodj
          ! in the encounter list for massive body pbodi
          exit inner
	  
        end if
	
      end do inner

    end if
  
  end do

  return
  end subroutine discard_mass_merge5
!!
  subroutine discard_mass_peri(param, time, nbod, pbod, iwhy)
  !---------------------------------------------------------------------------------------
  !         DISCARD_MASS_PERI.F90
  !---------------------------------------------------------------------------------------
  ! Check if a body should be discarded because its perihelion distance gets too small
  !
  ! Input:  param      ==> Global parameters (See swift module & io_init_param)
  !         time       ==> Current time
  !         nbod       ==> Number of bodies
  !         pbod(nbod) ==> Database entries for bodies (See swift module)
  !         iwhy(nbod) ==> Status for bodies
  !
  ! Output: iwhy(nbod) ==> Status for bodies
  !         pbod.iperi ==> Perhelion status flag
  !                        - iperi = -1 if body is before perhelion
  !                        - iperi =  0 if body went through perihelion
  !                        - iperi = +1 if body is after perihelion
  !
  ! Remarks: Based on discard_peri
  ! Authors: Hal Levison
  ! Date: 12/30/96
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !---------------------------------------------------------------------------------------
  implicit none

  ! Inputs variables
  integer(ik) :: nbod
  real(rk) :: time
  type(body_t) :: pbod(nbod)
  type(param_t) :: param

  ! Input and Output variables
  integer(ik) :: iwhy(nbod)

  ! Internal variables
  integer(ik) :: i
  logical(lk), save :: lfirst = .true.
  logical(lk) :: lperi(nbod)
  real(rk) :: peri(nbod)

  !-----------------!
  ! Executable code !
  !-----------------!

  ! If first time through, set things up
  if(lfirst) then

    ! Determine whether perihelion of a body has taken place
    call util_mass_peri(nbod, pbod, peri, lperi)

    lfirst = .false. ! So we don't enter here again

  else

    ! Determine whether perihelion of a body has taken place
    call util_mass_peri(nbod, pbod, peri, lperi)

    ! If a body went through perihelion, did not have an encounter and
    ! its perihelion distance is less than qmin then flag for removal
    !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i) FIRSTPRIVATE(param, time) &
    !$OMP SHARED(nbod, pbod, lperi, peri, iwhy)
    do i = 2, nbod

      if(lperi(i) .and. (.not. pbod(i).lenc) .and. (peri(i) <= param.qmin)) then

        write(*,'(a,i7,a,1pe14.6)') ' Particle ', i, ' perihelion distance too small at t = ', time
        iwhy(i) = -4

      end if

    end do
    !$OMP END PARALLEL DO

  end if

  return
  end subroutine discard_mass_peri
!!
  subroutine discard_massive5(param, time, nbod, nbodm, pbod, pmerge, eoffset)
  !---------------------------------------------------------------------------------------------------
  !         DISCARD_MASSIVE5.F90
  !---------------------------------------------------------------------------------------------------
  ! Check if any bodies should be discarded, either due to merger or outside simulation space
  !
  ! Input:  param      ==> Global parameters (See swift module & io_init_param)
  !         time       ==> Current time
  !         nbod       ==> Current number of bodies
  !         nbodm      ==> Current index of last massive body
  !         pbod(nbod) ==> Current database entries for bodies
  !         pmerge     ==> Current list of mergers (See swift module & symba5_merge)
  !         eoffset    ==> Current energy offset
  !
  ! Output: nbod       ==> Updated number of bodies
  !         nbodm      ==> Updated index of last massive body
  !         pbod(nbod) ==> Updated database entries for bodies
  !         eoffset    ==> Updated energy offset
  !
  ! Authors: Hal Levison
  ! Date: 12/30/96
  ! Last revision: 05/13/99
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !---------------------------------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: time
  type(param_t) :: param

  ! Input and Output variables
  integer(ik) :: nbod, nbodm
  integer(ik), pointer :: pmerge(:)
  real(rk) :: eoffset
  type(body_t) :: pbod(nbod)

  ! Internal variables
  integer(ik) :: i, j, k, nmerge, iwhy(nbod)
  real(rk) :: rmin2, rmax2, rmaxu2, energy
  real(rk) :: ei, ef, ke, pot, l(NDIM), vdotr
  real(rk) :: rh2, rb2, vb2, msys
  type(body_t) :: p(nbod)

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Update perihelion status flag for all massive bodies
  do i = 2, nbodm

    ! Dot product between heliocentric position and velocity
    vdotr = dot_product(pbod(i).r, pbod(i).v)

    if(vdotr > 0.0_rk) then

      pbod(i).iperi = 1

    else

      pbod(i).iperi = -1

    end if

  end do

  ! Update the entries for the particles that have already been merged
  if(associated(pmerge)) then

    ! Sort the merger list in ascending order
    pmerge => util_sort_pointer(pmerge)

    do i = size(pmerge), 1, -1

      j = pmerge(i) ! Contains *only* merged bodies

      ! Remove the entry for body j in the given variables, shifting all lower entries up by one
      call discard_mass_reorder5(j, nbod, pbod)

    end do

    ! Set the number of mergers to zero, and deallocate the merger list
    deallocate(pmerge)

  end if

  ! Initialize the removal status flag, then make a copy of pbod and param
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i) SHARED(nbod, pbod, p, iwhy)
  do i = 1, nbod

    iwhy(i) = 0
    p(i) = pbod(i)

  end do
  !$OMP END PARALLEL DO

  ! Now check if any other bodies need to be removed because of their distance, or energetics
  if((param.rmin >= 0.0_rk) .or. (param.rmax >= 0.0_rk) .or. (param.rmaxu >= 0.0_rk)) then

    rmin2 = param.rmin**2
    rmax2 = param.rmax**2
    rmaxu2 = param.rmaxu**2

    ! Convert from heliocentric => barycentric
    call coord_h2b(param, nbod, p, msys)

    ! Restore coordinate flags to their previous values
    param.lhelio = .true.
    param.lbary = .false.
    param.lcanon = .false.

    ! For each body, check each of the three criteria for removal
    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, rh2, rb2, vb2, energy) &
    !$OMP FIRSTPRIVATE(param, time, rmin2, rmax2, rmaxu2, msys) SHARED(nbod, p, pbod, iwhy)
    do i = 2, nbod

      ! Square of heliocentric distance
      rh2 = sum(pbod(i).r**2)

      ! Too far from central body
      if((param.rmax >= 0.0_rk) .and. (rh2 > rmax2)) then

        write(*,'(a,i7,a,1pe14.6)') ' Particle ', i, ' too far from the central body at t = ', time
        write(*,'(4(1pe14.6))') pbod(i).r, sqrt(rh2)
        iwhy(i) = -3

      end if

      ! Too close to the central body
      if((param.rmin >= 0.0_rk) .and. (rh2 < rmin2)) then

        write(*,'(a,i7,a,1pe14.6)') ' Particle ', i, ' too close to the central body at t = ', time
        write(*,'(4(1pe14.6))') pbod(i).r, sqrt(rh2)
        iwhy(i) = 1

      end if

      ! If body i did not have an encounter, or wasn't removed for either of the two reasons above
      if((.not. p(i).lenc) .and. (param.rmaxu >= 0.0_rk) .and. (iwhy(i) == 0)) then

        ! Compute distance, velocity and energy with respect to the barycentre
        rb2 = sum(p(i).r**2)
        vb2 = sum(p(i).v**2)
        energy = 0.5_rk*vb2 - msys/sqrt(rb2)

        ! Unbound and too far from barycentre
        if((energy > 0.0_rk) .and. (rb2 > rmaxu2)) then

          write(*,'(a,i7,a,1pe14.6)') ' Particle ', i, ' is unbound and too far from barycentre at t = ', time
          write(*,'(5(1pe14.6))') p(i).r, sqrt(rb2), energy
          iwhy(i) = -2

        end if

      end if

    end do
    !$OMP END PARALLEL DO

  end if

  ! Check if any body should be removed if its perihelion distance is less than qmin
  if(param.qmin >= 0.0_rk) call discard_mass_peri(param, time, nbod, pbod, iwhy)

  ! If any particles have been flagged for removal
  if(any(iwhy /= 0)) then
  
    ! If we are following the energy of the system
    ! Compute the total energy of the system *before* removal
    ! This will not include the merger entries, since they have already been removed
    if(param.lenergy) call anal_energy(param, nbod, nbodm, pbod, ke, pot, ei, l)

    ! Locate each body flagged for removal
    do i = 2, nbod

      if(iwhy(i) /= 0) then

        ! Write discard information to discard_mass.dat
        call io_discard_mass(time, pbod(i), iwhy(i))

        ! Append body i to the merger list for simpler bookkeeping
        pmerge => util_append_index(pmerge, i)

      end if

    end do

    ! If we have any particles for removal
    if(associated(pmerge)) then
	
      ! Sort the merger list in ascending order
      !pmerge => util_sort_pointer(pmerge)

      ! Remove all the bodies in the merger list in reverse order
      do i = size(pmerge), 1, -1

        j = pmerge(i) ! Contains *only* discarded bodies
        call discard_mass_reorder5(j, nbod, pbod)

      end do

      ! Deallocate the merger list
      deallocate(pmerge)

    end if

  end if

  ! If any particles have been flagged for removal,
  ! and if we are following the energy of the system
  if(any(iwhy /= 0) .and. param.lenergy) then

    ! Compute the total energy of the system *after* removal
    call anal_energy(param, nbod, nbodm, pbod, ke, pot, ef, l)

    ! Update energy offset
    eoffset = eoffset + ei - ef

  end if

  return
  end subroutine discard_massive5
!!
!!
end module discard
