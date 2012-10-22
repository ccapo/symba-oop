module util
! Module for routines contained in the util directory
use swift
use orbel

! Maximum vector size to use in interchange_sort
integer(ik), parameter :: MAX_SIMPLE_SORT_SIZE = 6

contains
!!
!!
  subroutine util_version
  !-------------------------------------------------------------------------
  !         UTIL_VERSION.F90
  !-------------------------------------------------------------------------
  ! Prints version of SWIFT and contact information
  !
  ! Authors: Hal Levison
  ! Date: 02/21/94
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------
  implicit none

  !-----------------!
  ! Executable code !
  !-----------------!

  write(*,'(a)')        "!---------------------------------------------------------!"
  write(*,'(a)')        "!                                                         !"
  write(*,'(a,f3.1,a)') "! SWIFT (Version: ", VER_NUM, ")                                    !"
  write(*,'(a)')        "!                                                         !"
  write(*,'(a)')        "!---------------------------------------------------------!"
  write(*,'(a)')        "!                                                         !"
  write(*,'(a)')        "! Authors:                                                !"
  write(*,'(a)')        "!  Martin Duncan: Queen's University                      !"
  write(*,'(a)')        "!  Hal Levison: Southwest Research Institute              !"
  write(*,'(a)')        "!                                                         !"
  write(*,'(a)')        "! Please address any comments or questions to:            !"
  write(*,'(a)')        "!  Hal Levison                                            !"
  write(*,'(a)')        "!  Geophysical, Astrophysical & Planetary Sciences        !"
  write(*,'(a)')        "!  Southwest Research Institute                           !"
  write(*,'(a)')        "!  1050 Walnut St.                                        !"
  write(*,'(a)')        "!  Suite 429                                              !"
  write(*,'(a)')        "!  Boulder, Co 80302                                      !"
  write(*,'(a)')        "!  (303) 546-0290                                         !"
  write(*,'(a)')        "!  Fax: (303) 546-9687                                    !"
  write(*,'(a)')        "!  (D)  swri::levison                                     !"
  write(*,'(a)')        "!  (I)  hal@gort.space.swri.edu                           !"
  write(*,'(a)')        "!                                                         !"
  write(*,'(a,/)')      "!---------------------------------------------------------!"

  return
  end subroutine util_version
!!
  subroutine util_exit(lflag)
  !---------------------------------------------------------------------
  !         UTIL_EXIT.F90
  !---------------------------------------------------------------------
  ! Exits program
  !
  ! Input:  lflag ==> Exit status flag
  !                   - lflag = .TRUE. if normal termination
  !                   - lflag = .FALSE. if terminating because of error
  !
  ! Authors: Hal Levison
  ! Date: 08/06/93
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !---------------------------------------------------------------------
  implicit none

  ! Input variable
  logical(lk) :: lflag

  !-----------------!
  ! Executable code !
  !-----------------!

  if(lflag) then

    write(*,'(/,a,f3.1,a)') "Normal termination of SWIFT (Version: ", VER_NUM, ")"
    write(*,'(a)')          "------------------------------------------"

  else

    write(*,'(/,a,f3.1,a)') "Terminating SWIFT (Version: ", VER_NUM, ") due to ERROR!!!"
    write(*,'(a)')          "------------------------------------------------"

  end if

  stop
  end subroutine util_exit
!!
  function util_disk_mass(sigma_0, a_in, a_out, p, a_snow, f_snow) result(m_disk)
  !--------------------------------------------------------------------------
  !         UTIL_DISK_MASS.F90
  !--------------------------------------------------------------------------
  ! Computes the disk mass in solids
  !
  ! Input:  sigma_0 ==> Surface density @ 1.0 AU [MSun/AU^2]
  !         a_in    ==> Inner radius of the disk [AU]
  !         a_out   ==> Outer radius of the disk [AU]
  !         p       ==> Power-law exponent
  !         a_snow  ==> Location of the snow line [AU]
  !         f_snow  ==> Material enhancement factor
  !
  ! Output: m_disk ==> Mass of disk in solids [MSun]
  !
  ! By: Chris Capobianco
  ! Date: 05/27/07
  !--------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: sigma_0, a_in, a_out, p, a_snow, f_snow

  ! Output variables
  real(rk) :: m_disk

  ! Internal variables
  real(rk) :: r_0, p2, a_is, a_os

  ! Stop and print warning if p >= 1.0 and a_in = 0.0
  if((p >= 1.0_rk) .and. (a_in == 0.0_rk)) then

    write(*,'(a)')       "SWIFT Error:: util_disk_mass"
    write(*,'(a)')       " Disk mass does not converge!"
    write(*,'(a,g14.6)') " Power-law = ", p
    write(*,'(a,g14.6)') " a_in (AU) = ", a_in
    call util_exit(FAILURE)

  end if

  r_0 = 1.0_rk
  p2 = 2.0_rk - p
  a_is = a_in/a_snow
  a_os = a_out/a_snow

  if(p /= 2.0_rk) then

    if((a_in < a_snow) .and. (a_snow < a_out)) then

      m_disk = twopi*sigma_0*r_0**p*a_snow**p2*(1.0_rk + f_snow*(a_os**p2 - 1.0_rk) - a_is**p2)/p2

    else ! a_in < a_snow and a_out < a_snow

      m_disk = twopi*sigma_0*r_0**p*(a_out**p2 - a_in**p2)/p2

    endif

  else

    if((a_in < a_snow) .and. (a_snow < a_out)) then

      m_disk = twopi*sigma_0*r_0**p*(f_snow*log(a_os) - log(a_is))

    else ! a_in < a_snow and a_out < a_snow

      m_disk = twopi*sigma_0*r_0**p*log(a_out/a_in)

    endif

  endif

  if((a_in > a_snow) .and. (a_out > a_snow)) m_disk = f_snow*m_disk

  return
  end function util_disk_mass
!!
  subroutine util_hills(nbod, pbod, rhill)
  !-------------------------------------------------------------------------
  !         UTIL_HILLS.F90
  !-------------------------------------------------------------------------
  ! Calculates the radius of the Hill sphere for all bodies
  !
  ! Input:  nbod        ==> Number of bodies
  !         pbod(nbod)  ==> Database entries for bodies (See swift module)
  !
  ! Output: rhill(nbod) ==> Hill sphere radius
  !
  ! Authors: Hal Levison
  ! Date: 02/19/93
  ! Last revision: 01/06/97
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized using OpenMP
  !-------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  type(body_t) :: pbod(nbod)

  ! Output variable
  real(rk) :: rhill(nbod)

  ! Internal variables
  integer(ik) :: i
  real(rk) :: mstar, mu, energy, a, r, v2

  !-----------------!
  ! Executable code !
  !-----------------!

  rhill(1) = 0.0_rk
  mstar = pbod(1).mass

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i, mu, r, v2, energy, a) &
  !$OMP FIRSTPRIVATE(mstar) SHARED(nbod, pbod, rhill)
  do i = 2, nbod

    if(pbod(i).mass > 0.0_rk) then

      mu = mstar*pbod(i).mass/(mstar + pbod(i).mass)
      r = sqrt(sum(pbod(i).r**2))
      v2 = sum(pbod(i).v**2)
      energy = -mstar*pbod(i).mass/r + 0.5_rk*mu*v2
      a = -mstar*pbod(i).mass/(2.0_rk*energy)
      rhill(i) = a*(mu/(3.0_rk*mstar))**(1.0_rk/3.0_rk)

    else

      rhill(i) = 0.0_rk

    end if

  end do
  !$OMP END PARALLEL DO

  return
  end subroutine util_hills
!!
  function util_hills_one(mstar, pbod) result(rhill)
  !-------------------------------------------------------------------------
  !         UTIL_HILLS_ONE.F90
  !-------------------------------------------------------------------------
  ! Calculates the radius of the Hill sphere for a *single* body
  !
  ! Input:  mstar ==> Mass of central body
  !         pbod  ==> Database entry for a *single* body
  !
  ! Output: rhill ==> Hill sphere radius
  !
  ! Remarks: Based on util_hill
  ! Authors: Hal Levison
  ! Date: 01/08/97
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !-------------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: mstar
  type(body_t) :: pbod

  ! Output variable
  real(rk) :: rhill

  ! Internal variables
  real(rk) :: mu, energy, a, r, v2

  !-----------------!
  ! Executable code !
  !-----------------!

  mu = mstar*pbod.mass/(mstar + pbod.mass)
  r = sqrt(sum(pbod.r**2))
  v2 = sum(pbod.v**2)
  energy = -mstar*pbod.mass/r + 0.5_rk*mu*v2
  a = -mstar*pbod.mass/(2.0_rk*energy)
  rhill = a*(mu/(3.0_rk*mstar))**(1.0_rk/3.0_rk)

  return
  end function util_hills_one
!!
  subroutine util_ir_ir3(nbod, istart, pbod, ir, ir3)
  !------------------------------------------------------------------------------
  !         UTIL_IR_IR3.F90
  !------------------------------------------------------------------------------
  ! Calculate inverse of distance (1/r) and inverse of distance cubed
  ! (1/r^3) for bodies from istart to nbod
  !
  ! Input:  nbod       ==> Number of bodies
  !         istart     ==> Index of body with which to start
  !         pbod(nbod) ==> Database entries for bodies (See swift module)
  !
  ! Output: ir(nbod)   ==> Inverse of distance, 1/r
  !         ir3(nbod)  ==> Inverse of distance cubed, 1/r^3
  !
  ! Author: Hal Levison
  ! Date: 02/02/93
  ! Last revision: 02/24/94
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !------------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, istart
  type(body_t) :: pbod(nbod)

  ! Output variables
  real(rk) :: ir(nbod), ir3(nbod)

  ! Internal variables
  integer(ik) :: i
  real(rk) :: r2

  !-----------------!
  ! Executable code !
  !-----------------!

  !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i, r2) SHARED(nbod, istart, pbod, ir, ir3)
  do i = istart, nbod

    r2 = sum(pbod(i).r**2)
    ir(i) = 1.0_rk/sqrt(r2)
    ir3(i) = ir(i)/r2

  end do
  !$OMP END PARALLEL DO

  return
  end subroutine util_ir_ir3
!!
  function util_kahan_sum(xsum_current, xi, xerror) result(xsum_new)
  !-------------------------------------------------------------------------------------------
  !         UTIL_KAHAN_SUM.F90
  !-------------------------------------------------------------------------------------------
  ! Sums two floating point scalars more accurately utilitizing the Kahan summation formula
  ! This function is designed to be used inside a do or while loop, where the initial value of
  ! of xsum_current is initialized appropriately and the initial value of xerror is 0.0_rk
  !
  ! N.B. Use this function if the summation is being performed for more than *three* terms
  !
  ! Input:  xsum_current - Current value of the sum
  !         xi           - i-th term to be added to the sum
  !         xerror       - Error term from the previous term of the sum
  !
  ! Output: xsum_new     - The updated value of the sum
  !         xerror       - The error term for this term of the sum
  !
  ! By: Chris Capobianco
  ! Date: 05/04/09
  !-------------------------------------------------------------------------------------------
  implicit none

  ! Input/Output variables
  real(rk), intent(in) :: xsum_current, xi
  real(rk), intent(inout) :: xerror
  real(rk) :: xsum_new

  ! Internal variables
  real(rk) :: low_bits

  low_bits = xi - xerror
  xsum_new = xsum_current + low_bits
  xerror = (xsum_new - xsum_current) - low_bits

  return
  end function util_kahan_sum
!!
  subroutine util_mass_peri(nbod, pbod, peri, lperi)
  !-------------------------------------------------------------------------
  !         UTIL_MASS_PERI.F90
  !-------------------------------------------------------------------------
  ! Determines whether perihelion passage of a body has taken place
  !
  ! Input:  nbod        ==> Number of bodies
  !         pbod(nbod)  ==> Database entries for bodies (See swift module)
  !
  ! Output: pbod.iperi  ==> Pericentre status
  !                         - iperi = -1 if body before perihelion
  !                         - iperi =  0 if body went through perihelion
  !                         - iperi = +1 if body after perihelion
  !         peri(nbod)  ==> Pericentre distance if iperi = 0
  !         lperi(nbod) ==> Pericentre flag
  !                         - lperi = .TRUE. if iperi = 0
  !
  ! Remarks: Based on util_peri.f
  ! Authors: Hal Levison
  ! Date: 12/30/96
  ! Last revision: 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !                             - Parallelized code using OpenMP
  !-------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  type(body_t) :: pbod(nbod)

  ! Output variables
  logical(lk) :: lperi(nbod)
  real(rk) :: peri(nbod)

  ! Internal variables
  logical(lk), save :: lfirst = .true.
  integer(ik) :: i, ialpha
  real(rk) :: mstar, gm, a, e, vdotr

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Enter here if we are setting things up
  if(lfirst) then

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, vdotr) SHARED(nbod, pbod)
    do i = 2, nbod

      vdotr = dot_product(pbod(i).r, pbod(i).v)

      if(vdotr > 0.0_rk) then

        pbod(i).iperi = 1

      else

        pbod(i).iperi = -1

      end if

    end do
    !$OMP END PARALLEL DO

    lfirst = .false. ! So we don't enter here again

  else

    mstar = pbod(1).mass

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i, vdotr, gm, ialpha, a, e) &
    !$OMP FIRSTPRIVATE(mstar) SHARED(nbod, pbod, lperi, peri)
    do i = 2, nbod

      lperi(i) = .false.           ! Initialize perihelion flag

      vdotr = dot_product(pbod(i).r, pbod(i).v)

      if(pbod(i).iperi == -1) then ! was coming in

        if(vdotr < 0.0_rk) then    ! still coming in

          pbod(i).iperi = -1

         else                      ! turned around

          pbod(i).iperi = 0
          lperi(i) = .true.
          gm = mstar + pbod(i).mass
          call orbel_xv2aeq(pbod(i).r, pbod(i).v, gm, ialpha, a, e, peri(i))

        end if

      else

        if(vdotr < 0.0_rk) then    ! coming in

          pbod(i).iperi = -1

        else

          pbod(i).iperi = 1        ! going out

        end if

      end if

    end do
    !$OMP END PARALLEL DO

  end if

  return
  end subroutine util_mass_peri
!!
  function util_randomu() result(ran)
  !--------------------------------------------------------------
  !         UTIL_RANDOMU.F90
  !--------------------------------------------------------------
  ! Generates a uniform random number between (0.0, 1.0) using
  ! Intel Fortran 90 random number generator
  !--------------------------------------------------------------
  implicit none

  ! Output variable
  real(rk) :: ran

  !-----------------!
  ! Executable code !
  !-----------------!

  call random_number(ran)

  return
  end function util_randomu
!!
  function util_power_law(p, xmin, xmax) result(x)
  !------------------------------------------------------------------
  !         UTIL_POWER_LAW.F90
  !------------------------------------------------------------------
  ! Returns a random number drawn from a power-law distribution
  ! with exponent p, where xmin and xmax are the minimum and
  ! maximum range for the random number.
  !
  ! Input:  p    ==> Exponent of power-law distribution
  !         xmin ==> Minimum value for random number
  !         xmax ==> Maximum value for random number
  !
  ! Output: x    ==> Random number drawn from power-law distrubution
  !
  ! By: Chris Capobianco
  ! Date: 02/04/09
  !------------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: p, xmin, xmax

  ! Output variable
  real(rk) :: x

  ! Internal variables
  real(rk) :: p1, ip1, xi

  !-----------------!
  ! Executable code !
  !-----------------!

  p1 = 1.0_rk - p
  ip1 = 1.0_rk/p1
  xi = util_randomu()

  if(p /= 1.0_rk) then

    x = (xmin**p1 + xi*(xmax**p1 - xmin**p1))**ip1

  else

    x = xmin*(xmax/xmin)**xi

  end if

  return
  end function util_power_law
!!
  function util_rayleigh(rms, xmin, xmax) result(x)
  !--------------------------------------------------------------
  !         UTIL_RAYLEIGH.F90
  !--------------------------------------------------------------
  ! Returns a random number drawn from a Rayleigh distribution
  ! with dispersion rms, where xmin and xmax are the minimum and
  ! maximum range for the random number.
  !
  ! Input:  rms  ==> RMS of Rayleigh distribution
  !         xmin ==> Minimum value for random number
  !         xmax ==> Maximum value for random number
  !
  ! Output: x    ==> Random number drawn from Rayleigh distrubution
  !
  ! By: Chris Capobianco
  ! Date: 02/04/09
  !--------------------------------------------------------------
  implicit none

  ! Input variables
  real(rk) :: rms, xmin, xmax

  ! Output variable
  real(rk) :: x

  ! Internal variables
  real(rk) :: fmin, fmax, xi

  !-----------------!
  ! Executable code !
  !-----------------!

  fmin = -exp(-0.5_rk*(xmin/rms)**2)
  fmax = -exp(-0.5_rk*(xmax/rms)**2)
  xi = util_randomu()

  x = rms*sqrt(-2.0_rk*log((fmin - fmax)*xi - fmin))

  return
  end function util_rayleigh
!!
  function util_append_index(pold, iappend) result(pnew)
  !-------------------------------------------------------------------
  !         UTIL_APPEND_INDEX.F90
  !-------------------------------------------------------------------
  ! Appends the integer iappend to the end of an integer pointer list
  !
  ! Input:  pold    ==> Current integer pointer list
  !         iappend ==> Index to append to pold
  !
  ! Output: pnew    ==> Updated integer pointer list
  !
  ! By: Chris Capobianco
  ! Date: 09/04/08
  !-------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik), pointer :: pold(:)
  integer(ik), intent(in) :: iappend

  ! Output variable
  integer(ik), pointer :: pnew(:)

  ! Internal variable
  integer(ik) :: nold, ierr

  !-----------------!
  ! Executable code !
  !-----------------!

  ! If pold is already associated
  if(associated(pold)) then

    ! Number of entries in pold
    nold = size(pold)

    ! Create a new pointer list that is one entry larger than pold
    allocate(pnew(nold + 1), stat = ierr)
    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: util_append_index"
      write(*,'(a)') " Cannot allocate memory to new pointer"
      call util_exit(FAILURE)

    end if

    ! Copy contents of pold into pnew
    pnew(1:nold) = pold

    ! Add iappend to final entry
    pnew(nold + 1) = iappend

    ! Destroy pold if associated
    deallocate(pold)

  else

    ! Create a pointer list that has one entry
    allocate(pnew(1), stat = ierr)
    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: util_append_index"
      write(*,'(a)') " Cannot allocate memory to new pointer"
      call util_exit(FAILURE)

    end if

    ! Put iappend in first entry
    pnew(1) = iappend

  end if

  return
  end function util_append_index
!!
  function util_remove_index(pold, iremove) result(pnew)
  !-------------------------------------------------------------------
  !         UTIL_REMOVE_INDEX.F90
  !-------------------------------------------------------------------
  ! Appends the integer iappend to the end of an integer pointer list
  !
  ! Input:  pold    ==> Current integer pointer list
  !         iappend ==> Index to append to pold
  !
  ! Output: pnew    ==> Updated integer pointer list
  !
  ! By: Chris Capobianco
  ! Date: 09/04/08
  !-------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik), pointer :: pold(:)
  integer(ik), intent(in) :: iremove

  ! Output variable
  integer(ik), pointer :: pnew(:)

  ! Internal variables
  integer(ik) :: nold, ierr

  ! Check if pold is associated
  if(associated(pold)) then

    ! Number of entries in pold
    nold = size(pold)

    ! Create a new pointer list that is one element smaller than pold
    if(nold > 1) then

      allocate(pnew(nold - 1), stat = ierr)
      if(ierr /= 0) then

        write(*,'(a)') "SWIFT Error:: util_remove_index"
        write(*,'(a)') " Cannot allocate memory to new pointer"
        call util_exit(FAILURE)

      end if

      if(iremove == 1) then

        ! Copy the remaining entries into pnew
        pnew = pold(2:nold)

      else if(iremove == nold) then

        ! Copy the remaining entries into pnew
        pnew = pold(1:nold - 1)

      else

        ! Copy the remaining entries into pnew
        pnew = (/ pold(1:iremove - 1), pold(iremove + 1:nold) /)

      end if

    else

      ! Set pnew as NULL() since we have removed the final entry
      deallocate(pnew)

    end if

    ! Destroy pold
    deallocate(pold)

  end if

  return
  end function util_remove_index
!!
  recursive subroutine util_sort_pointer_recursive(left_end, right_end, p)
  !-----------------------------------------------------------------
  !         UTIL_SORT_POINTER_RECURSIVE.F90
  !-----------------------------------------------------------------
  ! Recursive sorting routine for integer pointers
  !-----------------------------------------------------------------
  implicit none
  
  ! Input variables
  integer(ik), intent(in) :: left_end, right_end
  
  ! Input and Output variable
  integer(ik), pointer :: p(:)
  
  ! Internal variables
  integer(ik) :: i, j
  integer(ik) :: reference

  !-----------------!
  ! Executable code !
  !-----------------!

  if(right_end < left_end + MAX_SIMPLE_SORT_SIZE) then

    call util_interchange_sort_pointer(left_end, right_end, p)

  else

    reference = p((left_end + right_end)/2)
    i = left_end - 1; j = right_end + 1

    do

      do

        i = i + 1
        if(p(i) >= reference) exit

      end do

      do

        j = j - 1
        if(p(j) <= reference) exit

      end do

      if(i < j) then

        call util_swap_integer(p(i), p(j))

      else if(i == j) then

        i = i + 1
        exit

      else

        exit

      end if

    end do

    if(left_end < j) call util_sort_pointer_recursive(left_end, j, p)
    if(i < right_end) call util_sort_pointer_recursive(i, right_end, p)

  end if

  end subroutine util_sort_pointer_recursive
!!
  recursive subroutine util_sort_body_recursive(left_end, right_end, p)
  !-----------------------------------------------------------------
  !         UTIL_SORT_BODY_RECURSIVE.F90
  !-----------------------------------------------------------------
  ! Recursive sorting routine for body entries
  !-----------------------------------------------------------------
  implicit none
  
  ! Input variables
  integer(ik), intent(in) :: left_end, right_end
  
  ! Input and Output variable
  type(body_t) :: p(:)
  
  ! Internal variables
  integer(ik) :: i, j
  real(rk) :: reference

  !-----------------!
  ! Executable code !
  !-----------------!

  if(right_end < left_end + MAX_SIMPLE_SORT_SIZE) then

    call util_interchange_sort_body(left_end, right_end, p)

  else

    reference = p((left_end + right_end)/2).mass
    i = left_end - 1; j = right_end + 1

    do

      do

        i = i + 1
        if(p(i).mass <= reference) exit

      end do

      do

        j = j - 1
        if(p(j).mass >= reference) exit

      end do

      if(i < j) then

        call util_swap_body(p(i), p(j))

      else if(i == j) then

        i = i + 1
        exit

      else

        exit

      end if

    end do

    if(left_end < j) call util_sort_body_recursive(left_end, j, p)
    if(i < right_end) call util_sort_body_recursive(i, right_end, p)

  end if

  end subroutine util_sort_body_recursive
!!
  subroutine util_interchange_sort_pointer(left_end, right_end, p)
  !---------------------------------------------------------------
  !         UTIL_INTERCHANGE_SORT_POINTER.F90
  !---------------------------------------------------------------
  ! Interchange sort for pointers
  !---------------------------------------------------------------
  implicit none
  
  ! Input variables
  integer(ik), intent(in) :: left_end, right_end
  
  ! Input and Output variable
  integer(ik), pointer :: p(:)
  
  ! Internal variables
  integer(ik) :: i, j

  !-----------------!
  ! Executable code !
  !-----------------!

  do i = left_end, right_end - 1

    do j = i + 1, right_end

      if(p(i) > p(j)) call util_swap_integer(p(i), p(j))

    end do

  end do

  end subroutine util_interchange_sort_pointer
!!
  subroutine util_interchange_sort_body(left_end, right_end, p)
  !---------------------------------------------------------------
  !         UTIL_INTERCHANGE_SORT_BODY.F90
  !---------------------------------------------------------------
  ! Interchange sort for particles
  !---------------------------------------------------------------
  implicit none
  
  ! Input variables
  integer(ik), intent(in) :: left_end, right_end
  
  ! Input and Output variable
  type(body_t) :: p(:)
  
  ! Internal variables
  integer(ik) :: i, j

  !-----------------!
  ! Executable code !
  !-----------------!

  do i = left_end, right_end - 1

    do j = i + 1, right_end

      if(p(j).mass > p(i).mass) call util_swap_body(p(i), p(j))

    end do

  end do

  end subroutine util_interchange_sort_body
!!
  subroutine util_swap_integer(a, b)
  !---------------------------------------------------------------
  !         UTIL_SWAP_INTEGER.F90
  !---------------------------------------------------------------
  ! Swaps the integers a and b
  !---------------------------------------------------------------
  implicit none
  
  ! Input and Output variables
  integer(ik), intent(inout) :: a, b
  
  ! Internal variables
  integer(ik) :: temp

  !-----------------!
  ! Executable code !
  !-----------------!

  temp = a
  a = b
  b = temp

  return
  end subroutine util_swap_integer
!!
  subroutine util_swap_body(pa, pb)
  !---------------------------------------------------------------
  !         UTIL_SWAP_BODY.F90
  !---------------------------------------------------------------
  ! Swaps the particles pa and pb
  !---------------------------------------------------------------
  implicit none
  
  ! Input and Output variables
  type(body_t), intent(inout) :: pa, pb
  
  ! Internal variable
  type(body_t) :: ptmp

  !-----------------!
  ! Executable code !
  !-----------------!

  ptmp = pa
  pa = pb
  pb = ptmp

  return
  end subroutine util_swap_body
!!
  function util_sort_pointer(pold) result(pnew)
  !---------------------------------------------------------------
  !         UTIL_SORT_POINTER.F90
  !---------------------------------------------------------------
  ! Sorts the integer pointer p in place
  !
  ! Input:  pold ==> Unsorted integer pointer list
  !
  ! Output: pnew ==> Sorted integer pointer list
  !
  ! By: Chris Capobianco
  ! Date: 02/04/09
  !---------------------------------------------------------------
  implicit none

  ! Input variable
  integer(ik), pointer :: pold(:)

  ! Output variable
  integer(ik), pointer :: pnew(:)

  ! Internal variable
  integer(ik) :: nold, ierr

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Check if pointer pold is even associated
  if(associated(pold)) then

    ! Number of entries in pold
    nold = size(pold)

    ! Create a new pointer list that has the same number of entries as pold
    allocate(pnew(nold), stat = ierr)
    if(ierr /= 0) then

      write(*,'(a)') "SWIFT Error:: util_sort_pointer"
      write(*,'(a)') " Cannot allocate memory for new pointer"
      call util_exit(FAILURE)

    end if

    ! Copy contents of pold into pnew
    pnew = pold

    ! Destroy pold
    deallocate(pold)

  else

    write(*,'(a)') "SWIFT Error:: util_sort_pointer"
    write(*,'(a)') " Input pointer not associated"
    call util_exit(FAILURE)

  end if

  ! Sort pointer pnew
  call util_sort_pointer_recursive(1, nold, pnew)

  end function util_sort_pointer
!!
  subroutine util_sort_body(nbod, nbodm, pbod)
  !---------------------------------------------------------------
  !         UTIL_SORT_body.F90
  !---------------------------------------------------------------
  ! Sorts the particle entries in order of decreasing mass, but
  ! *only* for entries above nbodm
  !
  ! Input:  pbod ==> Unsorted integer pointer list
  !
  ! Output: pnew ==> Sorted integer pointer list
  !
  ! By: Chris Capobianco
  ! Date: 02/04/09
  !---------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod, nbodm

  ! Input and Output variable
  type(body_t) :: pbod(nbod)

  !-----------------!
  ! Executable code !
  !-----------------!

  ! Sort particle entries pbod
  call util_sort_body_recursive(1, nbodm, pbod)

  end subroutine util_sort_body
!!
  subroutine util_nbodm(param, nbod, nbodm, pbod)
  !----------------------------------------------------------------------------
  !         UTIL_NBODM.F90
  !----------------------------------------------------------------------------
  ! Returns the location of the last massive body greater than mtiny
  !
  ! Input:  param      ==> Global parameters (See swift module & io_init_param)
  !         nbod       ==> Number of bodies
  !         pbod(nbod) ==> Database entries for bodies (See swift module)
  !
  ! Output: nbodm      ==> Index of the last massive body
  !
  ! Remarks: If all the objects are massive, then nbodm = nbod - 1 so that
  !          the do loops will have the correct limits
  !
  ! Authors: Hal Levison
  ! Date: 03/20/97
  ! Last revision: 01/29/06
  !              : 02/09/09 CCC - Converted to Fortran 90/95 syntax
  !----------------------------------------------------------------------------
  implicit none

  ! Input variables
  integer(ik) :: nbod
  type(body_t) :: pbod(nbod)
  type(param_t), intent(in) :: param

  ! Output variable
  integer(ik) :: nbodm

  ! Internal variables
  integer(ik) :: i
  real(rk) :: mtiny

  !-----------------!
  ! Executable code !
  !-----------------!

  ! This assumes that the particle identifiers are in order of decreasing mass
  !nbodm = count(pbod(1:nbod - 1).mass > param.mtiny)

  if(pbod(nbod).mass > param.mtiny) then

    nbodm = nbod - 1

    ! Sort the massive bodies in order of decreasing mass
    call util_sort_body(nbod, nbodm, pbod)

    write(*,'(a,i6,a)') " Out of ", nbod, " bodies, all are massive."

  else

    ! Initialize the value of nbodm
    nbodm = 1

    ! Make a private copy of mtiny available for each thread
    mtiny = param.mtiny

    !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i) &
    !$OMP FIRSTPRIVATE(mtiny) SHARED(nbod, pbod) REDUCTION(MAX : nbodm)
    do i = 2, nbod - 1

      if(pbod(i).mass > mtiny) nbodm = i

    end do
    !$OMP END PARALLEL DO

    ! Sort the massive bodies in order of decreasing mass
    call util_sort_body(nbod, nbodm, pbod)

    write(*,'(2(a,i6),a)') " Out of ", nbod, " bodies, ", nbodm, " are massive."

  end if

  return
  end subroutine util_nbodm
!!
!!
end module util
