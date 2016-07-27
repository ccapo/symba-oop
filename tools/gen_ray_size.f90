program gen_ray_size
!---------------------------------------------------------------------------------
!				GEN_RAY_SIZE.F90
!---------------------------------------------------------------------------------
!
! Generates a SyMBA-style pl.in file with planetesimals distributed
! according to a power-law surface density and Rayleigh e-i distn.
! Also has option for one or more 'big guys' with orbital elements read in.
!
! May 22: This has feature forbidding
! plsmls to be initially within a sphere of radius 2 r_Hill from any 'big' guy
!
! March08: This version needs a file called 'size.dat' (generated by running
! 'drag_size_gen') with bins for distribution of drag radii
!
! April 20/2008:  Outputs file called 'dragradius.dat' with particle number,
! radius bin, and physical radius (km)
!
!---------------------------------------------------------------------------------
use swift
use util
use orbel
implicit none

! Parameters
integer(ik), parameter :: kbinmax = 1000

! Internal variables
integer(ik) :: nseed, seed(12), time(8)

integer(ik) :: nbig, nplsml, npl
integer(ik) :: i, j, iseed, ialpha
real(rk) :: amin, amax, power, p2, factor, frac
real(rk) :: rmse, rmsi, random

real(rk) :: rhopl, mpl, rhpl, rpl
real(rk) :: apl, epl, ipl, capom, omega, capm
real(rk) :: gm, r(NDIM), v(NDIM)

real(rk) :: mbig(NPLMAX), rhsqbig(NPLMAX), r2(NPLMAX), rbig(NPLMAX,NDIM)
real(rk) :: mplsml, rhplsml, rplsml, aplsml, eplsml, iplsml

integer(ik) :: nplk(kbinmax)
real(rk), dimension(kbinmax) :: mplk, rhok, rdragk, rmsek, rmsik

integer(ik) :: k, kk, kbin, kplsml

! Random seed based on current date and time
call random_seed(size = nseed)
call date_and_time(values = time)
seed = 79867_ik*sum(abs(time)) + 37_ik * (/ (k - 1, k = 1, nseed) /)
!write(*,*) 'seed = ', seed
call random_seed(put = seed)

! Start with data for the plsml distn in semimajor axis
write(*,'(/,a)', advance = 'no') 'For plsml distn, enter amin (AU), amax (AU), power law exponent: '
read(*,*) amin, amax, power
write(*,'(/,3(1pe10.3))') amin, amax, power
p2 = 2.0_rk - power

! Now data for 'big guys'
write(*,'(/,a)', advance = 'no') 'Input number of big guys and their density (g/cm^3): '
read(*,*) nbig, rhopl
write(*,'(/,i6,1pe10.3)') nbig, rhopl

! Read in data from 'size.dat'
! Note we need to figure out how many plsmls before outputting any pl.in data
nplsml = 0
open(unit = 10, file = 'size.dat', status = 'old')

read(10,*) kbin

do k = 1, kbin

  read(10,*) kk, nplk(k), mplk(k), rhok(k), rdragk(k), rmsek(k), rmsik(k)
  nplsml = nplsml + nplk(k)

end do

! Now open the output file and write out vbles, starting with sun
npl = nbig + nplsml + 1
write(*,'(/,a,1x,i7,/)') 'Total number of bodies including Sun: ', npl

open(unit = 8, file = 'pl.in', status = 'unknown')
write(8,'(1x,i7,1x,1pe22.15)') npl, MSUN

! Now add the big guys, reading in their masses and orb elements
ialpha = -1

if(nbig > 0) then

  do i = 1, nbig

    write(*,*) 'Enter mpl (M_EARTH), apl (AU), epl, ipl, capom, omega, capm (degrees)'
    read(*,*) mpl, apl, epl, ipl, capom, omega, capm
    write(*,*) mpl, apl, epl, ipl, capom, omega, capm
    rpl = 5.221e-5_rk*(3.0_rk*mpl/rhopl)**(1.0_rk/3.0_rk)
    mpl = mpl*MEARTH
    rhpl = apl*(mpl/(3.0_rk*MSUN))**(1.0_rk/3.0_rk)
    ipl = ipl/DEGRAD
    capom = capom/DEGRAD
    omega = omega/DEGRAD
    capm =  capm/DEGRAD
    gm = MSUN + mpl
    call orbel_el2xv(gm, ialpha, apl, epl, ipl, capom, omega, capm, r, v)

    write(8,'(1x,i7)')          i + 1
    write(8,'(5(1x,1pe22.15))') mpl, rhpl, rpl, 0.0_rk, rhopl
    write(8,'(3(1x,1pe22.15))') r
    write(8,'(3(1x,1pe22.15))') v
    mbig(i) = mpl
    rhsqbig(i) = rhpl*rhpl
    rbig(i,:) = r

  end do

end if

! Now write out plsml data and the rdrag data
!open(unit = 12, file = 'dragradius.dat', status = 'unknown')

kplsml = nbig + 2

do k = 1, kbin

  ! Compute radius of plsml in AU and convert mass to sim units
  rplsml = 5.221e-5_rk*(3.0_rk*mplk(k)/rhok(k))**(1.0_rk/3.0_rk)
  mplsml = mplk(k)*MEARTH
  rmse = rmsek(k)
  rmsi = rmsik(k)

  factor = 1.0_rk/max(1.0_rk, real(nplk(k) - 1, rk))

  do i = 1, nplk(k)

    frac = real(i - 1, rk)*factor
    aplsml = (amin**p2 + frac*(amax**p2 - amin**p2))**(1.0_rk/p2)
    rhplsml = aplsml*(mplsml/(3.0_rk*MSUN))**(1.0_rk/3.0_rk)
    eplsml = util_rayleigh(rmse, 0.0_rk, 1.0_rk)
    !eplsml = 1.0_rk - 30.0_rk/aplsml ! eccentricity based on a q = 30 AU
    iplsml = util_rayleigh(rmsi, 0.0_rk, pi)
    gm = MSUN + mplsml

    if(nbig > 0) then

      r2 = 0.0_rk
      do while(any(r2(1:nbig) < 4.0_rk*rhsqbig(1:nbig)))

        capom = twopi*util_randomu()
        omega = twopi*util_randomu()
        capm = twopi*util_randomu()
        call orbel_el2xv(gm, ialpha, aplsml, eplsml, iplsml, capom, omega, capm, r, v)

        ! Now check for proximity to big guys
        do j = 1, nbig

          r2(j) = sum((rbig(j,:) - r)**2)

        end do

      end do

    else

      capom = twopi*util_randomu()
      omega = twopi*util_randomu()
      capm = twopi*util_randomu()
      call orbel_el2xv(gm, ialpha, aplsml, eplsml, iplsml, capom, omega, capm, r, v)

    end if

    write(8,'(1x,i7)')          kplsml
    write(8,'(5(1x,1pe22.15))') mplsml, rhplsml, rplsml, rdragk(k), rhok(k)
    write(8,'(3(1x,1pe22.15))') r
    write(8,'(3(1x,1pe22.15))') v

    !write(12,'(1x,i7,2(1x,1pe22.15))') kplsml, rdragk(k), 0.84423_rk/(rdragk(k)*rhok(k))
    kplsml = kplsml + 1

  end do

end do

write(*,'(/,a)') 'Completed successfully'

stop
end program gen_ray_size
