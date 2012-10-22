program drag_size_gen
!---------------------------------------------------------------------------------------
!				DRAG_SIZE_GEN.F90
!---------------------------------------------------------------------------------------
!
! DRAG_SIZE_GEN is used to generate a file size.dat that
! contains input data to gen_ray_size.
! DRAG_SIZE_GEN generates relevant data for particle bins with a specified
! power-law distribution in planetesimal drag size 'r' s.t.
! dN = const r^(-q) dr between rdragmin and rdragmax
! in 'k' logarithmically distributed bins. The drag radius of each
! particle in a given bin is assigned to be the geometric mean of
! the bin's inner and outer radius.
! The mass in each bin is normalized so that sum of masses
! in bins equals prescribed value mtot.
!
! The desired number of superparticles in each bin is read in, along
! with the physical density of the real particle and
! rms ecc and inc for the real particles. The superparticle mass and superparticle
! radius are then computed and the relevant data for each bin is written to the file
! 'size.dat'. NOTE that the superparticle radius is needed and used by the original
! SYMBA routines (to compute the necessity and consequences of a merger for example).
! The plsml size assumed for gas drag calculations is called rdrag.
!
!---------------------------------------------------------------------------------------
use swift
implicit none

integer(ik), parameter :: kbinmax = 1000
integer(ik) :: i, kbin, nplk, n, nplsml,j
real(rk) :: rdragmin, rdragmax, q, rat
real(rk) :: rdrag(kbinmax), mbin(kbinmax), rhopl(kbinmax)
real(rk) :: mtot, mbintot
real(rk) :: rho, rmse, rmsi

! Input parameters
write(*,'(/,a)', advance = 'no') 'Enter rdragmin, rdragmax, differential index q: '
read(*,*) rdragmin, rdragmax, q
write(*,'(/,a,1pe10.3,a,1pe10.3,a,1pe9.2)') 'rdragmin = ', rdragmin, ',  rdragmax = ', rdragmax, ',  q = ', q

write(*,'(/,a)', advance = 'no') 'Enter number of size bins and summed total mass (Earth masses): '
read(*,*) kbin, mtot
write(*,'(/,a,i2,a,1pe10.3,/)') 'no. of bins = ', kbin, ', mtot = ', mtot

! Now assign drag radii and fractional masses to bins
rat = (rdragmax/rdragmin)**(1.0_rk/real(kbin - 1, rk))
mbintot = 0.0_rk

rdrag(1) = rdragmin
mbin(1) = rdrag(1)**(4.0_rk - q)
mbintot = mbintot + mbin(1)

do i = 2, kbin

  rdrag(i) = rat*rdrag(i - 1)
  mbin(i) = rdrag(i)**(4.0_rk - q)
  mbintot = mbintot + mbin(i)

end do

do i = 1, kbin

  mbin(i) = mtot*mbin(i)/mbintot
  write(*,'(a,i2,a,1pe10.3,a,1pe10.3)') 'i = ', i, ', rdrag(i) = ', rdrag(i), ', mbin(i) = ', mbin(i)

end do

! Now open output file for binned data and fill it
! by reading how many super-plsmls desired in each bin
write(*,'(/,a)') 'Enter sequence with bin number, number of plsmls in bin, '
write(*,'(a)')   'true physical density, rms ecc and rms i'

! Open size.dat
open(unit = 10, file = 'size.dat', status = 'unknown')
write(10,'(i5)') kbin

do i = 1, kbin

  read(*,*) j, nplk, rho, rmse, rmsi
  write(*,'(i6,1x,i6,1x,3(1x,1pe10.3))') j, nplk, rho, rmse, rmsi
  write(10,'(i6,1x,i6,5(1x,1pe10.3))') i, nplk, mbin(i)/real(nplk, rk), rho, rdrag(i), rmse, rmsi

end do

! Close size.dat
close(unit = 10)

end program drag_size_gen
