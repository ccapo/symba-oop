program gen_size_dist
!--------------------------------------------------------------
! Generates particle bins with a specified power-law
! distribution in planetesimal size: dN/dr = N_0(r_0/r)^q
! between rmin and rmax in nr logarithmically distributed bins
!
! The number of particles in each bin is scaled to match its
! fraction of the total times the desired number of particles
!
! The mass in each bin is scaled to match its fraction of the
! total times the desired mass of the disk
!--------------------------------------------------------------
implicit none
integer, parameter :: ik = selected_int_kind(9)
integer, parameter :: rk = selected_real_kind(15)
integer(ik) :: arg_num, iargc
integer(ik) :: i, nr, ntotal
real(rk) :: rmin, rmax, q, iq, membryo
real(rk) :: disk_mass, dn
real(rk), dimension(:), allocatable :: r, mass, fraction_mass, npart, npart_old, part_mass, part_mass_old
character(len = 80) :: buf

! Get the number of command line arguments
arg_num = iargc()

! Read in the command line arguments
if(arg_num == 7) then

  call getarg(1, buf)
  read(buf,*) nr
  call getarg(2, buf)
  read(buf,*) rmin
  call getarg(3, buf)
  read(buf,*) rmax
  call getarg(4, buf)
  read(buf,*) q
  call getarg(5, buf)
  read(buf,*) disk_mass
  call getarg(6, buf)
  read(buf,*) membryo
  call getarg(7, buf)
  read(buf,*) ntotal

else

  write(*,'(a)', advance = 'no') 'Enter number of size bins, minimum size, maximum size, and differential size exponent: '
  read(*,*) nr, rmin, rmax, q

  write(*,'(a)', advance = 'no') 'Enter the desired disk mass, the embryo mass, and total number of particles: '
  read(*,*) disk_mass, membryo, ntotal

end if

! Print summary of input parameters
write(*,'(/,a)') 'Input parameter summary:'
write(*,'(a,i2,3(a,1pg12.5))') 'nr = ', nr, '; rmin = ', rmin, ' km ; rmax = ', rmax, ' km ; q = ', q
write(*,'(2(a,1pg12.5),a,i7,/)') 'disk_mass = ', disk_mass, ' M_Earth ; membryo = ', membryo, ' M_Earth ; ntotal = ', ntotal

! Allocate vectors
allocate(r(nr), mass(nr), npart(nr), npart_old(nr), part_mass(nr), part_mass_old(nr), fraction_mass(nr))

! Define size of each bin
r = rmin*(rmax/rmin)**(/ ( real(i - 1, rk)/real(nr - 1, rk), i = 1, nr ) /)

! Compute the number of particles in each size bin, w.r.t. one particle at r = rmax
npart = (rmax/r)**q

! Define the mass in each bin
mass = npart*r**3

! Compute the fraction of mass in the bins
fraction_mass = mass/sum(mass)

! Scale the mass in each bin to disk_mass times its mass fraction
mass = disk_mass*fraction_mass

! Put an equal number of particles in each size bin
npart = ntotal/real(nr, rk)
npart_old = npart

! Compute the mass of each particle in each size bin
part_mass = mass/npart
part_mass_old = part_mass

! Make sure the mass of each particle does not exceed membryo/750.0
where(part_mass >= membryo/750.0_rk)
  part_mass = membryo/750.0_rk
  npart = mass/part_mass
end where

! To conserve particle number, remove extra particles from smaller particle bins
do i = nr, 2, -1
  if(part_mass(i) /= part_mass_old(i)) then
    dn = abs(npart_old(i) - npart(i))
    npart(1:i - 1) = npart(1:i - 1) - dn/real(i - 1, rk)
  end if
end do

! Ensure the total number of particles matches desired total
if(nint(sum(npart)) /= ntotal) then
  dn = ntotal - nint(sum(npart))
  npart = npart + sign(dn, dn)/real(nr, rk)
end if

! Adjust mass per particle to reflect the new distribution
where(npart > 0.0_rk) part_mass = mass/npart

! In case something undesirable happens
if(any(npart <= 0.0_rk)) then
  write(*,'(a)') '*** Error! ***'
  write(*,'(a)') 'There are one or more bins with no particles!'
  write(*,'(a,/)') 'Try increasing ntotal, decreasing disk_mass or the number of size bins'
  stop
end if

! Print warning message
if(any(npart <= 1.0e3_rk)) then
  write(*,'(a)') '*** Warning! ***'
  write(*,'(a)') 'There are less than 1000 particles in one or more bins'
  write(*,'(a)') 'There may not be enough particles in these bins to be statistically significant'
  write(*,'(a,/)') 'Try increasing ntotal, decreasing disk_mass or the number of size bins'
end if

! Print results
write(*,'(a)') 'Size Bin (km); Particles per Bin; Mass per Bin (M_Earth); Mass per Particle (M_Earth)'
do i = 1, nr
  write(*,'(1pe12.5,5x,i7,2(12x,1pe12.5))') r(i), nint(npart(i)), mass(i), part_mass(i)
end do
write(*,'(/,2(a,1pg12.5))') 'Sum of particles = ', nint(sum(npart)), '; Sum of masses = ', sum(npart*part_mass)

! Deallocate vectors
deallocate(r, mass, npart, npart_old, part_mass, part_mass_old, fraction_mass)

! *** Martin's code *** !
!rpl(1) = rpl_start
!mpl(1) = rpl(1)**(4.0d0 - q)
!mpltot = mpl(1)

!if(nbins > 1) then
!  ratio = (rpl_end/rpl_start)**(1.0d0/dble(nbins - 1))
!  do i = 2, nbins
!    rpl(i) = rpl(i - 1)*ratio
!    mpl(i) = rpl(i)**(4.0d0 - q)
!    mpltot = mpltot + mpl(i)
!  end do
!end if

!do i = 1, kbins
!  mpl(i) = mtot*mpl(i)/mpltot
!  write(*,'(a,i2,a,1pe10.3,a,1pe10.3)') 'i = ', i, ', rpl(i) = ', rpl(i), ', mpl(i) = ', mpl(i)
!end do
! *** Martin's code *** !

stop
end program gen_size_dist