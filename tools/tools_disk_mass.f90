program tools_disk_mass
!
! Computes disk mass in solids
!
use swift
use util
implicit none

! Parameters
real(rk), parameter :: a_snow = 2.7_rk ! Snow line [AU]
real(rk), parameter :: f_snow = 4.2_rk ! Material enhancement beyond snow line

! Input variables
integer(ik) :: narg, command_argument_count
real(rk) :: sigma_0, a_in, a_out, p
real(rk) :: sigma, mdisk
character(len = 80) :: buffer

!-----------------!
! Executable Code !
!-----------------!

! Get the number of command line arguments
narg = command_argument_count()

! Read in the command line arguments
if(narg == 7) then

  call get_command_argument(1, buffer)
  read(buffer,*) sigma_0

  call get_command_argument(2, buffer)
  read(buffer,*) a_in

  call get_command_argument(3, buffer)
  read(buffer,*) a_out

  call get_command_argument(4, buffer)
  read(buffer,*) p

else

  write(*, '(a)', advance = 'no') 'Input the surface density @ 1.0 AU [g/cm^2]: '
  read(*,*) sigma_0

  write(*, '(a)', advance = 'no') 'Input the inner edge of the disk [AU]: '
  read(*,*) a_in

  write(*, '(a)', advance = 'no') 'Input the outer edge of the disk [AU]: '
  read(*,*) a_out

  write(*, '(a)', advance = 'no') 'Input the power-law exponent: '
  read(*,*) p

end if


!sigma_0 = 7.0_rk ! [g/cm^2]
!sigma = 3.749e-2_rk*sigma_0 ! [M_Earth/AU^2]
!a_in = 5.0_rk
!a_out = 15.0_rk
!p = 1.5_rk

! Convert sigma_0 from [g/cm^2] to [M_Earth/AU^2]
sigma = 3.749e-2_rk*sigma_0

! Compute the value of mdisk given the parameters
mdisk = util_disk_mass(sigma, a_in, a_out, p, a_snow, f_snow)

! Print values to the screen
write(*,'(/,a,1pe13.5)') "sigma_0 [g/cm^2] = ", sigma_0
write(*,'(a,1pe13.5)') "a_in [AU]        = ", a_in
write(*,'(a,1pe13.5)') "a_out [AU]       = ", a_out
write(*,'(a,1pe13.5)') "mdisk [M_Earth]  = ", mdisk

! Done
write(*,'(/,a)') '*** Done ***'

end program tools_disk_mass
