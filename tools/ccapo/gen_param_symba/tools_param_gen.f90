	program tools_param_gen
	! param_gen
	implicit none
	integer :: i, iargc, narg
	real :: apl, mpl, mdisk, fgs, logr
	real :: r, amin, amax, rho_g, sigma_s
	character(len = 80) :: arg

	! Check for correct number of arguments
	narg = iargc()
	if(narg /= 5) then
	  write(*,'(a)') "Incorrect number of arguments: apl, mpl, mdisk, fgs, logr"
	  stop
	end if

	! Extract each argument
	call getarg(1,arg)
	read(arg,*) apl
	call getarg(2,arg)
	read(arg,*) mpl
	call getarg(3,arg)
	read(arg,*) mdisk
	call getarg(4,arg)
	read(arg,*) fgs
	call getarg(5,arg)
	read(arg,*) logr

	! Compute dependent quantities
	r = 10.0**logr
	amin = apl*(1.0 - 7.0*0.01*mpl**(1.0/3.0))
	amax = apl*(1.0 + 7.0*0.01*mpl**(1.0/3.0))
	rho_g = 1.4e-9*mdisk
	sigma_s = 7.0*mdisk/fgs

	! Print input arguments and computed values
	write(*,'(a)') "Input arguments:"
	write(*,'(a)')       "  apl       mpl       mdisk     fgs       logr"
	write(*,'(5(1pe10.2),/)') apl, mpl, mdisk, fgs, logr
	write(*,'(a)') "Computed quantities:"
	write(*,'(a)')       "  rplsml    amin      amax      rho_g     sigma_s"
	write(*,'(5(1pe10.2))') r, amin, amax, rho_g, sigma_s

	open(unit = 10, file = "param.in", status = "unknown")

	! Simulation parameters
	write(10,'(a)') "&sim_params"
	write(10,'(a)') " t0          = 0.0E0        ! Initial time (years)"
	write(10,'(a)') ",tstop       = 2.0E4        ! Final time (years)"
	write(10,'(a)') ",dt          = 0.5E0        ! Time step (years)"
	write(10,'(a)') ",dtout       = 1.0E2        ! Time interval between outputs (years)"
	write(10,'(a)') ",dtdump      = 1.0E2        ! Time interval between dumps (years)"
	write(10,'(a)') ",outfile     = xdr.dat      ! binary output filename"
	write(10,'(a)') ",outtype     = XDR          ! Output format [ REAL | XDR ]"
	write(10,'(a)') ",fopenstat   = UNKNOWN      ! Output status [ NEW | UNKNOWN | APPEND ]"
	write(10,'(a)') ",lenergy     = F            ! Flag to compute the system's energy and angular momentum [ T | F ]"
	write(10,'(a)') ",ljacobi     = F            ! Flag to compute the system's Jacobi integral [ T | F ]"
	write(10,'(a,/)') "/"

	! Stellar parameters
	write(10,'(a)') "&stellar_params"
	write(10,'(a)') " loblate     = F            ! Flag to include the J2 and J4 stellar oblateness terms [ T | F ]"
	write(10,'(a)') ",j2          = 0.0          ! Stellar J2 oblateness term"
	write(10,'(a)') ",j4          = 0.0          ! Stellar J4 oblateness term"
	write(10,'(a)') ",lclose      = T            ! Flag to check for planetary close encounters [ T | F ]"
	write(10,'(a)') ",rmin        = 3.0E0        ! Remove particle if r < rmin  (AU) (Too close to star, set rmin < 0 to ignore)"
	write(10,'(a)') ",rmax        = 1.0E3        ! Remove particle if r > rmax  (AU) (Too distant from star, set rmax < 0 to ignore)"
	write(10,'(a)') ",rmaxu       = 1.0E3        ! Remove particle if r > rmaxu (AU) (Particle unbound, set rmaxu < 0 to ignore)"
	write(10,'(a)') ",qmin        = -1.0         ! Remove particle if q < qmin  (AU) (q too close to star, set qmin < 0 to ignore)"
	write(10,'(a,/)') "/"

	! Symba parameters
	write(10,'(a)') "&symba_params"
	write(10,'(a)') "mtiny = 0.1                 ! Smallest mass to self-gravitate (M_Earth)"
	write(10,'(a,/)') "/"

	! Planetary parameters
	write(10,'(a)') "&pl_params"
	write(10,'(a)') " npl         = 1            ! Number of embryos"
	write(10,'(a)') ",rhopl       = 1.0          ! Embryo mass density (g/cm^3)"
	write(10,'(a,1pe8.2,a)') ",mpl         = ", mpl, "     ! Embryo mass (M_Earth)"
	write(10,'(a,1pe8.2,a)') ",apl         = ", apl, "     ! Embryo semimajor axis (AU)"
	write(10,'(a)') ",epl         = 0.0          ! Embryo eccentricity"
	write(10,'(a)') ",wpl         = 46.0         ! Embryo longitude of pericentre (degrees)"
	write(10,'(a)') ",ipl         = 0.0          ! Embryo inclination (degrees)"
	write(10,'(a)') ",opl         = 72.0         ! Embryo longitude of ascending node (degrees)"
	write(10,'(a)') ",mapl        = 56.0         ! Embryo mean anomaly (degrees)"
	write(10,'(a,/)') "/"

	! Model parameters
	write(10,'(a)') "&model_params"
	write(10,'(a)') " ms_model    = CONSTANT     ! Planetesimal radius model [ CONSTANT | POWER ]"
	write(10,'(a)') ",ei_model    = RAYLEIGH     ! Planetesimal eccentricity and inclination model [ CONSTANT | POWER | RAYLEIGH ]"
	write(10,'(a)') ",lplsml      = T            ! Flag to include planetesimals [ T | F ]"
	write(10,'(a)') ",lgas        = T            ! Flag to include aerodynamic and Type-I gas drag [ T | F ]"
	write(10,'(a)') ",lshear      = T            ! Flag to ensure planetesimals are in the shear-dominated regime [ T | F ]"
	write(10,'(a,/)') "/"

	! Planetesimal semi-major axis parameters
	write(10,'(a)') "&aplsml_params"
	write(10,'(a,1pe8.2,a)') " amin        = ", amin, "     ! Planetesimal minimum semimajor axis (AU)"
	write(10,'(a,1pe8.2,a)') ",amax        = ", amax, "     ! Planetesimal maximum semimajor axis (AU)"
	write(10,'(a)') ",apower      = 1.5          ! Planetesimal surface number density power-law exponent"
	write(10,'(a,/)') "/"

	! Planetesimal mass parameters
	write(10,'(a)') "&msplsml_params"
	write(10,'(a)') " smin        = -1.0         ! Planetesimal minimum radius (km)"
	write(10,'(a)') ",smax        = -1.0         ! Planetesimal maximum radius (km)"
	write(10,'(a)') ",spower      = 1.0          ! Planetesimal radius power-law exponent"
	write(10,'(a,1pe8.2,a)') ",sigma0      = ", sigma_s, "     ! Solid surface density at 1 AU (g/cm^2), fiducial: 7.0 (Hayashi MMSN)"
	write(10,'(a)') ",asnow       = 2.7          ! Snow line distance (AU)"
	write(10,'(a)') ",fsnow       = 4.2          ! Solid material enhancement beyond the snow line"
	write(10,'(a)') ",rhoplsml    = 0.5          ! Planetesimal mass density (g/cm^3)"
	write(10,'(a)') ",mratio      = 2.5E19       ! Ratio of mpl(1) to mplsmls (assuming mpl(1) is the most massive)"
	write(10,'(a,/)') "/"

	! Planetesimal eccentricity and inclination parameters
	write(10,'(a)') "&eiplsml_params"
	write(10,'(a)') " e0          = -1.0         ! Planetesimal constant eccentricity"
	write(10,'(a)') ",i0          = -1.0         ! Planetesimal constant inclination (radians)"
	write(10,'(a)') ",eipower     = -1.0         ! Planetesimal eccentricity and inclination power-law exponent"
	write(10,'(a)') ",emin        = 0.0          ! Planetesimal minimum eccentricity"
	write(10,'(a)') ",emax        = 1.0          ! Planetesimal maximum eccentricity"
	write(10,'(a)') ",imin        = 0.0          ! Planetesimal minimum inclination (radians)"
	write(10,'(a)') ",imax        = 1.570796     ! Planetesimal maximum inclination (radians)"
	write(10,'(a)') ",rmse        = 0.02         ! Planetesimal RMS eccentricity"
	write(10,'(a)') ",rmsi        = 0.01         ! Planetesimal RMS inclination (radians)"
	write(10,'(a,/)') "/"

	! Gas parameters
	write(10,'(a)') "&gas_params"
	write(10,'(a,1pe8.2,a)') " deng0s     = ", rho_g, "      ! Gas volume density at 1 AU (g/cm^3), fiducial: 1.4e-9 (Hayashi MMSN)"
	write(10,'(a)') ",gpower     = 2.75          ! Gas density power-law exponent, fiducial: 2.75 (Hayashi MMSN)"
	write(10,'(a,1pe8.2,a)') ",rcomet     = ", r, "      ! Comet radius used for aerodynamic drag (km)"
	write(10,'(a)') ",dcomet     = 0.5           ! Comet density used for aerodynamic drag (g/cm^3)"
	write(10,'(a)') ",zscale     = 0.05          ! Gas disk scale height (AU) at 1 AU, fiducial: 0.047 (Hayashi MMSN)"
	write(10,'(a,1pe8.2,a)') ",rgi        = ", amin, "      ! Gas disk inner edge distance (AU)"
	write(10,'(a,1pe8.2,a)') ",rgf        = ", amax, "      ! Gas disk outer edge distance (AU)"
	write(10,'(a)') ",tgdecay    = 1.0e29        ! Gas disk decay time (years)"
	write(10,'(a)') ",ca         = 1.0e-29       ! Semimajor axis Type-I drag efficiency, typically: 1.0"
	write(10,'(a)') ",ce         = 1.0e-29       ! Eccentricity Type-I drag efficiency, typically: 1.0"
	write(10,'(a)') ",ngap       = 0             ! Number of gaps in the gas disk"
	write(10,'(a)') ",igap       = 2             ! Particle ID(s) that are situated in a gap"
	write(10,'(a)') ",fwgap      = 0.05          ! Fractional width of the gap"
	write(10,'(a)') "/"

	close(unit = 10)

	end program tools_param_gen