c--------------------------------------------------------------
c	TOOLS_GEN_PL_SYMBA.F
c--------------------------------------------------------------
c Generates a SyMBA-style pl.in file with embryos and
c planetesimals.  The planetesimals are distributed in
c semimajor axis according to a power-law.  The mass/radius and
c the eccentricity/inclinations of the planetesimals can be
c given a constant value, distributed according to a power-law
c or a Rayleigh distribution (e and i only).
c--------------------------------------------------------------
	program tools_gen_pl_symba

	include 'swift.inc'

c... Parameters
	integer ialpha
	real*8 msun, mearth
	parameter(ialpha = -1)
	parameter(msun = twopi**2)
	parameter(mearth = 3.0d-6*msun)

c... Internal variables
	integer i, ierr, nrand, iseed, time(8)

c... Stellar variables
	logical*2 loblate, lclose
	real*8 j2, j4, rmin, rmax, rmaxu, qmin

c... Embryo variables
	logical*2 lgas
	integer npl
	real*8 rhopl, rhpl, rpl
	real*8 mpl(NPLMAX), apl(NPLMAX), epl(NPLMAX), wpl(NPLMAX), ipl(NPLMAX), opl(NPLMAX), mapl(NPLMAX)
	real*8 gm, x, y, z, vx, vy, vz

c... Planetesimal variables
	logical*2 lplsml, lshear
	character*24 ms_model, ei_model
	integer nplsml, ntotal, ndummy, util_plsml_number
	real*8 amin, amax, apower
	real*8 smin, smax, spower, sigma0, asnow, fsnow, rhoplsml, mratio
	real*8 e0, i0, eipower, emin, emax, imin, imax, rmse, rmsi
	real*8 mplsml, rhplsml, rplsml, mtotal, mtiny
	real*8 aplsml, eplsml, wplsml, iplsml, oplsml, maplsml
	real*8 util_randomu, util_power_law, a2, ia2, factor, ximin

c... Namelists
	namelist / stellar_params / loblate, j2, j4, lclose, rmin, rmax, rmaxu, qmin
	namelist / symba_params / mtiny
	namelist / pl_params / npl, rhopl, mpl, apl, epl, wpl, ipl, opl, mapl
	namelist / model_params / ms_model, ei_model, lplsml, lgas, lshear
	namelist / aplsml_params / amin, amax, apower
	namelist / msplsml_params / smin, smax, spower, sigma0, asnow, fsnow, rhoplsml, mratio ! Add oligarch option here
	namelist / eiplsml_params / e0, i0, eipower, emin, emax, imin, imax, rmse, rmsi

c-----
c...  Executable code

c...  Open and read in parameters
	open(unit = 7, file = "param.in", status = 'old')
	read(unit = 7, nml = stellar_params) ! Stellar parameters
	read(unit = 7, nml = symba_params)   ! SyMBA parameters
	read(unit = 7, nml = pl_params)      ! Embryo parameters
	read(unit = 7, nml = model_params)   ! Model parameters
	read(unit = 7, nml = aplsml_params)  ! Planetesimal semimajor axis parameters
	read(unit = 7, nml = msplsml_params) ! Planetesimal mass/radius parameters
	read(unit = 7, nml = eiplsml_params) ! Planetesimal eccentricity and inclination parameters
	close(unit = 7)

c... Program header
	write(*,'(/,a)') ' !-----------------------------------------------------------!'
	write(*,'(a)')   ' !                   TOOLS_GEN_PL_SYMBA                      !'
	write(*,'(a,/)') ' !-----------------------------------------------------------!'
	write(*,'(a,/)') ' Generates a SyMBA-style pl.in file with embryos and planetesimals'

c... Random seed based on current date and time
	call date_and_time(values = time)
	iseed = -13*sum(abs(time))
	call random_seed(size = nrand)
	call random_seed(put = (/ ( i*iseed, i = 1, nrand ) /))

c... Print the value of mtiny, then convert to sim. units
	write(*,'(a,1pe11.5,/)') ' The smallest mass to self gravitate (M_Earth) = ', mtiny
	mtiny = mtiny*mearth

c... Compute the number, mass and radius of the planetesimals
	if(lplsml) then

	  call util_ms_model(mpl(1), nplsml, mplsml, rplsml)

	  write(*,'(a,/)')              ' ------------------------'
	  write(*,'(a)')                ' Planetesimal parameters:'
	  write(*,'(a,i6)')             ' nplsml                = ', nplsml
	  write(*,'(a,1pe12.5,/)')      ' rhoplsml (g/cm^3)     = ', rhoplsml

	  write(*,'(a26,a24)')          ' Mass/Radius Model     =  ', ms_model

	  if((ms_model .eq. "constant") .or. (ms_model .eq. "CONSTANT")) then ! Constant mass and radius

	    write(*,'(a,1pe12.5)')      ' rplsml (AU)           = ', rplsml
	    write(*,'(a,1pe12.5,/)')    ' mplsml (M_Earth)      = ', mplsml

	  else if((ms_model .eq. "power") .or. (ms_model .eq. "POWER")) then  ! Power-law mass and radius distributions

	    write(*,'(a,2(1pe12.5))')   ' smin, smax (AU)       = ', smin, smax
	    write(*,'(a,1pe12.5)')      ' spower                = ', spower
	    rhoplsml = rhoplsml*5.60d11 ! Convert from [g/cm^3] to [M_Earth/AU^3]
	    smin = smin*6.685d-9        ! Convert from [km] to [AU]
	    smax = smax*6.685d-9        ! Convert from [km] to [AU]
	    write(*,'(a,1pe12.5)')      ' mplsml_min (M_Earth)  = ', (4.0d0/3.0d0)*pi*rhoplsml*smin**3
	    write(*,'(a,1pe12.5,/)')    ' mplsml_max (M_Earth)  = ', (4.0d0/3.0d0)*pi*rhoplsml*smax**3

	  endif

	  write(*,'(a,2(1pe12.5))')     ' amin, amax (AU)       = ', amin, amax
	  write(*,'(a,1pe12.5)')        ' power-law             = ', apower
	  write(*,'(a,1pe12.5)')        ' asnow (AU)            = ', asnow
	  write(*,'(a,1pe12.5,/)')      ' fsnow                 = ', fsnow

	  write(*,'(a26,a24)')          ' Ecc. and Inc. Model   =  ', ei_model
	  write(*,'(a,g2.0)')           ' Shear Regime          = ', lshear

c... Print eccentricity and inclinations in the shear regime
	  if(npl >= 1) then
	    ximin = 1.0d-2*minval(mpl(1:npl))**(1.0d0/3.0d0)
	  else
	    ximin = 1.0d-2*(mtiny/mearth)**(1.0d0/3.0d0)
	  endif
	  if(lshear) then
	    e0 = ximin/sqrt(twopi)
	    rmse = ximin/sqrt(twopi)
	  endif
	  if(i0 .gt. e0/2.0) i0 = e0/2.0
	  if(rmsi .gt. rmse/2.0) rmsi = rmse/2.0

	  if((ei_model .eq. "constant") .or. (ei_model .eq. "CONSTANT")) then ! Constant eccentricity and inclination

	    write(*,'(a,1pe12.5)')      ' e_0                   = ', e0
	    write(*,'(a,1pe12.5,/)')    ' i_0 (radians)         = ', i0

	  else if((ei_model .eq. "power") .or. (ei_model .eq. "POWER")) then  ! Power-law eccentricity and inclination distributions

	    write(*,'(a,2(1pe12.5))')   ' emin, emax            = ', emin, emax
	    write(*,'(a,2(1pe12.5))')   ' imin, imax (radians)  = ', imin, imax
	    write(*,'(a,1pe12.5,/)')    ' eipower               = ', eipower

	  else if((ei_model .eq. "rayleigh") .or. (ei_model .eq. "RAYLEIGH")) then ! Rayleigh eccentricity and inclination distributions

	    write(*,'(a,2(1pe12.5))')   ' emin, emax            = ', emin, emax
	    write(*,'(a,2(1pe12.5))')   ' imin, imax (radians)  = ', imin, imax
	    write(*,'(a,2(1pe12.5),/)') ' rmse, rmsi (radians)  = ', rmse, rmsi

	  endif

	else

	  nplsml = 0 ! No planetesimals

	endif

c... Now for the embryos
	if(lgas) then ! Adds an embryo for book-keeping purposes in the gas drag routines

	  npl       = npl + 1
	  mpl(npl)  = mtiny/mearth
c	  apl(npl)  = 10.0d0*amax
	  apl(npl)  = 0.999d0*rmax
	  epl(npl)  = 0.0d0
	  wpl(npl)  = 0.0d0
	  ipl(npl)  = 0.0d0
	  opl(npl)  = 0.0d0
	  mapl(npl) = 0.0d0

	endif

	write(*,'(a,/)')         ' ------------------------'
	write(*,'(a)')           ' Embryo parameters:      '
	write(*,'(a,i6)')        ' npl                   = ', npl
	write(*,'(a,1pe12.5)')   ' rhopl (g/cm^3)        = ', rhopl

c... Total number of particles
	ntotal = npl + nplsml + 1

c... Open the output file and write out variables, starting with the sun
	open(unit = 8, file = 'pl.in', status = 'unknown')

	x = 0.0d0
	y = 0.0d0
	z = 0.0d0
	vx = 0.0d0
	vy = 0.0d0
	vz = 0.0d0

	write(8,'(i6)') ntotal
	if(loblate) then ! Include the J2 and J4 stellar oblateness terms
	  write(8,'(3(1pe23.15))') msun, j2, j4
	else
	  write(8,'(1pe23.15)') msun
	endif
	write(8,'(3(1pe23.15))') x, y, z
	write(8,'(3(1pe23.15))') vx, vy, vz

c... Now add the embryos
	if(npl .ge. 1) then
	  write(*,'(a,51(1pe12.5))') ' mpl (M_Earth)         = ', ( mpl(i), i = 1, npl )
	  write(*,'(a,51(1pe12.5))') ' apl (AU)              = ', ( apl(i), i = 1, npl )
	  write(*,'(a,51(1pe12.5))') ' epl                   = ', ( epl(i), i = 1, npl )
	  write(*,'(a,51(1pe12.5))') ' ipl (degress)         = ', ( ipl(i), i = 1, npl )
	  write(*,'(a,51(1pe12.5))') ' capom (degress)       = ', ( opl(i), i = 1, npl )
	  write(*,'(a,51(1pe12.5))') ' omega (degress)       = ', ( wpl(i), i = 1, npl )
	  write(*,'(a,51(1pe12.5))') ' capm (degress)        = ', ( mapl(i), i = 1, npl )
	  rhopl = rhopl*5.60d11 ! Convert from [g/cm^3] to [M_Earth/AU^3]
	  do i = 1, npl
	    rpl = (3.0d0*mpl(i)/(4.0d0*pi*rhopl))**(1.0d0/3.0d0)
	    mpl(i) = mpl(i)*mearth
	    rhpl = apl(i)*(mpl(i)/(3.0d0*msun))**(1.0d0/3.0d0)
	    ipl(i) = ipl(i)/degrad
	    opl(i) = opl(i)/degrad
	    wpl(i) = wpl(i)/degrad
	    mapl(i) = mapl(i)/degrad
	    gm = msun + mpl(i)
	    call orbel_el2xv(gm, ialpha, apl(i), epl(i), ipl(i), opl(i), wpl(i), mapl(i), x, y, z, vx, vy, vz)
	    if(lclose) then ! Include planetary radius for close encounter detection
	      write(8,'(3(1pe23.15))') mpl(i), rhpl, rpl
	    else
	      write(8,'(2(1pe23.15))') mpl(i), rhpl
	    endif
	    write(8,'(3(1pe23.15))') x, y, z
	    write(8,'(3(1pe23.15))') vx, vy, vz
	  enddo
	endif
	write(*,'(/,a)') ' ------------------------'

c... Now add the planetesimals
	if(nplsml .ge. 1) then
	  do i = 1, nplsml
	    aplsml = util_power_law(apower, amin, amax)
	    call util_ms_model(mpl(1), ndummy, mplsml, rplsml) ! ndummy is used to prevent any changes to nplsml
	    mplsml = mplsml*mearth
	    rhplsml = aplsml*(mplsml/(3.0d0*msun))**(1.0d0/3.0d0)
	    call util_ei_model(eplsml, iplsml)
c	    eplsml = 1.0d0 - 30.0d0/aplsml ! eccentricity based on a q = 30 AU
	    oplsml = twopi*util_randomu()
	    wplsml = twopi*util_randomu()
	    maplsml = twopi*util_randomu()
	    gm = msun + mplsml
	    call orbel_el2xv(gm, ialpha, aplsml, eplsml, iplsml, oplsml, wplsml, maplsml, x, y, z, vx, vy, vz)
	    if(lclose) then ! Include planetesimal radius for close encounter detection
	      write(8,'(3(1pe23.15))') mplsml, rhplsml, rplsml
	    else
	      write(8,'(2(1pe23.15))') mplsml, rhplsml
	    endif
	    write(8,'(3(1pe23.15))') x, y, z
	    write(8,'(3(1pe23.15))') vx, vy, vz
	  enddo
	endif
	close(unit = 8)
	write(*,'(/,a)') '*** Done ***'

	end ! tools_gen_pl_symba.f