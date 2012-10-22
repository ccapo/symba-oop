c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer*4 function ixdrclose( ixdr )

	implicit none

	integer*4 ixdr, ixdrm1, ierr

	ixdrm1 = ixdr - 1
	call cxdrclose( ixdrm1, ierr )

	ixdrclose = ierr

	return
	end
