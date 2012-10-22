c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer*4 function ixdrimat( ixdrs, nels, ival )

	implicit none

	integer*4   ixdrs, ixdrsm1, ierr, nels
	integer*4   ival(nels)

	ixdrsm1 = ixdrs - 1
	call cxdrimat( ixdrsm1, nels, ival, ierr )

	ixdrimat = ierr

	return
	end


