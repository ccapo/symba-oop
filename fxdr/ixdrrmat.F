c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer*4 function ixdrrmat( ixdrs, nels, rval )

	implicit none

	integer*4   ixdrs, ixdrsm1, ierr, nels
	real	  rval(nels)

	ixdrsm1 = ixdrs - 1

	call cxdrrmat( ixdrsm1, nels, rval, ierr )

	ixdrrmat = ierr

	return
	end


