c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer*4 function ixdrreal( ixdrs, rval )

	implicit none

	integer*4   ixdrs, ixdrsm1, ierr
	real	  rval

	ixdrsm1 = ixdrs - 1

	call cxdrreal( ixdrsm1, rval, ierr )

	ixdrreal = ierr

	return
	end


