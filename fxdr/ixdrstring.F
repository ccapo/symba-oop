c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer*4 function ixdrstring( ixdrs, string )

	implicit none

	integer*4   ixdrs, ixdrsm1, ierr, ilen
	character*(*) string

	ixdrsm1 = ixdrs - 1

	ilen = len(string)
	call cxdrstring( ixdrsm1, ilen, string, ierr )

	ixdrstring = ierr

	return
	end


