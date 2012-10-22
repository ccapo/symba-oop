c	-----------------------------
c	Returns < 0 if error, 0 if OK
c	-----------------------------
	integer*4 function ixdrint( ixdrs, ival )

	implicit none

	integer*4   ixdrs, ixdrsm1, ierr
	integer*4   ival

	ixdrsm1 = ixdrs - 1
	call cxdrint( ixdrsm1, ival, ierr )

	ixdrint = ierr

	return
	end


