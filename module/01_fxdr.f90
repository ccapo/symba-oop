module fxdr
! Definition of interfaces of functions in the FXDR (Fortran eXternal Data Representation) library
!
! FXDR is a library, written and maintained by David W. Pierce, that enables calls to the XDR
! (eXternal Data Representation) routines from Fortran.
!
! Reference : http://meteora.ucsd.edu/~pierce/fxdr_home_page.html
implicit none

interface
  function initxdr(filename, mode, returnonerror)
  character(*), intent(in) :: filename
  character(1), intent(in) :: mode
  logical, intent(in)      :: returnonerror
  integer                  :: initxdr
  end function initxdr
end interface

interface
  function ixdrclose(ixdr)
  integer, intent(in) :: ixdr
  integer             :: ixdrclose
  end function ixdrclose
end interface

interface
  function ixdrdmat(ixdrs, nels, dval)
  integer, intent(in)                           :: ixdrs, nels
  double precision, dimension(nels), intent(in) :: dval
  integer                                       :: ixdrdmat
  end function ixdrdmat
end interface

interface
  function ixdrdouble(ixdrs, dval)
  integer, intent(in)          :: ixdrs
  double precision, intent(in) :: dval
  integer                      :: ixdrdouble
  end function ixdrdouble
end interface

interface
  function ixdrimat(ixdrs, nels, ival)
  integer, intent(in)                  :: ixdrs, nels
  integer, dimension(nels), intent(in) :: ival
  integer                              :: ixdrimat
  end function ixdrimat
end interface

interface
  function ixdrint(ixdrs, ival)
  integer, intent(in) :: ixdrs, ival
  integer             :: ixdrint
  end function ixdrint
end interface

interface
  function ixdrreal(ixdrs, rval)
  integer, intent(in) :: ixdrs
  real, intent(in)    :: rval
  integer             :: ixdrreal
  end function ixdrreal
end interface

interface
  function ixdrrmat(ixdrs, nels, rval)
  integer, intent(in)               :: ixdrs, nels
  real, dimension(nels), intent(in) :: rval
  integer                           :: ixdrrmat
  end function ixdrrmat
end interface

end module fxdr
