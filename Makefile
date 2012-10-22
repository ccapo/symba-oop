#******************************************************************************
#
#  Unit Name   : Makefile
#  Unit Type   : makefile
#  Project     : SWIFT
#  Package     : N/A
#  Language    : GNU makefile syntax
#
#  Description : Controls, via the make program, the building of the swift
#                modules, library, drivers, and tools, as well as initiating
#                the build of the FXDR library by means of its own makefile
#
#  Input
#    Arguments : Zero or more of the following targets:
#                (1) all     : builds modules, entire swift library, FXDR
#                              library, swift drivers and tools
#                (2) mod     : builds modules
#                (3) lib     : builds entire swift library
#                (4) libdir  : compiles local directory source and adds the
#                              resulting objects to the swift library
#                (5) fxdr    : builds FXDR library by invoking its makefile
#                (6) drivers : builds swift drivers
#                (7) tools   : builds swift tools
#                (8) bin     : compiles local directory source and installs
#                              resulting executables to $(SWIFT_HOME)/bin
#                (9) clean   : removes all soft links to Makefile and
#                              Makefile.Defines from subdirectories of
#                              $(SWIFT_HOME), removes the entire contents
#                              of $(SWIFT_HOME)/lib and $(SWIFT_HOME)/bin,
#                              and removes the include file installed by the
#                              FXDR makefile
#    Terminal  : none
#    File      : Makefile.Defines
#
#  Output
#    Arguments : none
#    Terminal  : status messages
#    File      : none
#
#  Invocation  : make [all|mod|lib|libdir|fxdr|fxdr_install|drivers|tools|bin|clean]
#
#  Notes       : The use of the above arguments as phony targets inside the
#                makefile precludes their use as base names of swift drivers
#                or tools
#
#******************************************************************************

include Makefile.Defines

.PHONY : all mod lib libdir fxdr drivers tools bin clean force

% : %.f90 force
	$(FORTRAN) $(FFLAGS) -I$(SWIFT_HOME)/include $< -o $@ \
	  -L$(SWIFT_HOME)/lib -lswift -lfxdr
	$(INSTALL_PROGRAM) $@ $(SWIFT_HOME)/bin
	rm -f $@

all:
	cd $(SWIFT_HOME); \
	  make fxdr; \
	  make mod; \
	  make drivers; \
	  make tools

mod:
	cd $(SWIFT_HOME)/module; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFT_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFT_HOME)/Makefile .; \
	  $(FORTRAN) $(FFLAGS) -I$(SWIFT_HOME)/include -c *.f90; \
	  $(AR) r $(SWIFT_HOME)/lib/libswift.a *.o; \
	  $(RANLIB) $(SWIFT_HOME)/lib/libswift.a; \
	  $(INSTALL_DATA) *.mod $(SWIFT_HOME)/include; \
	  rm -f *.o *.mod

fxdr:
	cd $(SWIFT_HOME)/fxdr; \
	  rm -f Makefile.Defines; \
	  rm -f Makefile.fxdr; \
	  ln -s $(SWIFT_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFT_HOME)/Makefile.fxdr .; \
	  make -f Makefile.fxdr; \
	  make -f Makefile.fxdr test; \
	  make -f Makefile.fxdr install; \
	  make -f Makefile.fxdr clean

drivers:
	cd $(SWIFT_HOME)/main; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFT_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFT_HOME)/Makefile .; \
	  make bin

tools:
	cd $(SWIFT_HOME)/tools; \
	  rm -f Makefile.Defines Makefile; \
	  rm -f $(SWIFT_HOME)/bin/@plot_frame; \
	  ln -s $(SWIFT_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFT_HOME)/Makefile .; \
	  ln -s $(SWIFT_HOME)/tools/@plot_frame $(SWIFT_HOME)/bin; \
	  make bin

bin: *.f90
	make $(basename $^)

clean:
	cd $(SWIFT_HOME)/module;      rm -f Makefile.Defines Makefile
	cd $(SWIFT_HOME)/fxdr;        rm -f Makefile.Defines
	cd $(SWIFT_HOME)/main;        rm -f Makefile.Defines Makefile
	cd $(SWIFT_HOME)/tools;       rm -f Makefile.Defines Makefile
	cd $(SWIFT_HOME)/lib;         rm -f lib*.a
	cd $(SWIFT_HOME)/include;     rm -f *.mod
	cd $(SWIFT_HOME)/bin;         rm -f *

force:

#******************************************************************************
#
#  Author(s)   : David E. Kaufmann
#
#  Revision Control System (RCS) Information
#
#  Source File : $RCSfile: Makefile,v $
#  Full Path   : $Source: /d1/kaufmann/swift/RCS/Makefile,v $
#  Revision    : $Revision: 0.2 $
#  Date        : $Date: 2003/04/17 21:05:27 $
#  Programmer  : $Author: kaufmann $
#  Locked By   : $Locker:  $
#  State       : $State: Exp $
#
#  Modification History:
#
#  $Log: Makefile,v $
#  Revision 0.2  2003/04/17 21:05:27  kaufmann
#  modified Makefile to install modules in $(SWIFT_HOME)/include
#  rather than $(SWIFT_HOME)/lib, and added -I$(SWIFT_HOME)/include
#  switch to compile statements so that compiler can access modules
#  from that directory rather than having to be copied locally
#
#  Revision 0.1  2003/04/15 22:56:34  kaufmann
#  Initial implementation
#
#
#******************************************************************************
