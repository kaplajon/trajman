# hello emacs this is a -*- makefile -*-
# $Id: Makefile,v 1.1 2006/03/10 22:48:47 johns Exp $
########################################################################
ARCH           = LINUXAMD64
ARCHDIR        = ../../compile/lib_$(ARCH)/molfile
#FC            = gfortran
FC             = gfortran #-fno-second-underscore
CC             = gcc
LD             = $(FC)
OPT            = -O2 -Wall
CPPFLAGS       = -I../../include -I$(ARCHDIR)# -D_F77_F2C_UNDERSCORE
CFLAGS         = -c $(CPPFLAGS) $(OPT)
FFLAGS         = -c $(OPT)
NETCDFLIB      = -L../../../vmd/lib/netcdf/lib_LINUX
NETCDFLDFLAGS  = -lnetcdf
TCLLIB         = -L../../../vmd/lib/tcl/lib_LINUX
TCLLDFLAGS     = -ltcl8.5
LDFLAGS        = -L$(ARCHDIR) $(TCLLIB) $(NETCDFLIB)
LDLIBS         =  -lmolfile_plugin $(NETCDFLDFLAGS) $(TCLLDFLAGS) -lstdc++ -ldl
LDLIBS         =  -lmolfile_plugin -lstdc++ -ldl
########################################################################
#OBJECTS  = tester.o f77_molfile.o 
OBJECTS  = trr_test.o f77_molfile.o 

default: trr_test

trr_test: $(OBJECTS)
	$(FC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

clean:
	rm -f tester trr_test $(OBJECTS)
########################################################################
# pattern rules
.SUFFIXES:
.SUFFIXES: .c .f .f90 .o

.c.o:
	$(CC) $(CFLAGS) $< -o $@

.f90.o:
	$(FC) $(FFLAGS) $< -o $@

# dependencies
trr_test.o: trr_test.f90
#tester.o     : tester.f
f77_molfile.o: f77_molfile.c
