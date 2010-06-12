ARCH = LINUXAMD64

PROG =	trajman
VMD_PLUGINS = /home/jon/src/vmd-1.8.7.src/plugins
SRCS =	module_input.F90 module_kinds.f90 module_readtraj.f90 \
	module_trajop.f90 module_util.f90 trajman.f90

OBJS =	module_input.o module_apl.o module_version.o module_kinds.o module_readtraj.o module_trajop.o \
	module_util.o trajman.o 

LIBS =  #/usr/lib64/liblapack.so.3

#INCLUDES = -I/data/jon/src/LAPACK95/lapack95_modules/ 
#INCLUDES = -I/home/jon/src/vmd-1.8.7.src/plugins/compile/lib_$(ARCH)/molfile -I/home/jon/src/vmd-1.8.7.src/plugins/include
#-I/home/jon/src/LAPACK95/lapack95_modules/ 
#LIBSPATH = -L/usr/common/sprng2.0/lib
#LIBSPATH = /data/jon/src/LAPACK95/lapack95.a /usr/lib64/liblapack.so.3 
#LIBSPATH =# /home/jon/src/LAPACK95/lapack95.a /usr/lib/liblapack.so.3 
CC = gcc
CFLAGS = -03
F90 = gfortran 
#F90 = ifort 
FFLAGS =  -O3 
F90FLAGS =  -O3 -g 
FC90FLAGS = -D "CINFO='$$(date)'" $(F90FLAGS)
NETCDFLIB      = -L/home/jon/src/vmd-1.8.7.src/vmd/lib/netcdf/lib_LINUX
NETCDFLDFLAGS  = -lnetcdf
TCLLIB         = -L/home/jon/src/vmd-1.8.7.src/vmd/lib/tcl/lib_LINUX
TCLLDFLAGS     = -ltcl8.5
#F77MOLFILEPATH = /home/jon/src/vmd-1.8.7.src/plugins/molfile_plugin/f77/
LDFLAGS        = $(VMD_PLUGINS)/molfile_plugin/f77/f77_molfile.o -L/home/jon/src/vmd-1.8.7.src/plugins/compile/lib_$(ARCH)/molfile $(TCLLIB) $(NETCDFLIB)
LDLIBS         = -lmolfile_plugin $(NETCDFLDFLAGS) $(TCLLDFLAGS) -lstdc++ -ldl

#LDFLAGS = 
ifeq ("$(shell bzr version-info --custom --template="{revno}")","$(shell if [ -f module_version.f90 ];then sed -ne "/^.[0-9][0-9]*[0-9]*/p" module_version.f90|cut -c 2-5;else echo '0';fi)")
A := 
    else
A := $(shell if [ -f module_version.f90 ];then rm module_version.f90;fi)
B := $(shell echo "A new revision detected, creating new version module.")
C := $(shell if [ -f module_version.o ];then rm module_version.o;fi)
D := $(shell if [ -f version.mod ];then rm version.mod;fi)
endif

#all: version $(PROG)
all: $(PROG)


$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBSPATH) $(LDLIBS) $(LIBS)

#.PHONY: version
#version:

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90 

.f90.o:
	$(F90) $(F90FLAGS) -c $< $(LIBS) $(INCLUDES)

#.SUFFIXES: $(SUFFIXES) .f90
#
#.f.o:
#	$(F90) $(F90FLAGS) -c $< $(LIBS)

.SUFFIXES: $(SUFFIXES) .F90 

.F90.o:
	$(F90) $(FC90FLAGS) -c $<  $(LIBS) $(INCLUDES)

module_version.f90:
	@echo $(B)
	bzr version-info --custom --template="!{revno}\nmodule version\n    use kinds\n    character(kind=1,len=*),parameter :: branch=\"{branch_nick}\",revision=\"{revno}\",revdate=\"{date}\"\nend module version\n" >module_version.f90

#module_input.o: module_kinds.o module_readtraj.o module_util.o module_version.o
module_input.o: module_util.o module_version.o module_apl.o module_trajop.o
module_trajop.o: module_apl.o
module_readtraj.o module_version.o module_trajop.o module_util.o: module_kinds.o
#module_trajop.o: module_kinds.o
module_readtraj.o: module_util.o 
module_apl.o: module_readtraj.o
trajman.o: module_input.o module_kinds.o module_readtraj.o module_trajop.o \
	module_util.o
