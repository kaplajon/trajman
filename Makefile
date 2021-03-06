PROG =	trajman

#---VMD MOLFILE PLUGIN---
VMD_ARCH = LINUXAMD64
#VMD_PLUGINS = /media/local/jon/src/vmd-1.8.7.src/plugins
VMD_PLUGINS=vmd_plugins1.9
#------------------------

SRCS =	module_input.F90 module_kinds.f90 module_readtraj.f90 \
	module_trajop.f90 module_util.f90 module_statistics.f90 trajman.f90 

OBJS =	module_input.o module_apl.o module_version.o module_kinds.o module_readtraj.o module_trajop.o \
	module_util.o module_statistics.o trajman.o 

LIBS =  #/usr/lib64/liblapack.so.3

#INCLUDES = -I/data/jon/src/LAPACK95/lapack95_modules/ 
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
NETCDFLDFLAGS  = -lnetcdf
TCLLDFLAGS     = -ltcl8.5
LDFLAGS        = $(VMD_PLUGINS)/molfile_plugin/f77/f77_molfile.o -L$(VMD_PLUGINS)/compile/lib_$(VMD_ARCH)/molfile
LDLIBS         = -lmolfile_plugin $(NETCDFLDFLAGS) $(TCLLDFLAGS) -lstdc++ -ldl

#LDFLAGS = 
ifeq ("$(shell echo $$(git log --max-count=1|grep commit))","$(shell if [ -f module_version.f90 ];then sed -ne "/^.commit.*/p" module_version.f90|cut -c 2-;else echo '0';fi)")
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
	echo $(shell echo ""\!"$$(git log --max-count=1|grep commit)\nmodule version\n    use kinds\n    character(kind=1,len=*),parameter :: branch=\"$$(git branch|awk '{print $2}')\",revision=\"$$(git log --max-count=1|grep commit)\",&\nrevdate=\"$$(git log --max-count=1|grep Date|cut -d ' ' -f 4-)\"\nend module version\n" >module_version.f90)
	git log > changelog.txt

#module_input.o: module_kinds.o module_readtraj.o module_util.o module_version.o
module_input.o: module_util.o module_version.o module_apl.o module_trajop.o
module_trajop.o: module_apl.o module_statistics.o
module_statistics.o: module_util.o module_readtraj.o
module_readtraj.o module_version.o module_trajop.o module_util.o module_statistics.o: module_kinds.o
#module_trajop.o: module_kinds.o
module_readtraj.o: module_util.o 
module_apl.o: module_readtraj.o
trajman.o: module_input.o module_kinds.o module_readtraj.o module_trajop.o \
	module_util.o module_statistics.o
