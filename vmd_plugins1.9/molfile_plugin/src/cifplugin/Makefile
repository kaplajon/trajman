# Makefile for molecule file readers

.SILENT:

.SUFFIXES:

PLUGINAPI = molfile_plugin.h vmdplugin.h
COMPILEDIR = ../compile
ARCHDIR=${COMPILEDIR}/lib_${ARCH}/molfile

SRCDIR = src
INCDIR = -I../include -I${SRCDIR}

VPATH = src ../include ${ARCHDIR}

SCCFLAGS = $(CCFLAGS) $(DEF)"STATIC_PLUGIN"
SCXXFLAGS = $(CCFLAGS) $(DEF)"STATIC_PLUGIN"

#
# Rules
#

STATICPLUGINS = biomoccaplugin cifplugin cpmdplugin psfplugin pdbplugin dcdplugin babelplugin gromacsplugin parmplugin crdplugin namdbinplugin binposplugin graspplugin msmsplugin stlplugin cubeplugin edmplugin ccp4plugin dsn6plugin brixplugin pltplugin raster3dplugin parm7plugin rst7plugin tinkerplugin uhbdplugin vaspchgcarplugin vaspposcarplugin vaspoutcarplugin vaspxdatcarplugin vaspxmlplugin dlpolyplugin lammpsplugin vtfplugin xyzplugin corplugin moldenplugin phiplugin pbeqplugin grdplugin situsplugin dxplugin spiderplugin mapplugin avsplugin fs4plugin pqrplugin mol2plugin gridplugin carplugin mdfplugin gamessplugin xsfplugin bgfplugin xbgfplugin
PLUGINS = biomoccaplugin.so cifplugin.so cpmdplugin.so psfplugin.so pdbplugin.so dcdplugin.so babelplugin.so gromacsplugin.so parmplugin.so crdplugin.so namdbinplugin.so binposplugin.so graspplugin.so msmsplugin.so stlplugin.so cubeplugin.so edmplugin.so ccp4plugin.so dsn6plugin.so brixplugin.so pltplugin.so raster3dplugin.so parm7plugin.so rst7plugin.so tinkerplugin.so uhbdplugin.so vaspchgcarplugin.so vaspposcarplugin.so vaspoutcarplugin.so vaspxdatcarplugin.so vaspxmlplugin.so dlpolyplugin.so lammpsplugin.so vtfplugin.so xyzplugin.so corplugin.so moldenplugin.so phiplugin.so pbeqplugin.so grdplugin.so situsplugin.so dxplugin.so spiderplugin.so mapplugin.so avsplugin.so fs4plugin.so pqrplugin.so mol2plugin.so gridplugin.so carplugin.so mdfplugin.so gamessplugin.so xsfplugin.so bgfplugin.so xbgfplugin.so

#
# Check to see if we should build the Tcl-based plugins
# XXX non-portable GNU make syntax used here...
#
ifdef TCLLIB
ifdef TCLINC
ifdef TCLLDFLAGS
STATICPLUGINS += webpdbplugin
PLUGINS       += webpdbplugin.so
endif
endif
endif

#
# Check to see if we should build the optional NetCDF-based plugins
# XXX non-portable GNU make syntax used here...
#
ifdef NETCDFLIB
ifdef NETCDFINC
ifdef NETCDFLDFLAGS
STATICPLUGINS += netcdfplugin
PLUGINS       += netcdfplugin.so
endif
endif
endif

# list of all optional plugins for use by distrib target
OPTPLUGINS = webpdbplugin.so netcdfplugin.so

STATICS = libmolfile_plugin.a libmolfile_plugin.h
WIN32STATICS = libmolfile_plugin.lib libmolfile_plugin.h
DISTFILES = $(PLUGINS) $(OPTPLUGINS) $(STATICS) $(WIN32STATICS)

bins:
win32bins:
dynlibs: ${ARCHDIR} mesg $(PLUGINS)
staticlibs: ${ARCHDIR} $(STATICS) 
win32staticlibs: ${ARCHDIR} $(WIN32STATICS)

distrib:
	@echo "Copying molfile plugins to $(PLUGINDIR) destination area"
	for file in $(DISTFILES) ; do \
		echo "    $$file ..."; \
		for localname in `find ../compile -name $$file -print`; do \
			pluginname=`echo $$localname | sed s/..\\\/compile\\\/lib_// `; \
			dir=`dirname $(PLUGINDIR)/$$pluginname`; \
			mkdir -p $$dir; \
			cp -p $$localname $(PLUGINDIR)/$$pluginname; \
		done; \
	done;

	
mesg:
	@echo "Building Molecule File Reader plugins"


#
# plugin library rules
#

biomoccaplugin.so: ${ARCHDIR}/biomoccaplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

cifplugin.so: ${ARCHDIR}/cifplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

cpmdplugin.so: ${ARCHDIR}/cpmdplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

psfplugin.so: ${ARCHDIR}/psfplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

pdbplugin.so: ${ARCHDIR}/pdbplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

dcdplugin.so: ${ARCHDIR}/dcdplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

babelplugin.so: ${ARCHDIR}/babelplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

gromacsplugin.so: ${ARCHDIR}/gromacsplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

parmplugin.so: ${ARCHDIR}/parmplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

netcdfplugin.so: ${ARCHDIR}/netcdfplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(NETCDFLIB) $(NETCDFLDFLAGS) $(LDFLAGS)

crdplugin.so: ${ARCHDIR}/crdplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

namdbinplugin.so: ${ARCHDIR}/namdbinplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

binposplugin.so: ${ARCHDIR}/binposplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

webpdbplugin.so: ${ARCHDIR}/webpdbplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(TCLLIB) $(TCLLDFLAGS) $(LDFLAGS)

graspplugin.so: ${ARCHDIR}/graspplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

msmsplugin.so: ${ARCHDIR}/msmsplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

stlplugin.so: ${ARCHDIR}/stlplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

cubeplugin.so: ${ARCHDIR}/cubeplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

edmplugin.so: ${ARCHDIR}/edmplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

ccp4plugin.so: ${ARCHDIR}/ccp4plugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

dsn6plugin.so: ${ARCHDIR}/dsn6plugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

brixplugin.so: ${ARCHDIR}/brixplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

pltplugin.so: ${ARCHDIR}/pltplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

raster3dplugin.so: ${ARCHDIR}/raster3dplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

rst7plugin.so: ${ARCHDIR}/rst7plugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

parm7plugin.so: ${ARCHDIR}/parm7plugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

vtfplugin.so: ${ARCHDIR}/vtfplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

xyzplugin.so: ${ARCHDIR}/xyzplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

lammpsplugin.so: ${ARCHDIR}/lammpsplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

tinkerplugin.so: ${ARCHDIR}/tinkerplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

uhbdplugin.so: ${ARCHDIR}/uhbdplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

vaspchgcarplugin.so: ${ARCHDIR}/vaspchgcarplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

vaspposcarplugin.so: ${ARCHDIR}/vaspposcarplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

vaspoutcarplugin.so: ${ARCHDIR}/vaspoutcarplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

vaspxdatcarplugin.so: ${ARCHDIR}/vaspxdatcarplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

vaspxmlplugin.so: ${ARCHDIR}/vaspxmlplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

dlpolyplugin.so: ${ARCHDIR}/dlpolyplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

corplugin.so: ${ARCHDIR}/corplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

moldenplugin.so: ${ARCHDIR}/moldenplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

phiplugin.so: ${ARCHDIR}/phiplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

pbeqplugin.so: ${ARCHDIR}/pbeqplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

grdplugin.so: ${ARCHDIR}/grdplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

dxplugin.so: ${ARCHDIR}/dxplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

situsplugin.so: ${ARCHDIR}/situsplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

spiderplugin.so: ${ARCHDIR}/spiderplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

mapplugin.so: ${ARCHDIR}/mapplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

avsplugin.so: ${ARCHDIR}/avsplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

fs4plugin.so: ${ARCHDIR}/fs4plugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

pqrplugin.so: ${ARCHDIR}/pqrplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

mol2plugin.so: ${ARCHDIR}/mol2plugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

bgfplugin.so: ${ARCHDIR}/bgfplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

xbgfplugin.so: ${ARCHDIR}/xbgfplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

gridplugin.so: ${ARCHDIR}/gridplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

carplugin.so: ${ARCHDIR}/carplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

mdfplugin.so: ${ARCHDIR}/mdfplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

gamessplugin.so: ${ARCHDIR}/gamessplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

xsfplugin.so: ${ARCHDIR}/xsfplugin.o 
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)

#
# object files
#
${ARCHDIR}/biomoccaplugin.o: biomoccaplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/cifplugin.o: cifplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/cpmdplugin.o: cpmdplugin.c ${PLUGINAPI} 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/psfplugin.o: psfplugin.c fortread.h ${PLUGINAPI} 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/pdbplugin.o: pdbplugin.c readpdb.h ${PLUGINAPI}
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/dcdplugin.o: dcdplugin.c ${PLUGINAPI}
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/babelplugin.o: babelplugin.c readpdb.h vmddir.h ${PLUGINAPI}
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/gromacsplugin.o: gromacsplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/parmplugin.o: parmplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/netcdfplugin.o: netcdfplugin.c ${PLUGINAPI}
	$(CC) $(CCFLAGS) $(NETCDFINC) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/crdplugin.o: crdplugin.c ${PLUGINAPI}
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/namdbinplugin.o: namdbinplugin.c ${PLUGINAPI}
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/binposplugin.o: binposplugin.c ${PLUGINAPI}
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/webpdbplugin.o: webpdbplugin.c readpdb.h ${PLUGINAPI} 
	$(CC) $(CCFLAGS) $(TCLINC) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/graspplugin.o: graspplugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/msmsplugin.o: msmsplugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/stlplugin.o: stlplugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/cubeplugin.o: cubeplugin.C ${PLUGINAPI} periodic_table.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/edmplugin.o: edmplugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/ccp4plugin.o: ccp4plugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/dsn6plugin.o: dsn6plugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/brixplugin.o: brixplugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/pltplugin.o: pltplugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/raster3dplugin.o: raster3dplugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/rst7plugin.o: rst7plugin.c ${PLUGINAPI} 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/parm7plugin.o: parm7plugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/vtfplugin.o: vtfplugin.c ${PLUGINAPI}
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/xyzplugin.o: xyzplugin.c ${PLUGINAPI} periodic_table.h
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/lammpsplugin.o: lammpsplugin.c ${PLUGINAPI} periodic_table.h
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/tinkerplugin.o: tinkerplugin.c ${PLUGINAPI}
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/uhbdplugin.o: uhbdplugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/vaspchgcarplugin.o: vaspchgcarplugin.c vaspplugin.h ${PLUGINAPI} periodic_table.h 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/vaspposcarplugin.o: vaspposcarplugin.c vaspplugin.h ${PLUGINAPI} periodic_table.h 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/vaspoutcarplugin.o: vaspoutcarplugin.c vaspplugin.h ${PLUGINAPI} periodic_table.h 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/vaspxdatcarplugin.o: vaspxdatcarplugin.c vaspplugin.h ${PLUGINAPI} periodic_table.h 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/vaspxmlplugin.o: vaspxmlplugin.c vaspplugin.h ${PLUGINAPI} periodic_table.h 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/dlpolyplugin.o: dlpolyplugin.c ${PLUGINAPI} periodic_table.h
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/corplugin.o: corplugin.c ${PLUGINAPI} 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/moldenplugin.o: moldenplugin.c ${PLUGINAPI} 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/phiplugin.o: phiplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/pbeqplugin.o: pbeqplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/grdplugin.o: grdplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/dxplugin.o: dxplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/situsplugin.o: situsplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/spiderplugin.o: spiderplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/mapplugin.o: mapplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/avsplugin.o: avsplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/fs4plugin.o: fs4plugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/pqrplugin.o: pqrplugin.c ${PLUGINAPI} 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/mol2plugin.o: mol2plugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/xbgfplugin.o: xbgfplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/bgfplugin.o: bgfplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/gridplugin.o: gridplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/carplugin.o: carplugin.c ${PLUGINAPI} 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/mdfplugin.o: mdfplugin.C ${PLUGINAPI} 
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/gamessplugin.o: gamessplugin.c ${PLUGINAPI} 
	$(CC) $(CCFLAGS) $(INCDIR) -c $< $(COPTO)$@

${ARCHDIR}/xsfplugin.o: xsfplugin.C ${PLUGINAPI} periodic_table.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $< $(COPTO)$@

#
# archive rules
#

ARCHIVEOBJS = ${ARCHDIR}/biomoccaplugin-s.o ${ARCHDIR}/cifplugin-s.o ${ARCHDIR}/cpmdplugin-s.o ${ARCHDIR}/psfplugin-s.o ${ARCHDIR}/pdbplugin-s.o ${ARCHDIR}/dcdplugin-s.o ${ARCHDIR}/babelplugin-s.o ${ARCHDIR}/gromacsplugin-s.o ${ARCHDIR}/parmplugin-s.o ${ARCHDIR}/crdplugin-s.o ${ARCHDIR}/namdbinplugin-s.o ${ARCHDIR}/binposplugin-s.o ${ARCHDIR}/graspplugin-s.o ${ARCHDIR}/msmsplugin-s.o ${ARCHDIR}/stlplugin-s.o ${ARCHDIR}/edmplugin-s.o ${ARCHDIR}/ccp4plugin-s.o ${ARCHDIR}/dsn6plugin-s.o ${ARCHDIR}/brixplugin-s.o ${ARCHDIR}/pltplugin-s.o ${ARCHDIR}/cubeplugin-s.o ${ARCHDIR}/raster3dplugin-s.o ${ARCHDIR}/rst7plugin-s.o ${ARCHDIR}/parm7plugin-s.o ${ARCHDIR}/vtfplugin-s.o ${ARCHDIR}/xyzplugin-s.o ${ARCHDIR}/tinkerplugin-s.o ${ARCHDIR}/uhbdplugin-s.o ${ARCHDIR}/vaspchgcarplugin-s.o ${ARCHDIR}/vaspposcarplugin-s.o ${ARCHDIR}/vaspoutcarplugin-s.o ${ARCHDIR}/vaspxdatcarplugin-s.o ${ARCHDIR}/vaspxmlplugin-s.o ${ARCHDIR}/dlpolyplugin-s.o ${ARCHDIR}/lammpsplugin-s.o ${ARCHDIR}/corplugin-s.o ${ARCHDIR}/moldenplugin-s.o ${ARCHDIR}/phiplugin-s.o ${ARCHDIR}/pbeqplugin-s.o ${ARCHDIR}/grdplugin-s.o ${ARCHDIR}/situsplugin-s.o ${ARCHDIR}/dxplugin-s.o ${ARCHDIR}/spiderplugin-s.o ${ARCHDIR}/mapplugin-s.o ${ARCHDIR}/avsplugin-s.o ${ARCHDIR}/fs4plugin-s.o ${ARCHDIR}/pqrplugin-s.o ${ARCHDIR}/mol2plugin-s.o ${ARCHDIR}/gridplugin-s.o ${ARCHDIR}/carplugin-s.o ${ARCHDIR}/mdfplugin-s.o ${ARCHDIR}/gamessplugin-s.o ${ARCHDIR}/xsfplugin-s.o ${ARCHDIR}/bgfplugin-s.o ${ARCHDIR}/xbgfplugin-s.o


#
# Check to see if we should build the Tcl-based plugins
# XXX non-portable GNU make syntax used here...
#
ifdef TCLLIB
ifdef TCLINC
ifdef TCLLDFLAGS
ARCHIVEOBJS += ${ARCHDIR}/webpdbplugin-s.o
endif
endif
endif


#
# Check to see if we should build the optional NetCDF-based plugins
# XXX non-portable GNU make syntax used here...
#
ifdef NETCDFLIB
ifdef NETCDFINC
ifdef NETCDFLDFLAGS
ARCHIVEOBJS += ${ARCHDIR}/netcdfplugin-s.o
endif
endif
endif

libmolfile_plugin.a: ${ARCHIVEOBJS}
	rm -f ${ARCHDIR}/$@ 
	$(AR) cr ${ARCHDIR}/$@ ${ARCHIVEOBJS} 
	$(RANLIB) ${ARCHDIR}/$@

libmolfile_plugin.lib: ${ARCHIVEOBJS}
	rm -f ${ARCHDIR}/$@ 
	lib /OUT:${ARCHDIR}/$@ ${ARCHIVEOBJS}

libmolfile_plugin.h: ${ARCHIVEOBJS}
	rm -f ${ARCHDIR}/$@
	touch ${ARCHDIR}/$@
	../create_static_header.sh MOLFILE molfile ${ARCHDIR}/$@ ${STATICPLUGINS}

#
# object files suitable for static linking
#
${ARCHDIR}/biomoccaplugin-s.o: biomoccaplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_biomoccaplugin" -c $< $(COPTO)$@

${ARCHDIR}/cifplugin-s.o: cifplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_cifplugin" -c $< $(COPTO)$@

${ARCHDIR}/cpmdplugin-s.o: cpmdplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_cpmdplugin" -c $< $(COPTO)$@

${ARCHDIR}/psfplugin-s.o: psfplugin.c fortread.h ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_psfplugin" -c $< $(COPTO)$@

${ARCHDIR}/pdbplugin-s.o: pdbplugin.c readpdb.h ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_pdbplugin" -c $< $(COPTO)$@

${ARCHDIR}/dcdplugin-s.o: dcdplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_dcdplugin" -c $< $(COPTO)$@

${ARCHDIR}/babelplugin-s.o: babelplugin.c readpdb.h vmddir.h ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_babelplugin" -c $< $(COPTO)$@

${ARCHDIR}/gromacsplugin-s.o: gromacsplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_gromacsplugin" -c $< $(COPTO)$@

${ARCHDIR}/parmplugin-s.o: parmplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_parmplugin" -c $< $(COPTO)$@

${ARCHDIR}/netcdfplugin-s.o: netcdfplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(NETCDFINC) $(INCDIR) $(DEF)"VMDPLUGIN=molfile_netcdfplugin" -c $< $(COPTO)$@

${ARCHDIR}/crdplugin-s.o: crdplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_crdplugin" -c $< $(COPTO)$@

${ARCHDIR}/namdbinplugin-s.o: namdbinplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_namdbinplugin" -c $< $(COPTO)$@

${ARCHDIR}/binposplugin-s.o: binposplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_binposplugin" -c $< $(COPTO)$@

${ARCHDIR}/webpdbplugin-s.o: webpdbplugin.c readpdb.h ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(TCLINC) $(DEF)"VMDPLUGIN=molfile_webpdbplugin" -c $< $(COPTO)$@

${ARCHDIR}/graspplugin-s.o: graspplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_graspplugin" -c $< $(COPTO)$@

${ARCHDIR}/msmsplugin-s.o: msmsplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_msmsplugin" -c $< $(COPTO)$@

${ARCHDIR}/stlplugin-s.o: stlplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_stlplugin" -c $< $(COPTO)$@

${ARCHDIR}/cubeplugin-s.o: cubeplugin.C ${PLUGINAPI} periodic_table.h
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_cubeplugin" -c $< $(COPTO)$@

${ARCHDIR}/edmplugin-s.o: edmplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_edmplugin" -c $< $(COPTO)$@

${ARCHDIR}/ccp4plugin-s.o: ccp4plugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_ccp4plugin" -c $< $(COPTO)$@

${ARCHDIR}/dsn6plugin-s.o: dsn6plugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_dsn6plugin" -c $< $(COPTO)$@

${ARCHDIR}/brixplugin-s.o: brixplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_brixplugin" -c $< $(COPTO)$@

${ARCHDIR}/pltplugin-s.o: pltplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_pltplugin" -c $< $(COPTO)$@

${ARCHDIR}/raster3dplugin-s.o: raster3dplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_raster3dplugin" -c $< $(COPTO)$@

${ARCHDIR}/rst7plugin-s.o: rst7plugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_rst7plugin" -c $< $(COPTO)$@

${ARCHDIR}/parm7plugin-s.o: parm7plugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_parm7plugin" -c $< $(COPTO)$@

${ARCHDIR}/vtfplugin-s.o: vtfplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_vtfplugin" -c $< $(COPTO)$@

${ARCHDIR}/xyzplugin-s.o: xyzplugin.c ${PLUGINAPI} periodic_table.h
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_xyzplugin" -c $< $(COPTO)$@

${ARCHDIR}/lammpsplugin-s.o: lammpsplugin.c ${PLUGINAPI} periodic_table.h
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_lammpsplugin" -c $< $(COPTO)$@

${ARCHDIR}/tinkerplugin-s.o: tinkerplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_tinkerplugin" -c $< $(COPTO)$@

${ARCHDIR}/uhbdplugin-s.o: uhbdplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_uhbdplugin" -c $< $(COPTO)$@

${ARCHDIR}/vaspchgcarplugin-s.o: vaspchgcarplugin.c vaspplugin.h ${PLUGINAPI} periodic_table.h 
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_vaspchgcarplugin" -c $< $(COPTO)$@

${ARCHDIR}/vaspposcarplugin-s.o: vaspposcarplugin.c vaspplugin.h ${PLUGINAPI} periodic_table.h 
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_vaspposcarplugin" -c $< $(COPTO)$@

${ARCHDIR}/vaspoutcarplugin-s.o: vaspoutcarplugin.c vaspplugin.h ${PLUGINAPI} periodic_table.h 
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_vaspoutcarplugin" -c $< $(COPTO)$@

${ARCHDIR}/vaspxdatcarplugin-s.o: vaspxdatcarplugin.c vaspplugin.h ${PLUGINAPI} periodic_table.h 
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_vaspxdatcarplugin" -c $< $(COPTO)$@

${ARCHDIR}/vaspxmlplugin-s.o: vaspxmlplugin.c vaspplugin.h ${PLUGINAPI} periodic_table.h 
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_vaspxmlplugin" -c $< $(COPTO)$@

${ARCHDIR}/dlpolyplugin-s.o: dlpolyplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_dlpolyplugin" -c $< $(COPTO)$@

${ARCHDIR}/corplugin-s.o: corplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_corplugin" -c $< $(COPTO)$@

${ARCHDIR}/moldenplugin-s.o: moldenplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_moldenplugin" -c $< $(COPTO)$@

${ARCHDIR}/phiplugin-s.o: phiplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_phiplugin" -c $< $(COPTO)$@

${ARCHDIR}/pbeqplugin-s.o: pbeqplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_pbeqplugin" -c $< $(COPTO)$@

${ARCHDIR}/grdplugin-s.o: grdplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_grdplugin" -c $< $(COPTO)$@

${ARCHDIR}/dxplugin-s.o: dxplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_dxplugin" -c $< $(COPTO)$@

${ARCHDIR}/situsplugin-s.o: situsplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_situsplugin" -c $< $(COPTO)$@

${ARCHDIR}/spiderplugin-s.o: spiderplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_spiderplugin" -c $< $(COPTO)$@

${ARCHDIR}/mapplugin-s.o: mapplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_mapplugin" -c $< $(COPTO)$@

${ARCHDIR}/avsplugin-s.o: avsplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_avsplugin" -c $< $(COPTO)$@

${ARCHDIR}/fs4plugin-s.o: fs4plugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_fs4plugin" -c $< $(COPTO)$@

${ARCHDIR}/pqrplugin-s.o: pqrplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_pqrplugin" -c $< $(COPTO)$@

${ARCHDIR}/mol2plugin-s.o: mol2plugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_mol2plugin" -c $< $(COPTO)$@

${ARCHDIR}/xbgfplugin-s.o: xbgfplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_xbgfplugin" -c $< $(COPTO)$@

${ARCHDIR}/bgfplugin-s.o: bgfplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_bgfplugin" -c $< $(COPTO)$@

${ARCHDIR}/gridplugin-s.o: gridplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_gridplugin" -c $< $(COPTO)$@

${ARCHDIR}/carplugin-s.o: carplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_carplugin" -c $< $(COPTO)$@

${ARCHDIR}/mdfplugin-s.o: mdfplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_mdfplugin" -c $< $(COPTO)$@

${ARCHDIR}/gamessplugin-s.o: gamessplugin.c ${PLUGINAPI}
	${CC} ${SCCFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_gamessplugin" -c $< $(COPTO)$@

${ARCHDIR}/xsfplugin-s.o: xsfplugin.C ${PLUGINAPI} periodic_table.h
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_xsfplugin" -c $< $(COPTO)$@

${ARCHDIR} :
	mkdir -p ${ARCHDIR}

clean:
	find ${COMPILEDIR} \( -name *.o -o -name *.a -o -name *.so -o -name *.dll \) -print | xargs rm -f

