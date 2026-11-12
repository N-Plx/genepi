CCSRCS     = $(wildcard *.cc)
# The executable is linked directly in the top directory, next to gpd_table.dat
EXE        = ../genepi

CXX        = g++

# ROOT: use $ROOTSYS if set, otherwise root-config from the PATH
ifdef ROOTSYS
ROOTCONFIG = $(ROOTSYS)/bin/root-config
else
ROOTCONFIG = root-config
endif

# Local headers first: ROOT ships its own hepevt.h, which must not shadow ../include/hepevt.h
CCFLAGS    = -I../include -I. -g -Wall -fPIC -Wno-deprecated
CCFLAGS   += $(shell $(ROOTCONFIG) --cflags)

ROOTLIB    = $(shell $(ROOTCONFIG) --glibs)

# gfortran runtime: find its directory (needed on macOS / non-system compilers)
GFORTRAN_LIB := $(shell gfortran -print-file-name=libgfortran$(if $(filter Darwin,$(shell uname -s)),.dylib,.so))
GFORTRAN_DIR := $(dir $(GFORTRAN_LIB))

JLIB       = ../jetset_lib/libjetset.a
# ludata.o (JETSET BLOCK DATA) is linked explicitly: some linkers do not pull a
# BLOCK DATA member out of an archive, leaving the JETSET defaults at zero.
MYLIBS     = ../jetset_lib/ludata.o -L../jetset_lib -ljetset -L$(GFORTRAN_DIR) -Wl,-rpath,$(GFORTRAN_DIR) -lgfortran -lm
LIBS       = $(ROOTLIB) $(MYLIBS)
CCOBJS     = $(CCSRCS:%.cc=%.o)

all: $(JLIB) $(EXE)

$(JLIB):
	@echo compiling Jetset lib
	$(MAKE) -C ../jetset_lib

$(EXE): $(CCOBJS) $(JLIB)
	$(CXX) $(CCFLAGS) -o $@ $(CCOBJS) $(LIBS)

.cc.o :
	@echo Compiling $@
	@$(CXX) $(CCFLAGS) -c $< -o $@

clean:
	rm -f $(EXE) $(CCOBJS) genepi

cleanall: clean
	rm -f ../jetset_lib/*.o $(JLIB)
