CCSRCS     = $(wildcard *.cc)
EXE        = genepi.exe
FSRCS      = $(wildcard *.f)

CXX        = g++
#CXX        = g++-3.4

CCFLAGS    = -g -Wall -fPIC -Wno-deprecated
CCFLAGS    += `$(ROOTSYS)/bin/root-config --cflags`
CCFLAGS    += -I../include -I. -I$(ROOTSYS)/include

ROOTLIB    = `$(ROOTSYS)/bin/root-config --glibs`
CERNLIB	   = `cernlib packlib mathlib kernlib`
#CERNLIB	   = -L$C(CERN_LIB) -static -lmathlib -lpacklib -lkernlib -Wl,-dy -llapack -lm -lnsl -lcrypt -ldl -lgfortran
JLIB       = ../jetset_lib/libjetset.a
#MYLIBS     = -Wl,-rpath, -lm -lg2c -L../jetset_lib -ljetset
MYLIBS     = -Wl,-rpath, -lm -lgfortran -L../jetset_lib -ljetset
LIBS       = $(ROOTLIB) $(CERNLIB) $(MYLIBS)

CCOBJS     = $(CCSRCS:%.cc=%.o)

all: $(JLIB) $(EXE)

$(JLIB):
	@echo compiling Jetset lib
	make -C ../jetset_lib

$(EXE): $(CCOBJS)
	$(CXX) $(CCFLAGS) -o $@ $(CCOBJS) $(FOBJS) $(LIBS)
	cp $(EXE) ../

.cc.o :
	@echo Compiling $@
	@$(CXX) $(CCFLAGS) -c $< -o $@

clean:
	rm -f $(EXE) $(CCOBJS) ../$(EXE)

cleanall:
	rm -f $(EXE) $(CCOBJS) ../$(EXE) ../jetset_lib/*.o $(JLIB)

