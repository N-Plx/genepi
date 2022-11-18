EXE        = genepi
CXX        = g++

CCSRCS     = $(wildcard *.cc)
FSRCS      = $(wildcard *.f)

CCFLAGS    = -g -Wall -fPIC -Wno-deprecated
CCFLAGS    += `$(ROOTSYS)/bin/root-config --cflags`
CCFLAGS    += -I../include -I. -I$(ROOTSYS)/include

ROOTLIB    = `$(ROOTSYS)/bin/root-config --glibs`
JLIB       = ../jetset_lib/libjetset.a
MYLIBS     = -Wl,-rpath, -lm -lgfortran -L../jetset_lib -ljetset
LIBS       = $(ROOTLIB) $(MYLIBS)

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
	rm -f $(EXE) $(CCOBJS) ../bin/$(EXE)

cleanall:
	rm -f $(EXE) $(CCOBJS) ../bin/$(EXE) ../jetset_lib/*.o $(JLIB)

