# Builds the test programs against the objects of ../src (run `make` at the top level first).
ifdef ROOTSYS
ROOTCONFIG = $(ROOTSYS)/bin/root-config
else
ROOTCONFIG = root-config
endif
CXXFLAGS = -I../include -O1 $(shell $(ROOTCONFIG) --cflags)
LIBS     = $(shell $(ROOTCONFIG) --libs)
OBJS     = ../src/dvcs.o ../src/nucl_ff.o ../src/read_optFile.o ../src/utils.o ../src/get_kinem_limits.o

all: xsec_points

xsec_points: xsec_points.cc $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ xsec_points.cc $(OBJS) $(LIBS)

clean:
	rm -f xsec_points
	rm -rf out
