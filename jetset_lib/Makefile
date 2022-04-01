fSRCS = $(wildcard *.f)
fOBJS = ${fSRCS:.f=.o}
LIB=libjetset.a

FF = gfortran
#FF = g77

#FFLAGS = -g -O -fno-automatic -fno-second-underscore \
# -ffixed-line-length-none -funroll-loops -Wunused -ffortran-bounds-check \
# -I../include.pepsi -I../include.lund

FFLAGS = -g -O -fno-automatic -fno-second-underscore \
 -ffixed-line-length-none -funroll-loops -Wunused

AR=ar
AR_FLAGS=rv
RANLIB=ranlib

$(LIB): $(fOBJS)
	@echo "" ; echo "Creating library $(LIB) ..."
	@$(AR) ${ARFLAGS} $(LIB) $(FOBJS) $(fOBJS) > /dev/null
	@$(RANLIB) $(LIB)

.f.o :
	${FF} -c ${FFLAGS} $<


clean:
	/bin/rm -f $(LIB) *.o
