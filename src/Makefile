# Makefile for HFMM beam-beam package -- LINUX G77

.SUFFIXES: .f .c

F77 = gfortran -Dunix -Dg77 -fno-second-underscore -fno-automatic -w -c

CC = gcc -c

OBJS = fastion.o beambeamnew.o sin2.o
LIBS = -lgfortran -Lfftw -lfftw

### BASIC
FOPTS = -O2 -g -march=native -DNDEBUG
### PROFILE
#FOPTS = -O2 -g -march=native -DNDEBUG -fno-omit-frame-pointer \
        -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls
### PEDANTIC
#FOPTS = -O2 -g -march=native -DNDEBUG -Wall -Wshadow -Wcast-align

### BASIC
COPTS = -O2 -g -march=native -DNDEBUG \
        -I/afs/cern.ch/sw/lcg/external/GSL/1.14/x86_64-slc5-gcc44-opt/include \
        -I/afs/cern.ch/sw/lcg/external/GSL/1.14/x86_64-slc5-gcc44-opt/include/gsl
### PROFILE
#COPTS = -O2 -g -march=native -DNDEBUG -fno-omit-frame-pointer \
        -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls\
        -I/afs/cern.ch/sw/lcg/external/GSL/1.14/x86_64-slc5-gcc44-opt/include \
        -I/afs/cern.ch/sw/lcg/external/GSL/1.14/x86_64-slc5-gcc44-opt/include/gsl
### PEDANTIC
#COPTS = -O2 -g -march=native -DNDEBUG \
        -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
        -I/afs/cern.ch/sw/lcg/external/GSL/1.14/x86_64-slc5-gcc44-opt/include \
        -I/afs/cern.ch/sw/lcg/external/GSL/1.14/x86_64-slc5-gcc44-opt/include/gsl

### BASIC
LOPTS = -g -O2 -march=native -DNDEBUG
### PROFILE
#LOPTS = -g -O2 -march=native -DNDEBUG -lprofiler -fno-omit-frame-pointer -fno-inline-functions \
        -fno-inline-functions-called-once -fno-optimize-sibling-calls
### PEDANTIC
#LOPTS = -g -O2 -march=native -DNDEBUG -Wall -Wshadow -Wpointer-arith \
        -Wcast-qual -Wcast-align

# -Wall -> warning all

.F.o:
	$(F77) $(FOPTS) $<

.f.o:
	$(F77) $(FOPTS) $<

.c.o: Makefile
	$(CC) $(COPTS) $<

fastion: $(OBJS)
	g++ -lm $(LOPTS) -o ../fastion $(OBJS) $(LIBS)
