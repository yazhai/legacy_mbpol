CXX = icpc
CXXFLAGS = -xHost -fopenmp -std=c++11
LIBFOLD = ../mbpol_icpc_xHost
CPPFLAGS = -I. -I$(LIBFOLD)
LDFLAGS = -xHost -fopenmp -std=c++11
GSL_FLAGS =  -lgsl -lgslcblas

all: test-mbpol 

test-mbpol: test-mbpol.o io-xyz.o xyz-water-utils.o
	$(CXX) -o test-mbpol test-mbpol.o io-xyz.o xyz-water-utils.o \
                  $(LIBFOLD)/libmbpol.a $(LDFLAGS)

clean:
	rm -f test-mbpol *.o

.PHONY: all clean
