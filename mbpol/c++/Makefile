CXX = icc
CXXFLAGS = -O0 -fopenmp
CPPFLAGS = -I. -I../mbpol
LDFLAGS = -O0 -fopenmp
GSL_FLAGS =  -lgsl -lgslcblas

all: test-mbpol optimize normal-modes

test-mbpol: test-mbpol.o io-xyz.o xyz-water-utils.o
	$(CXX) -o test-mbpol test-mbpol.o io-xyz.o xyz-water-utils.o \
                  ../mbpol/libmbpol.a $(LDFLAGS)

optimize: optimize.o io-xyz.o xyz-water-utils.o io-utils.o xmin.o
	$(CXX) -o optimize optimize.o io-xyz.o xyz-water-utils.o \
	          io-utils.o xmin.o ../mbpol/libmbpol.a $(LDFLAGS)

normal-modes: normal-modes.o io-xyz.o xyz-water-utils.o io-utils.o xmin.o
	$(CXX) -o normal-modes normal-modes.o io-xyz.o xyz-water-utils.o \
	          io-utils.o xmin.o ../mbpol/libmbpol.a $(LDFLAGS) $(GSL_FLAGS)

clean:
	rm -f test-mbpol optimize normal-modes *.o

.PHONY: all clean
