CXX = icpc
CXXFLAGS = -std=c++11 -xHost
LIBFOLD = ../mbpol
CPPFLAGS = -I. -I$(LIBFOLD)
LDFLAGS = -fopenmp 

all: clean test-openmp

test-openmp: test-mbpol_openmp.o io-xyz.o xyz-water-utils.o timestamps.o
	$(CXX) -o test-openmp test-mbpol_openmp.o io-xyz.o xyz-water-utils.o timestamps.o \
                  $(LIBFOLD)/libmbpol.a $(CPPFLAGS) $(LDFLAGS) $(CXXFLAGS)

%.o: %.cpp 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -c -o $@ $<
	
clean:
	rm -f test-mbpol *.o

.PHONY: all clean
