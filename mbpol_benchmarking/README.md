# Benchmarking the function obtaining the potential energy in 4d polynomial 

## General
This is a sub-branch benchmarking the run time cost of the function that calculates the potential energy of a water dimer in 4d polynomial.
The timer is set at the part obtaining potential energy in the function `mbpol/x2b-v9x::eval(const double* w1, const double* w2)`.
A total of 42508 samples are provided, in which only 42105 samples are available for benchmarking. The filted out samples have binding energy exceeding 60 kcal/mol. 

## Folder list
- `mbpol`                          : Create mbpol library
- `c++`                            : Run mbpol benchmarking

## To Compile and run:
   - In folder `mbpol`, run `make` to create the library
        - The default setting in Makefile uses Intel Compiler `icpc` and flags `-xHost -fopenmp`
   - In folder `c++`, unzip the input data file `2b.xyz.bz2` by `bunzip2 2b.xyz.bz2`
   - Run `make` to create tester
   - Run `run_benchmarking.x` script to benchmark
        - The default number of samples to benchmark is 100, in which 98 are valid.


