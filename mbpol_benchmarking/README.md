# Benchmarking the function obtaining the potential energy in 4d polynomial 

## General
This sub-branch is to benchmark the run-time cost of the function that calculates the potential energy of a water dimer in 4d polynomial.
There are four inserted timers by default trying to time these functions: 
   - function `E_poly = poly_2b_v6x::eval()` in `mbpol\x2b-v9x.cpp:x2b_v9x::eval`
   - function `E2poly += x2b_v9x::eval()` in `mbpol\mbpol.cpp:mbpol::operator()`
   - function `E = pot(nw, crd, grd);` in `c++\test-mbpol_openmp.cpp:main()`
   - function `E_nogrd = pot(nw, crd);` in `c++\test-mbpol_openmp.cpp:main()`  
     
The code is designed valid for OPENMP parallel programming

A tester (test-mbpol_openmp.cpp) is provided with OPENMP-enabled.  
The tester reads in a total of 42508 sample dimers from file `2b.xyz`, sets up the above timers, and run the 4d_poly program with all dimers.   
Coming along with the tester, a bash script is used to configue threads for OPENMP and to make the run time statistics.


## Timer class
The timer class is in files `c++\timestamp.h/cpp`  
To use a timer, follow these steps:  
   - `timers_t timers` to declare a timer class instance
   - `timerid_t id` to declare an `id` as type `timerid_t`. `id` represents the unique ID for each timer.
   - Initialize a timer:
      - `timers.insert_random_timer(id, threadid=0, label="")` to initialize a timer. The `id` will be returned via reference. 
      - or, one may `timers.insert_timer(id, threadid=0, label="")` to initialize a timer with a specific id(=0~999)
      - `int threadid[=0]` is a user-defined labeling integer
      - `string label[=""]` is a user-defined labeling string 
   - Start a timer: `timers.start_timer(id)`
   - End a timer: `timers.end_timer(id, ifadd=true, ifsave=false)`
      - `bool ifadd[=true]` sets if the metered time will be cumulatively added. Only timers with the same `threadid` and `label` will be added cumulatively.
      - `bool ifsave[=false]` sets if the timer will be deleted after the end call.
   - Show the cumulative time: `timers.get_time_collections()` will print out the cumulative time, `threadid`, and `label`.
   - Show all saved timers: `timers.get_all_timers_info()` 


## Folders and files list:
- `mbpol`                           : Create mbpol library
- `c++`                             : Run mbpol benchmarking
   - `test-mbpol_openmp.cpp`        : The tester file
   - `run_bench_openmp.x`           : The script file 
   - `timestamps.h/cpp`             : The timer class files
   
## To Compile and run:
   - In folder `mbpol`, run `make` to create the library
        - The default setting in Makefile uses Intel Compiler `icpc` and flags `-xHost -fopenmp`
   - In folder `c++`, unzip the input data file `2b.xyz.bz2` by `bunzip2 2b.xyz.bz2`
   - Run `make` to create tester
   - Run `run_bench_openmp.x` script to benchmark
        - The default benchmark is to run with [1 2 4 8 16 24] OPENMP threads.
        - Default repeated runs for each OPENMP thread is set to 10
        - Output statistics shows an average/min/max run time of each function on each thread.

## Runtime result: 
The following test is run on Skylate with intel compiler `icpc -fopenmp -xHost` configuration

 Statistics of function runtime in [ms] for one thread :  
 
| OMP_THRDS |   E_poly  MEAN  |      E_2poly   MEAN  |     E_nogrd    MEAN    | 
| --------- | ------- | ------- | ------- |
|     1     | 2.2e+05  |2.7e+05 | 8.2e+05 |  
|     2     | 1.2e+05  | 1.5e+05  | 4.3e+05 | 
|     4     | 5.9e+04  |  7.4e+04 | 2.2e+05  |  
|     8     | 3.3e+04  | 4.2e+04 | 1.3e+05 |   
|    16     | 1.8e+04  |2.3e+04 |  6.9e+04 | 
|    24     | 1.4e+04  | 1.8e+04 | 5.2e+04 | 

