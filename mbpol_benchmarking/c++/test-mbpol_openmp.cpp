#include <cmath>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdexcept>
#include <omp.h>
#include <vector>
#include <time.h>
#include <string.h>

#include "io-xyz.h"
#include "xyz-water-utils.h"
#include "mbpol.h"
#include "timestamps.h"

using namespace std;

// all timers, set to global instance for use in other subroutines.
vector<timers_t> alltimers;  



int main(int argc, char** argv)
{
     //srand(time(NULL));
    if (argc < 2) {
        std::cerr << "usage: test-mbpol XYZ_FILENAME [2|3 body model, default=2]"
                  << std::endl;
        return 0;
    }

    std::cout << std::scientific << std::setprecision(9);
     
    std::vector< std::string > elements_all;
    std::vector< double      > crd_all;
    size_t natoms=6;
    if (argc ==3){
          natoms = 3* atoi(argv[2]) ;
    }

    // Read out all records at once and save them to the vector containers.
    try {
        std::ifstream ifs(argv[1]);

        if (!ifs)
            throw std::runtime_error("could not open the XYZ file");
          
        std::string comment;
        
        while (!ifs.eof()) {
             kit::io::load_xyz(ifs, comment, elements_all, crd_all, natoms);
        }
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    
    size_t nsys = elements_all.size()/natoms;
    std::cout<<"Total number of molecular=" <<  nsys << std::endl;
    
    double* E_grd = new double[nsys];
    double* E_nogrd = new double[nsys]; 
    double* grd_all = new double[elements_all.size()*3]; 
     
    int num_threads=1; 
#ifdef _OPENMP
#pragma omp parallel shared(elements_all, crd_all, alltimers)
{
#pragma omp master
{
     num_threads=omp_get_num_threads();
#endif
     for (int ii=0; ii<num_threads; ii++){
          timers_t timers;
          alltimers.push_back(timers);    // Initialize timers according to the number of threads
     }
#ifdef _OPENMP     
}
#pragma omp barrier
#pragma omp for 
#endif
for(int ii=0; ii<(nsys); ii++)
{        
         timerid_t timerid;
         int threadid=0;
#ifdef _OPENMP        
         threadid = omp_get_thread_num();
#endif           
         timers_t & timers_this_thread = alltimers[threadid];   
         
                  
         std::vector< std::string > elements;         
         std::vector< double      > crd;


	    for(int jj=0; jj<natoms; jj++){
		     elements.push_back(elements_all[ii*natoms+jj]);
		     crd.push_back(crd_all[ii*natoms*3+jj*3  ]) ;
		     crd.push_back(crd_all[ii*natoms*3+jj*3+1]) ;
		     crd.push_back(crd_all[ii*natoms*3+jj*3+2]) ;
	    }

         const size_t nw = elements.size()/3;

         x2o::mbpol pot;     
          
         double grd[9*nw]; 
         
         // Insert a timer and record its threadid and label. Return by reference the unique id.
         // The id is needed with the start and end timer functions are called.
         timers_this_thread.insert_random_timer(timerid,threadid,"E(PIP)w/grd");   
         timers_this_thread.timer_start(timerid); 
         E_grd[ii] = pot(nw, &(crd[0]), grd);
         timers_this_thread.timer_end(timerid);
         memcpy(&(grd_all[ii*9*nw]), &(grd[0]), 9*nw*sizeof(double));

         timers_this_thread.insert_random_timer(timerid,threadid,"E(PIP)w/ogrd");
         timers_this_thread.timer_start(timerid);
         E_nogrd[ii] = pot(nw, &(crd[0]));
         timers_this_thread.timer_end(timerid);
         
}  


#ifdef _OPENMP
}
#endif

     // 

      for(int ii=0; ii<alltimers.size(); ii++){
           timers_t & timer = alltimers[ii];
           timer.get_time_collections();
           //timer.get_all_timers_info();
      }

      std::ifstream infile("E_grd.rst");
      if( !infile.good() ){
            std::ofstream ofs1 ("E_grd.rst", std::ofstream::out);
            std::ofstream ofs2 ("E_nogrd.rst", std::ofstream::out);
            std::ofstream ofs3 ("grd_all.rst", std::ofstream::out);
            ofs1 << std::scientific << std::setprecision(16);
            ofs2 << std::scientific << std::setprecision(16);
            ofs3 << std::scientific << std::setprecision(16);
            for(int ii=0; ii<nsys; ii++){
               ofs1 << E_grd[ii] << std::endl;
               ofs2 << E_nogrd[ii] << std::endl;
               for(int jj=0; jj<natoms*3; jj++){
                    ofs3 << grd_all[ii*natoms*3+jj] << "\t";
               }
               ofs3 << std::endl;
            } 
            ofs1.close();
            ofs2.close();
            ofs3.close();
      }

     delete [] E_grd;
     delete [] E_nogrd;
     delete [] grd_all;

    return 0;
}
