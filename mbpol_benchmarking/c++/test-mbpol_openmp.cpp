#include <cmath>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <omp.h>
#include <vector>
#include <time.h>

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
    if (argc != 2) {
        std::cerr << "usage: test-mbpol water.xyz"
                  << std::endl;
        return 0;
    }

    std::cout << std::scientific << std::setprecision(9);
     
    std::vector< std::string > elements_all;
    std::vector< double      > crd_all;
    size_t natoms=6;

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

    std::cout<<"Total number of molecular=" << elements_all.size()/natoms << std::endl;
    
     
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
for(int ii=0; ii<(elements_all.size()/natoms); ii++)
{        
         timerid_t timerid;
         int threadid=0;
#ifdef _OPENMP        
         threadid = omp_get_thread_num();
#endif           
         timers_t & timers_this_thread = alltimers[threadid];   // Rename `alltimers[threadid]` to `timers_this_thread`
         
                  
         std::vector< std::string > elements;         
         std::vector< double      > crd;

         //for(int jj=0; jj<6; jj++){      
         //      elements.push_back(elements_all[ii*6+jj]);                                           
         //      crd.push_back(crd_all[ii*18+jj*3  ]) ;
         //      crd.push_back(crd_all[ii*18+jj*3+1]) ;
         //      crd.push_back(crd_all[ii*18+jj*3+2]) ;               
         //}

	for(int jj=0; jj<natoms; jj++){
		elements.push_back(elements_all[ii*natoms+jj]);
		crd.push_back(crd_all[ii*natoms*3+jj*3  ]) ;
		crd.push_back(crd_all[ii*natoms*3+jj*3+1]) ;
		crd.push_back(crd_all[ii*natoms*3+jj*3+2]) ;
	}

         const size_t nw = elements.size()/3;

         x2o::mbpol pot;     
          
         double grd[9*nw]; // where does 9 come from?
         
         // Insert a timer and record its threadid and label. Return by reference the unique id.
         // The id is needed with the start and end timer functions are called.
         timers_this_thread.insert_random_timer(timerid,threadid,"E(PIP)w/grd");   
         timers_this_thread.timer_start(timerid); 
         double E = pot(nw, &(crd[0]), grd);
         timers_this_thread.timer_end(timerid);
         

         timers_this_thread.insert_random_timer(timerid,threadid,"E(PIP)w/ogrd");
         timers_this_thread.timer_start(timerid);
         double E_nogrd = pot(nw, &(crd[0]));
         timers_this_thread.timer_end(timerid);
         
         /*
         // A part for testing the runtime of timer functions
         timerid_t randins, ins, startid, stopid;
         timers_this_thread.insert_random_timer(randins,threadid,"Test_Random_Insert");
         timers_this_thread.timer_start(randins);
         timers_this_thread.insert_random_timer(timerid,threadid);
         timers_this_thread.timer_end(randins);
         
         timers_this_thread.insert_random_timer(ins,threadid,"Test_Fix_Insert");
         timers_this_thread.timer_start(ins);
         timers_this_thread.insert_timer(1, threadid);
         timers_this_thread.timer_end(ins);
         
         timers_this_thread.insert_random_timer(startid,threadid,"Test_Start");
         timers_this_thread.timer_start(startid);
         timers_this_thread.timer_start(timerid);
         timers_this_thread.timer_end(startid);
         
         timers_this_thread.insert_random_timer(stopid,threadid,"Test_End");
         timers_this_thread.timer_start(stopid);
         timers_this_thread.timer_end(timerid);
         timers_this_thread.timer_end(stopid);
         */
         
         
         
       
/*
#ifdef _OPENMP                            
        #pragma omp critical 
        {
#endif
         std::cout << "       E = " << E << "\n"
                   << "E[nogrd] = " << E_nogrd << std::endl;
#ifdef _OPENMP               
        }
#endif
*/

// Tester for another part in the code
    /*      
         const double eps = 1.0e-4;
         for (int n = 0; n < 9*nw; ++n) {        
             
             double tmp[9*nw];
             const double x_orig = crd[n];

             crd[n] = x_orig + eps;
             const double Ep = pot(nw, &(crd[0]), tmp);

             crd[n] = x_orig - eps;
             const double Em = pot(nw, &(crd[0]), tmp);

             const double gfd = 0.5*(Ep - Em)/eps;
             crd[n] = x_orig;        
          
     //        std::cout << grd[n] << ' ' << gfd << '\n';
             //std::cout << n << ' ' << std::fabs(grd[n] - gfd) << '\n';
            
         }    
     */

     //#pragma omp critical
     //{
          //alltimers.get_all_timers_info();
     //}   
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




    return 0;
}
