#ifndef GLOBAL_INCLUDE_H
#define GLOBAL_INCLUDE_H

#include <vector> 
#include <omp.h>
#include "timestamps.h"
/*
extern timers_t alltimers;
#ifdef _OPENMP
#pragma omp threadprivate(alltimers)
#endif
*/

extern std::vector<timers_t> alltimers; 

#endif
