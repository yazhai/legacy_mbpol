#!/bin/bash

exefile="./test-openmp"
inpfile="2b.xyz"
logfile="runtime.log"
outfile="runtime_statistics.out"
omp_thread_list=(1 8)
itr_each_thread=1
interested_labels="E_poly   E_2poly  E_main  E_nogrd"  # interested timer sequence, separated by space



rm -f ${logfile}


## run bench with different threads

for ((j=0; j<${#omp_thread_list[@]} ; j++ ))
do
     omp_thread=${omp_thread_list[j]}
     export OMP_NUM_THREADS=${omp_thread}
     echo "OMP_NUM_THREADS =" $OMP_NUM_THREADS   >> $logfile
     for ((i=1; i<=$itr_each_thread; i++))
     do
          echo "Iteration No. = "  $i >> $logfile
          ${exefile} ${inpfile} >> ${logfile}
     done
done


## Readout logfile and make statstics
awk -v utfile=${outfile} -v allfds="${interested_labels}" '
BEGIN{     
##     PROCINFO["sorted_in"] = "@ind_str_asc";
     nfds=split(allfds,fd_mtx," ");     
     omp_threads=0;
     first_time=1;        
     
     printf ( " Statistics of function runtime in [ms] for one thread :\n           |") > utfile ;
    
     for (i=1; i<=nfds; i++){
          pos = int( (30 + length(fd_mtx[i]))/2 );
          printf ( "%*s%*s", pos , fd_mtx[i], 30-pos, "|") >> utfile;
     }
     printf ( "\n OMP_THRDS |" ) >> utfile;
     for (i=1; i<=nfds; i++){
          printf ( "    MEAN       MIN       MAX |")>> utfile;
     }
     printf ( "\n" )>> utfile;

     ##printf ( "   [ms]    |             E_poly          |        E_2poly                     \n"  );
     ##printf ( " OMP_THRDS |    MEAN       MIN       MAX |                                    \n"  );
     ##printf ( "    2      |     317  11111308       327 |                                    \n"  );     
}

$1~/OMP_NUM_THREADS/ {      
     if (first_time!=1) {     
          printf ( "%6d     |", omp_threads)  >> utfile;
          for ( i=1; i<=nfds; i++ ) {
               var = fd_mtx[i];
               printf ( "%8.2g  %8.2g  %8.2g |", time_mtx[var]/count[var]  , minT_mtx[var], maxT_mtx[var]  )>> utfile;
          }
          printf ( "\n" ) >> utfile;                    
     } else {
          first_time=0;
     }   
     delete time_mtx;
     delete maxT_mtx;
     delete minT_mtx;    
     delete count; 
     omp_threads=$3;  
}

##$1~/Iteration/ {  iter[omp_threads] = $4;}

(!($1~/OMP_NUM_THREADS/)) && (!($1~/Iteration/)){
     count[$4]++;
     time_mtx[$4]+= $6;
     if (maxT_mtx[$4] < $6) { maxT_mtx[$4] = $6; }
     if ( (minT_mtx[$4] > $6) || (minT_mtx[$4] == 0) ) { minT_mtx[$4] = $6; }  
}

END{
     printf ( "%6d     |", omp_threads)  >> utfile;
     for ( i=1; i<=nfds; i++ ) {
          var = fd_mtx[i];
          printf ( "%8.2g  %8.2g  %8.2g |", time_mtx[var]/count[var]  , minT_mtx[var], maxT_mtx[var]  )>> utfile;
     }
     printf ( "\n" ) >> utfile;
}' $logfile



cat ${outfile}



