#!/bin/bash

logfilename="runtime.log"
benchfilename="benchmarkresult.txt"
time_at_col_in_log=12
testcount=100
totalcount=42105


rm -f ${logfilename} ${benchfilename}
for i in $(seq 1 ${testcount})
do
    awk 'NR > '\"$i\"'*8 -8 && NR <= '\"$i\"' *8 ' 2b.xyz >${i}.xyz
    ./test-mbpol ${i}.xyz >> ${logfilename}
    rm ${i}.xyz 
done

awk -v ut=${benchfilename} -v timecol=${time_at_col_in_log} -v total=${totalcount} '
BEGIN{
     nline=0;
     time_ms=0;
}
{
     nline=nline+1;
     time_ms=$timecol+time_ms;
}
END{
     printf ("Number of run samples (no gradient)   : %i \n", nline) > ut;
     printf ("Total   run time                 [ms] : %i \n", time_ms) >> ut;
     printf ("Average run time per sample      [ms] : %f \n", time_ms/nline) >> ut;
     printf ("Total run time of %6i samples [ms] : %f \n", total, time_ms/nline*total) >> ut;
}
' ${logfilename}

cat ${benchfilename}
