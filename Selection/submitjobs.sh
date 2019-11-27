#!/bin/bash

baconFolder=/eos/uscms/store/user/lpcbacon/15/

outputFolder=/eos/uscms/store/user/klawhorn/WVJJTree_Nov25/2017

#xrdfs root://cmseos.fnal.gov/ mkdir ${outputFolder}
#xrdfs root://cmseos.fnal.gov/ mkdir ${outputFolder}/log

outputJDL=nov_25_submit.jdl

cat stub.jdl > ${outputJDL}

while read line
do
    outputName=${line%%/*}

    echo "Output = " `pwd`/log/${outputName}.stdout >> ${outputJDL}
    echo "Error = " `pwd`/log/${outputName}.stderr >> ${outputJDL}
    echo "Log = " `pwd`/log/${outputName}.log >> ${outputJDL}
    echo "Arguments = " ${baconFolder}/${line} ${outputFolder} ${outputName} 2017 >> ${outputJDL}
    echo "Queue" >> ${outputJDL}
   
done < tosubmit2017.dat