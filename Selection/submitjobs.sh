#!/bin/bash

#baconFolder=/eos/uscms/store/user/lpcbacon/15/
baconFolder=/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/
#baconFolder=/eos/uscms/store/user/klawhorn/BaconSamples/

outputFolder=/eos/uscms/store/user/klawhorn/WVJJTree_Dec5/2016

#xrdfs root://cmseos.fnal.gov/ mkdir ${outputFolder}
#xrdfs root://cmseos.fnal.gov/ mkdir ${outputFolder}/log

outputJDL=dec_5_submit_2.jdl

cat stub.jdl > ${outputJDL}

while read line
do
    outputName=${line%%/*}

    echo "Output = " `pwd`/log/${outputName}.stdout >> ${outputJDL}
    echo "Error = " `pwd`/log/${outputName}.stderr >> ${outputJDL}
    echo "Log = " `pwd`/log/${outputName}.log >> ${outputJDL}
    echo "Arguments = " ${baconFolder}/${line} ${outputFolder} ${outputName} 1 2016 >> ${outputJDL}
    echo "Queue" >> ${outputJDL}
   
done < tosubmit2016_2.dat