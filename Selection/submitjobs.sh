#!/bin/bash

#baconFolder=/eos/uscms/store/user/lpcbacon/15/
#baconFolder=/eos/uscms/store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/
#baconFolder=/eos/uscms/store/user/klawhorn/BaconSamples/
baconFolder=/eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018/

outputFolder=/eos/uscms/store/user/klawhorn/WVJJTree_nanoAOD/

outputJDL=apr12_2.jdl

cat stub_dat.jdl > ${outputJDL}

while read line
do
    outputName=${line%%/*}

    echo "Output = " `pwd`/log/${outputName}.stdout >> ${outputJDL}
    echo "Error = " `pwd`/log/${outputName}.stderr >> ${outputJDL}
    echo "Log = " `pwd`/log/${outputName}.log >> ${outputJDL}
    echo "Arguments = " ${baconFolder}/${line} ${outputFolder} ${outputName} 0 2018 >> ${outputJDL}
    echo "Queue" >> ${outputJDL}
   
done < Datav5.dat