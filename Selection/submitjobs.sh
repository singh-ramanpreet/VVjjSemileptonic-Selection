#!/bin/bash

baconFolder=/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v6_3May
#baconFolder=/store/user/lpcbacon/arapyanRun2017/

outputFolder=/store/user/klawhorn/WVJJTree_Jun4/2018

outputJDL=may3_MC.jdl

cat stub.jdl > ${outputJDL}

while read line
do
    outputName=${line%%/*}
    testvar="$( eos root://cmseos.fnal.gov ls ${baconFolder}/${line} 2>&1 > /dev/null )"
    if [[ ${testvar} != "" ]]; then
        echo "Input folder does not exist... skipping" ${line}
        continue
    fi
    filelist=(`eos root://cmseos.fnal.gov ls ${baconFolder}/${line} | grep root`)
    #size=${#filelist[@]}
    echo ${#filelist[@]}
    i=0
    for inputFile in "${filelist[@]}"
    do
        echo "Output = " `pwd`/log/${outputName}_${i}.stdout >> ${outputJDL}
        echo "Error = " `pwd`/log/${outputName}_${i}.stderr >> ${outputJDL}
        echo "Log = " `pwd`/log/${outputName}_${i}.log >> ${outputJDL}
        echo "Arguments = " ${baconFolder}/${line}/${inputFile} ${outputFolder} ${outputName}_${i} 1 2018 6 >> ${outputJDL}
        echo "Queue" >> ${outputJDL}
        i=$((i+1))
    done
   
done < submit_2018_3May.txt