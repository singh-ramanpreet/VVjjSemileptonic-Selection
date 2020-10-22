#!/bin/bash

outputFolder=/store/user/klawhorn/WVJJTree_Oct22/

eosmkdir ${outputFolder}
eosmkdir ${outputFolder}/2016
eosmkdir ${outputFolder}/2017
eosmkdir ${outputFolder}/2018

cd ../../../..
rm CMSSW_10_6_10.tgz
tar -zcvf CMSSW_10_6_10.tgz CMSSW_10_6_10
xrdcp CMSSW_10_6_10.tgz root://cmseos.fnal.gov/${outputFolder}

#### change
outputJDL=2016_Background.jdl
####
cat stub.jdl > ${outputJDL}

while read line
do
    echo ${line[1]}
    outputName=${line[1]}
    echo ${outputName}
    testvar="$( eos root://cmseos.fnal.gov ls ${line[2]} 2>&1 > /dev/null )"
    if [[ ${testvar} != "" ]]; then
        echo "Input folder does not exist... skipping" ${line}
        continue
    fi
    filelist=(`eos root://cmseos.fnal.gov ls ${line[2]} | grep root`)
    #size=${#filelist[@]}
    echo ${#filelist[@]}
    i=0
    for inputFile in "${filelist[@]}"
    do
        echo "Output = " `pwd`/log/${outputName}_${i}.stdout >> ${outputJDL}
        echo "Error = " `pwd`/log/${outputName}_${i}.stderr >> ${outputJDL}
        echo "Log = " `pwd`/log/${outputName}_${i}.log >> ${outputJDL}
	#### inputFile outputFolder outputFile isMC[0/1] era[2016/2017/2018] version[6/7]
        echo "Arguments = " ${line[2]}/${inputFile} ${outputFolder} ${outputName}_${i} 1 2016 7 >> ${outputJDL}
	####
        echo "Queue" >> ${outputJDL}
        i=$((i+1))
    done
#### change
done < submit/run2016_Background_v7.txt 
####
