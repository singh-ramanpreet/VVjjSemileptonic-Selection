#!/bin/bash

era=${1}
ntuples=${2}
files=${@:3}

source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc8-opt/setup.sh

[ -z "${files}" ] && exit

./run.sh ${era} ${files}

xrdcp Data_noDup_${era}.root root://cmseos.fnal.gov/${ntuples}/Data_noDup_${era}.root

rm -f docker_stderror

exit
