#!/bin/bash

year=${1}
ntuples=${2}

# Logs dir
mkdir -p condor_logs

# make tar
make_tar.sh
cmssw_tar=$(basename ${CMSSW_BASE}).tgz

eras=$(eos root://cmseos.fnal.gov find -name "*Run${year}*.root" ${ntuples} | sed -E 's/.*Run'${year}'([A-Z]).*\.root/\1/' | sort -u | tr '\n' ' ')
echo $eras

for era in ${eras}
  do
  files=$(eos root://cmseos.fnal.gov find -name *Run${year}${era}*.root ${ntuples} | sed -e 's|^|root://cmseos.fnal.gov/|' | tr '\n' ' ')
  echo $files
  condor_submit \
    universe=vanilla \
    executable=$(which run_remove_duplicate.sh) \
    transfer_input=True \
    transfer_output=True \
    log="/dev/null" \
    output="A" \
    error="B" \
    transfer_input_files="${cmssw_tar}" \
    transfer_output_files="\"\"" \
    -append "arguments=${cmssw_tar} ${era} ${ntuples} ${files}" \
    -append "queue 1" \
    /dev/null
done
exit
