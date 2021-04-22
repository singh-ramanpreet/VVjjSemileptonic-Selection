#!/bin/bash

base_ntuples_location=${1}

for year in 2016 2017 2018
do
  for era in A B C D E F G H
    do
    ntuples=${base_ntuples_location}/${year}
    files=$(eos root://cmseos.fnal.gov find -name *Run${year}${era}*.root ${ntuples} | sed -e 's|^|root://cmseos.fnal.gov/|' | tr '\n' ' ')
    condor_submit \
      universe=vanilla \
      executable=job.sh \
      transfer_input=True \
      transfer_output=True \
      transfer_input_files="RemoveDuplicateEvents.C,run.sh,job.sh" \
      log_filename="../condor_logs/postjob_${year}_${era}_\$(Cluster)_\$(Process)" \
      log="\$(log_filename).log" \
      output="\$(log_filename).out" \
      error="\$(log_filename).err" \
      -append "arguments=${era} ${ntuples} ${files}" \
      -append "queue 1" \
      /dev/null
  done
done
exit
