#!/bin/bash

input_list_dir=${1}
year=${2}
output_eos_dir=${3}
resubmit=${4:-"no"}
nano_ver=${5:-"7"}

# Logs dir
mkdir -p condor_logs

# make tar
make_tar.sh
cmssw_tar=$(basename ${CMSSW_BASE}).tgz

for dataset_txt in $(ls ${input_list_dir}/*.txt)
do

  if [[ "$resubmit" == "yes" ]]
  then
    check_file_exists=$(eos root://cmseos.fnal.gov ls ${output_eos_dir}/${year}/$(basename ${dataset_txt%.txt}).root &> /dev/null && echo true || echo false)
  else
    check_file_exists=false 
  fi
  $check_file_exists && continue

  condor_submit \
    universe=vanilla \
    executable="$(which run_selection.sh)" \
    transfer_input=True \
    transfer_output=True \
    stream_error=True \
    stream_output=True \
    log_filename="condor_logs/\$Fn(${year}_$(basename ${dataset_txt}))" \
    log="/dev/null" \
    output="\$(log_filename).out" \
    error="\$(log_filename).err" \
    transfer_input_files="${dataset_txt},${cmssw_tar}" \
    transfer_output_files="\"\"" \
    -append "arguments = ${cmssw_tar} $(basename ${dataset_txt}) ${year} ${nano_ver} ${output_eos_dir}/${year}/" \
    -append "queue" \
    /dev/null
done

exit
