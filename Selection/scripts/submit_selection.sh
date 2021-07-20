#!/bin/bash

input_list_dir=${1}
year=${2}
output_eos_dir=${3}
nano_ver=${4:-"7"}

# Logs dir
mkdir -p condor_logs

# make tar
make_tar.sh
cmssw_tar=$(basename ${CMSSW_BASE}).tgz

for dataset_txt in $(ls ${input_list_dir}/*.txt)
do
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
    -append "arguments = ${cmssw_tar} ${dataset_txt} ${year} ${nano_ver} ${output_eos_dir}/${year}/" \
    -append "queue" \
    /dev/null
done

exit
