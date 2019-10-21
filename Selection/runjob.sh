#!/bin/bash
echo "Starting job on " `date`
echo "Running on: `uname -a`"
echo "System software: `cat /etc/redhat-release`"
source /cvmfs/cms.cern.ch/cmsset_default.sh
#Arguments =  /eos/uscms/store/user/lpcbacon/15//SingleElectronRun2017B_31Mar2018_v1 /eos/uscms/store/user/klawhorn/WVJJTree_Oct21 SingleElectronRun2017B_31Mar2018_v1
### copy the input root files if they are needed only if you require local reading
xrdcp -s root://cmseos.fnal.gov//store/user/klawhorn/CMSSW_10_2_13.tgz  .
tar -xf CMSSW_10_2_13.tgz
rm CMSSW_10_2_13.tgz
cd CMSSW_10_2_13/src/WVJJAna/Selection/
echo "====> List files : " 
ls -alh
echo "====> Remove any file with name similar to WWTree*.root... " 
rm *.root
scramv1 b ProjectRename
eval `scram runtime -sh`
echo "====> List files : " 
ls -alh
/usr/bin/eos root://cmseos.fnal.gov ls ${1} | grep root > ${3}.txt
Selection ${3}.txt ${3}.root
echo "====> List files : " 
ls -alh
echo "====> List root files : " 
ls ${3}.root
echo "====> copying WWTree*.root file to stores area..." 
xrdcp -f ${3}.root root://cmseos.fnal.gov/${2}
rm ${3}.root
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_10_2_13
