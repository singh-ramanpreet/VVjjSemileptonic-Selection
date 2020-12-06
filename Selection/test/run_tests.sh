#!/bin/bash

# Args
recompile=${1:-true}

if $recompile
then
  (cd $CMSSW_BASE/src/WVJJAna/Selection; scram b)
fi

for txtfile in $(ls *.txt)
do
  testfile=$(head -n 1 ${txtfile})
  isMC=$(python -c "import ROOT;f=ROOT.TFile.Open('root://cmseos.fnal.gov/$testfile'); t=f.Get('Runs');print(int(t.GetBranchStatus('genEventSumw')))")
  year=$(grep -Eo "201[6-8]" <<< "${txtfile}")
  txtpath=$(readlink -f ${txtfile})
  output=${txtpath%.txt}.root

  echo "Testing ${txtfile}, this is MC=$isMC, year=$year, with following command"
  echo Selection $txtpath $output $isMC $year 7
  (cd $CMSSW_BASE/src/WVJJAna/Selection; Selection $txtpath $output $isMC $year 7;)
done
