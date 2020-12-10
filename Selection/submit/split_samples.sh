#!/bin/bash

sample=${1}
outFile=${2}
nSplits=${3}
removeOriginal=${4:-false}

nLines=$(( ( $(cat ${sample} | wc -l) / ${nSplits} ) + 1 ))

#echo $nLines

split --numeric-suffixes=1 --additional-suffix=.txt  -l $nLines $sample $outFile

echo "Total lines orginal $(cat $sample | wc -l)"

for s in $(ls ${outFile}*)
do
  echo "${s} has $(cat ${s} | wc -l) lines"
done

if $removeOriginal
then
  rm -v $sample
fi

