#!/bin/bash

echo "First doing the hadd of input list of root files"

rm Data.root

hadd Data.root $@

echo "Removing duplicate events now"

root -b -q -n -l RemoveDuplicateEvents.C\(\"Data.root\",\"noDuplicate.root\"\)

hadd Data_noDup.root noDuplicate.root

rootcp Data.root:/TotalEvents Data_noDup.root:/

rm noDuplicate.root

exit
