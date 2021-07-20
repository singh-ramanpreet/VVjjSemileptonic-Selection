#!/bin/bash

a=$(eos root://cmseos.fnal.gov/ find -name "*.root" /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2016_v7_2021-07-15/)
readarray -t mc < <(grep -vE "Run201[6-8][A-Z]" <<<"$a");
declare -p mc &> /dev/null

readarray -t data < <(grep -E "Run201[6-8][A-Z]" <<<"$a");
declare -p data &> /dev/null

echo ${mc[0]} > test_2016_mc.txt
echo ${data[0]} > test_2016_data.txt


a=$(eos root://cmseos.fnal.gov/ find -name "*.root" /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2017_v7_2021-07-15/)
readarray -t mc < <(grep -vE "Run201[6-8][A-Z]" <<<"$a");
declare -p mc &> /dev/null

readarray -t data < <(grep -E "Run201[6-8][A-Z]" <<<"$a");
declare -p data &> /dev/null

echo ${mc[0]} > test_2017_mc.txt
echo ${data[0]} > test_2017_data.txt


a=$(eos root://cmseos.fnal.gov/ find -name "*.root" /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v7_2021-07-15/)
readarray -t mc < <(grep -vE "Run201[6-8][A-Z]" <<<"$a");
declare -p mc &> /dev/null

readarray -t data < <(grep -E "Run201[6-8][A-Z]" <<<"$a");
declare -p data &> /dev/null

echo ${mc[0]} > test_2018_mc.txt
echo ${data[0]} > test_2018_data.txt
