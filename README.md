# VBS selection for semi-leptonic final states

## Input:

Uses `nanoAOD` skims produced by https://github.com/singh-ramanpreet/VVjjSemileptonic-NanoSkim

## Setup:

```bash
  cmsrel CMSSW_10_6_22
  cd CMSSW_10_6_22/src
  git clone https://github.com/singh-ramanpreet/VVjjSemileptonic-Selection.git VVjjSemileptonic/Selection
  scram b
  cd VVjjSemileptonic/Selection
  get_fake_rates.sh
```

## Run interactively: 

```bash
  Selection <list of input files> <output file> <1=MC, 0=data> <era(2016 or 2017)> <nanoaod version(7 only)>
```

## Submit jobs to the FNAL HTCondor

### 1. Make sample lists

It will make samples list and split.

Working directory:
```bash
cd $CMSSW_BASE/src/VVjjSemileptonic/Selection
```

Code:

```bash
make_sample_list.sh <eos dir> <out list dir> <year> <split, n files per list> 
```

Quick use (current skims):

```bash
make_sample_list.sh /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2016_v7_2021-07-15/ inputs 2016 10
make_sample_list.sh /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2017_v7_2021-07-15/ inputs 2017 10
make_sample_list.sh /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v7_2021-07-15/ inputs 2018 10
```

### 2. Submit condor jobs

Working directory:

```bash
cd $CMSSW_BASE/src/VVjjSemileptonic/Selection/submit
```

Code: Submits over all `.txt` files in `dir`, output will be sub `dir` of `<eos dir>`.

```bash
submit_selection.sh <input list dir> <year> <eos dir> <resubmit> <nano version>
```

Quick use:

```bash
submit_selection.sh ../inputs/2016 2016 /eos/uscms/store/user/rsingh/test/
submit_selection.sh ../inputs/2017 2017 /eos/uscms/store/user/rsingh/test/
submit_selection.sh ../inputs/2018 2018 /eos/uscms/store/user/rsingh/test/
```

### 3. Combine `Data` and remove duplicate events

It combines `Data` per `era` (A, B, C, ... ).

Code:

```bash
submit_remove_duplicate.sh <year> <ntuples location>
```

Quick use:

```bash
submit_remove_duplicate.sh 2016 /eos/uscms/store/user/rsingh/test/
submit_remove_duplicate.sh 2017 /eos/uscms/store/user/rsingh/test/
submit_remove_duplicate.sh 2018 /eos/uscms/store/user/rsingh/test/
```
