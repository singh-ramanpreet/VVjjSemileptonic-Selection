VBS selection for semi-leptonic final states on nanoAOD

Initial setup:

```
  cmsrel CMSSW_10_6_10
  cd CMSSW_10_6_10/src
  git clone https://github.com/singh-ramanpreet/VVjjSemileptonic-Selection.git VVjjSemileptonic
  cd VVjjSemileptonic
  scram b
  . getfakerates.sh
```

Run interactively: 

```
  Selection <list of input files> <output file> <1=MC, 0=data> <era(2016 or 2017)> <nanoaod version(7 only)>
```

Submit jobs to the fnal condor see these instructions [here](Selection/submit).
