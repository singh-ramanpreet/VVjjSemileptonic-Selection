VBS selection for semi-leptonic final states on nanoAOD

Initial setup:

```
  cmsrel CMSSW_10_2_13
  cd CMSSW_10_2_13/src
  git clone https://github.com/jaylawhorn/WVJJAna.git
  cd WVJJAna
  scram b
  . getfakerates.sh
```

Run interactively: 

```
  Selection test.dat output.root 1 2017 7
  Selection <list of input files> <output file> <1=MC, 0=data> <era(2016 or 2017)> <nanoaod version(7 only)>
```

Submit jobs to the fnal condor with submitjobs.sh


Notes: 

Scale factor file/histogram locations is defined in src/ScaleFactors.cc

submit folder has the locations of all nanoAOD v7 files

Full JEC uncertainties aren't saved

Event201*.h are MakeClass() dumps of a ttbar semileptonic file from each era for reference

Need to download fake rates from latino repository with getfakerates.sh


