VBS selection for semi-leptonic final states on 2016/2017 bacon files

* Initial setup:

  cmsrel CMSSW_10_2_13
  cd CMSSW_10_2_13/src
  git clone https://github.com/jaylawhorn/WVJJAna.git
  cd WVJJAna
  scram b

* Run interactively: 

  Selection test.dat output.root 1 2017
  Selection <list of input files> <output file> <1=MC, 0=data> <era(2016 or 2017)>


* Submit jobs to the fnal condor with submitjobs.sh