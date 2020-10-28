### How To

### Working conditions

1. All the commands should run from this directory.
2. By default for nanoAOD v7, for v6 make appropriate changes.

#### Make sample lists (current setup)

*Note*: Change path for v6

```bash
./make_sample_list.sh /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2016_v7_6July2020/ 2016
./make_sample_list.sh /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2016_v7_customNano_12Oct2020/ 2016
./make_sample_list.sh /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2017_v7_18Aug2020_Hadd/ 2017
./make_sample_list.sh /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2017_v7_customNano_12Oct2020/ 2017
./make_sample_list.sh /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v7_07Sep2020/ 2018
./make_sample_list.sh /eos/uscms/store/user/lnujj/VVjj_aQGC/nanoAOD_skim/Run2018_v7_customNano_12Oct2020/ 2018
```

### Make CMSSW sandbox

```bash
./make_tar.sh
```

### Submit condor jobs

*Note*: To run over v6, change the parameter `nano_ver` in `.jdl` file.

For 2017

```bash
condor_submit job_per_sample.jdl year=2017 sample_list_dir=2017 eos_output_path=/eos/uscms/store/user/rsingh/test/
```

This will loop over `*.txt` files in `sample_list_dir` and put output rootfiles in `eos_output_path/year` sub-directory.

For 2016, 2018

```bash
condor_submit job_per_sample.jdl year=2016 sample_list_dir=2016 eos_output_path=/eos/uscms/store/user/rsingh/test/
```

```bash
condor_submit job_per_sample.jdl year=2018 sample_list_dir=2018 eos_output_path=/eos/uscms/store/user/rsingh/test/
```
