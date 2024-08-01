# CSC Layer Alignment
## Check out package
```
scram p -n "CSCLayerAlignment" CMSSW CMSSW_14_0_3
cd CSCLayerAlignment/src
cmsenv
git cms-init
git clone https://github.com/mkizilov/CSCLayerAlignment
scram b -j 8
```

## Run Analyzer Part
To run analyzer part go to CSCLayerAlignmentProducer/test and run CSCAnalyzer_cfg.py or submit a corresponding crab job.

## Run Fitter
Hadd analyzer output and run fitter CSCLayerAlignmentProducer/scripts/run_CSCLayerAlignmentFitter.sh

## Apply shifts
After you run fitter you'll get .csv file with layer shifts. Use CSCLayerAlignmentProducer/plugins/CSCAlignmentDBWriter.cc to apply this shifts to geometry.

## Plot results
To plot results you can use https://github.com/mkizilov/MuonAlignmentAnalyser
