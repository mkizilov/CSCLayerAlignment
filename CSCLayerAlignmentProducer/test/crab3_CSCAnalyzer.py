from CRABClient.UserUtilities import config
config = config()

#section general
config.General.requestName = '2023-12-20_NOREFIT_v1'
config.General.workArea = 'crab3_CSCAnalyzer' #working dir
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'CSCAnalyzer_cfg.py'
config.JobType.numCores = 1


misalign = True
if misalign:
  config.JobType.inputFiles = ['./Run2023D_prompt_CSC_zeroGPR_zeroCSC_v1_03.db', './GlobalAlignment_Run2_Run3_v1_ZeroMuonGPR.db']

#section Data
config.Data.inputDataset = '/Muon1/Run2023D-MuAlCalIsolatedMu-PromptReco-v2/ALCARECO'



config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'

config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/mkizilov/crab3_out/CSCAnalyzer/2023-12-20_NOREFIT_v1/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Site.storageSite = 'T3_CH_CERNBOX'