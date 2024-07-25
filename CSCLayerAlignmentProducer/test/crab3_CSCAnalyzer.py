from CRABClient.UserUtilities import config
config = config()

#section general
config.General.requestName = '2024-07-23v3'
config.General.workArea = 'crab3_CSCAnalyzer' #working dir
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'CSCAnalyzer_cfg.py'
config.JobType.numCores = 1


misalign = True
if misalign:
  config.JobType.inputFiles = ['./2024-07-22-CSC_2024B_json_03.db', './GlobalAlignment_Run2_Run3_v1_ZeroMuonGPR.db']

#section Data
config.Data.inputDataset = '/Muon1/Run2024B-MuAlCalIsolatedMu-PromptReco-v1/ALCARECO'
config.Data.lumiMask = './2024B.json'



config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'

config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/mkizilov/crab3_out/CSCAnalyzer/2024-07-23v3/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Site.storageSite = 'T3_CH_CERNBOX'