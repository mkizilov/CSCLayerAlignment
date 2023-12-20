from CRABClient.UserUtilities import config
config = config()

#section general
config.General.requestName = '2023-12-13_v2'
config.General.workArea = 'crab3_ME11_ana'#working dir
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_ME11ana.py'
# config.JobType.maxMemoryMB = 6000
# config.JobType.maxJobRuntimeMin = 2880 # 1440min = 24hours
config.JobType.numCores = 1
# config.JobType.allowUndistributedCMSSW = True
#config.JobType.generator
#config.JobType.pyCfgParams
#config.JobType.inputFiles


misalign = True
if misalign:
  config.JobType.inputFiles = ['./Run2023D_prompt_CSC_zeroGPR_zeroCSC_v1_03.db', './GlobalAlignment_Run2_Run3_v1_ZeroMuonGPR.db']

#section Data
config.Data.inputDataset = '/Muon1/Run2023D-MuAlCalIsolatedMu-PromptReco-v2/ALCARECO'
# config.Data.userInputFiles = open('Run2023D_Pv2.list').readlines()
# config.Data.runRange = '348776,348773,349073'



#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
# config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 1
# config.Data.outLFNDirBase = '/store/user/mkizilov/CSC_layer_geom/Run2023BC_IdealGeometry'
config.Data.outLFNDirBase = '/store/user/mkizilov/crab3_out/ME11_ana/2023-12-13_v2/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Site.storageSite = 'T3_CH_CERNBOX'