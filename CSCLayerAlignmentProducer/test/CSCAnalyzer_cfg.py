import FWCore.ParameterSet.Config as cms
#from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from Configuration.Eras.Era_Run3_cff import Run3

#process = cms.Process('analyzer',Phase2C9)
process = cms.Process('CSCAnalyzer',Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.MagneticField_0T_cff') #0T for cruzet runs

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('TrackingTools.TrackRefitter.globalMuonTrajectories_cff')
process.load('TrackingTools.TrackFitters.TrackFitters_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
from CondCore.CondDB.CondDB_cfi import *
CondDBSetup = CondDB.clone()
CondDBSetup.__delattr__('connect')


process.CSCGeometryESModule.applyAlignment = cms.bool(True)
## This is the misalignment part
refit = True
debug = False
isCosmic = False
misalign = True
initial_CSC_geometry = '2024-07-22-CSC_2024B_json_03.db'
GPR_file = 'GlobalAlignment_Run2_Run3_v1_ZeroMuonGPR.db'
Global_Tag = '140X_dataRun3_Prompt_v4'
output_root_file = 'out_CSCAnalyzer.root'
dataset_list = 'Run2024B.list'


if misalign:
  db_file = 'sqlite_file:' + initial_CSC_geometry
  process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        connect = cms.string(db_file),
        record = cms.string('CSCAlignmentRcd'),
        tag = cms.string('CSCAlignmentRcd')
    ),
    cms.PSet(
        connect = cms.string(db_file),
        record = cms.string('CSCAlignmentErrorExtendedRcd'),
        tag = cms.string('CSCAlignmentErrorExtendedRcd')
    ),
    cms.PSet(
		connect = cms.string('sqlite_file:'+GPR_file),
		record=cms.string('GlobalPositionRcd'), 
		tag = cms.string('GlobalPositionRcd'))
  )
  process.CSCGeometryESModule.applyAlignment = cms.bool(False)

### Take GPR from file
# from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup

CondDBSetup = CondDB.clone()
CondDBSetup.__delattr__('connect')



process.globalPosition = cms.ESSource("PoolDBESSource", CondDBSetup,
                                      connect = cms.string('sqlite_file:'+GPR_file),
                                      toGet   = cms.VPSet(cms.PSet(record = cms.string("GlobalPositionRcd"), tag = cms.string("GlobalPositionRcd")))
                                      )
process.es_prefer = cms.ESPrefer("PoolDBESSource", "globalPosition")
################################


process.GlobalTag = GlobalTag(process.GlobalTag, Global_Tag, '')




process.MessageLogger.cerr.FwkReport.reportEvery = 5000

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('nEvents',
			-1, #Max number of events 
			VarParsing.multiplicity.singleton, 
			VarParsing.varType.int, 
			"Number of events")
options.parseArguments()

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(options.nEvents)
)
process.maxEvents.input = cms.untracked.int32(-1)


process.source = cms.Source("PoolSource", 
				fileNames = cms.untracked.vstring(options.inputFiles), 
				inputCommands = cms.untracked.vstring(
			"keep *", 
			"drop TotemTimingDigiedmDetSetVector_totemTimingRawToDigi_TotemTiming_reRECO", 
			"drop TotemTimingRecHitedmDetSetVector_totemTimingRecHits__reRECO"
			)
				)

if (len(dataset_list) > 0):
    with open(dataset_list, "r") as f:
        testfiles = f.readlines()
    for testfile in testfiles:
        process.source.fileNames.append('file:/eos/cms' + testfile.strip())


# process.options = cms.untracked.PSet(
#                         SkipEvent = cms.untracked.vstring('ProductNotFound')
#                         )

process.TFileService = cms.Service("TFileService", fileName = cms.string(output_root_file)) #variable name set above

from TrackingTools.TrackRefitter.globalMuonTrajectories_cff import *
process.MuonAlignmentFromReferenceGlobalMuonRefit = globalMuons.clone()
process.MuonAlignmentFromReferenceGlobalMuonRefit.Tracks = cms.InputTag("ALCARECOMuAlCalIsolatedMu:TrackerOnly")
process.MuonAlignmentFromReferenceGlobalMuonRefit.TrackTransformer.RefitRPCHits = cms.bool(False)

process.CSCAnalyzer = cms.EDAnalyzer('CSCAnalyzer', 
	process.MuonServiceProxy,
    csc2DRecHits = cms.InputTag("csc2DRecHits"), 
	gemRecHits = cms.InputTag("gemRecHits"), 
	gemSimHits = cms.InputTag("g4SimHits", "MuonGEMHits"), 
    muons = cms.InputTag("ALCARECOMuAlCalIsolatedMu:SelectedMuons"),
	ref_track = cms.InputTag("MuonAlignmentFromReferenceGlobalMuonRefit:Refitted"),
	vertexCollection = cms.InputTag("offlinePrimaryVerticies"),
    debug = cms.bool(debug),
    isCosmic = cms.bool(isCosmic),
    refitter = cms.bool(refit)
)

#process.p = cms.EndPath(process.analyzer)
process.p = cms.Path(process.MuonAlignmentFromReferenceGlobalMuonRefit + process.CSCAnalyzer)
