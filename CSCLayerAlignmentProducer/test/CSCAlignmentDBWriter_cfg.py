import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3
from Configuration.AlCa.GlobalTag import GlobalTag

# Input files
initial_CSC_geometry = '2023-10-07_2022D_03.db'
GPR_file = 'GlobalAlignment_Run2_Run3_v1_ZeroMuonGPR.db'
CSCLayerAlignmentFile = '../script/standAloneGemAlignment/CSC_layer_al_2023-11-08-v1_small_2DOF.csv'
output_DB_file = 'CSC_layer_al_2023-11-08-v1_2DOF.db'
Global_Tag = '131X_mcRun3_2023_design_v6'
#


process = cms.Process("TEST", Run3)
# Message logger service
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.GlobalTag = GlobalTag(process.GlobalTag, Global_Tag, '')
process.source = cms.Source("EmptySource")

import Geometry.DTGeometryBuilder.dtGeometryDB_cfi
process.DTGeometryMuonMisalignedProducer = Geometry.DTGeometryBuilder.dtGeometryDB_cfi.DTGeometryESModule.clone()
process.DTGeometryMuonMisalignedProducer.appendToDataLabel = 'idealForMuonMisalignedProducer'
process.DTGeometryMuonMisalignedProducer.applyAlignment = cms.bool(False)
import Geometry.CSCGeometryBuilder.cscGeometryDB_cfi
process.CSCGeometryMuonMisalignedProducer = Geometry.CSCGeometryBuilder.cscGeometryDB_cfi.CSCGeometryESModule.clone()
process.CSCGeometryMuonMisalignedProducer.appendToDataLabel = 'idealForMuonMisalignedProducer'
process.CSCGeometryMuonMisalignedProducer.applyAlignment = cms.bool(True)
process.GEMGeometryMuonMisalignedProducer = Geometry.GEMGeometryBuilder.gemGeometryDB_cfi.GEMGeometryESModule.clone()
process.GEMGeometryMuonMisalignedProducer.appendToDataLabel = 'idealForMuonMisalignedProducer'
process.GEMGeometryMuonMisalignedProducer.applyAlignment = cms.bool(False)

process.CSCAlignmentDBWriter = cms.EDAnalyzer("CSCAlignmentDBWriter",
                                       CSCLayerAlignmentFile = cms.untracked.string(CSCLayerAlignmentFile),
                                      )

# Database output service if you want to store soemthing in MisalignedMuon
from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup

process.muonCscAlignment = cms.ESSource("PoolDBESSource", CondDBSetup,
                                    connect = cms.string('sqlite_file:{initial_CSC_geometry}'.format(initial_CSC_geometry = initial_CSC_geometry)),
                                    toGet   = cms.VPSet(cms.PSet(record = cms.string("CSCAlignmentRcd"), tag = cms.string("CSCAlignmentRcd")))
                                    )
process.es_prefer_muonCscAlignment = cms.ESPrefer("PoolDBESSource","muonCscAlignment")

process.globalPosition = cms.ESSource("PoolDBESSource", CondDBSetup,
                                    connect = cms.string('sqlite_file:{GPR_file}'.format(GPR_file = GPR_file)),
                                    toGet   = cms.VPSet(cms.PSet(record = cms.string("GlobalPositionRcd"), tag = cms.string("GlobalPositionRcd")))
                                    )
process.es_prefer_globalPosition = cms.ESPrefer("PoolDBESSource","globalPosition")


process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    CondDBSetup,
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('DTAlignmentRcd'),
        tag = cms.string('DTAlignmentRcd')
        ),
        cms.PSet(
            record = cms.string('DTAlignmentErrorExtendedRcd'),
            tag = cms.string('DTAlignmentErrorExtendedRcd')
        ),
        cms.PSet(
            record = cms.string('CSCAlignmentRcd'),
            tag = cms.string('CSCAlignmentRcd')
        ),
        cms.PSet(
            record = cms.string('CSCAlignmentErrorExtendedRcd'),
            tag = cms.string('CSCAlignmentErrorExtendedRcd')
        ),
        cms.PSet(
            record = cms.string('GEMAlignmentRcd'),
            tag = cms.string('GEMAlignmentRcd')
        ),
        cms.PSet(
            record = cms.string('GEMAlignmentErrorExtendedRcd'),
            tag = cms.string('GEMAlignmentErrorExtendedRcd')
        ),
        cms.PSet(
            record = cms.string('GlobalPositionRcd'),
            tag = cms.string('GlobalPositionRcd')
        )
    ),

    connect = cms.string('sqlite_file:{output_DB_file}'.format(output_DB_file = output_DB_file))
)
process.p1 = cms.Path(process.CSCAlignmentDBWriter)
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    )
)
